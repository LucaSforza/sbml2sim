import math
from typing import Any
import opt2
import s2s
import proteomic as ptc
from bioutils import Proteomics, convert_ibaq_to_concentrations, ParameterId, SpeciesId
import nevergrad as ng


import json
import argparse
import time
import sys
import random
import itertools

def set_parameters_to_model(sim: s2s.Simulator, args,parameter: ng.p.Parameter):    
    if not args.unify:
        for (id,value) in parameter.items():
            if args.scalar:
                sim.set_parameter(id,10**value.value)
            else:
                sim.set_parameter(id, value.value)
    else:
        data1: dict[str,float] = parameter[args.file1]
        # TODO: sono logaritmici oppure no?
        visited = set()
        for (id,value) in data1.value.items():
            visited.add(id)
            sim.set_parameter(id,value)
        data2: dict[str,float] = parameter[args.file2]
        # TODO: che succede se una costante cinetica viene risettata?
        for (id,value) in data2.value.items():
            if id in visited:
                print("[ERROR] parameter overwritten")
            sim.set_parameter(id,value) 
        if args.file3:
            data3: dict[str,float] = parameter[args.file3]
            for (id,value) in data3.value.items():
                if id in visited:
                    print("[ERROR] parameter overwritten")
                sim.set_parameter(id,value) 
        else:
            for (k, v) in parameter.items():
                if k != args.file1 and k != args.file2:
                    sim.set_parameter(k,v.value)
        

def parameters_to_optimize(args, parameterId: ParameterId) -> ng.p.Parameter:
    if args.scalar:
       return ng.p.Scalar(lower=-6,upper=6) 
    elif args.unify:
        print("[FATAL ERROR] unexpected error")
        exit(1)
    else:
        if parameterId.startswith("k_output"):
            return ng.p.Choice([10**i for i in range(-48,-35)])
        return ng.p.Choice([10**i for i in range(-6, 7)])

def optimize(sbml: s2s.SBMLDoc, args, concentrations: dict[SpeciesId, float], seed: int) -> dict[ParameterId, float]:
    workers = int(args.workers)
    budget = int(args.budget)
    parallel_degree = int(args.parallel_degree)
    parallel_simulator = s2s.ParallelSimulator(workers)

    error_handler = s2s.error_sum_create()
    s2s.error_sum_add_handler(error_handler, s2s.stability_error_create(), 1e6)
    s2s.error_sum_add_handler(error_handler, s2s.transitorial_error_create(), 1)
    """
    for (species, conc) in concentrations.items():
        if sbml.is_input(species):
            sbml.set_initial_concentration(species, conc)
        elif not sbml.is_output(species):
            sbml.set_initial_concentration(species, conc)
            error_handler.add_real_concentration(species, conc)
        else:
            error_handler.add_output(species)
    # TODO: choose error
    error_handler.order_real_concentration()
    """
    
    kinetic_constants: list[ParameterId] = sbml.get_kinetic_constants()

    # parametri da ottimizzare    
    kinetic_param_dict = {}
    for kc in kinetic_constants:
        if not args.unify:
            kinetic_param_dict[kc] = parameters_to_optimize(args, kc)
        else:
            # leggi i file da cui prendere i dati
            # TODO: generalizza per un numero arbitrario di file
            data1: list[dict[ParameterId, float]] = None
            with open(args.file1) as f:
                data1 = json.load(f)
            kinetic_param_dict[args.file1] = ng.p.Choice(data1)
            data2: list[dict[ParameterId, float]] = None
            with open(args.file2) as f:
                data2 = json.load(f)
            kinetic_param_dict[args.file2] = ng.p.Choice(data2)
            if args.file3:
                data3: list[dict[ParameterId, float]] = None
                with open(args.file3) as f:
                    data3 = json.load(f)
                kinetic_param_dict[args.file3] = ng.p.Choice(data3)
            else:
                for kc in kinetic_constants:
                    if kc not in data1[0] and kc not in data2[0]:
                        kinetic_param_dict[kc] = ng.p.Choice([10**(-3), 1, 10**3])
        
        
    parametrization = ng.p.Dict(**kinetic_param_dict)
    optimizer = ng.optimizers.DiscreteDE(parametrization=parametrization, budget=budget, num_workers=parallel_degree)
    
    
    for _ in range(parallel_degree):
        sim = s2s.rr_simualtor(sbml)
        #for (species, _) in concentrations.items():
        #    s2s.rr_simulator_set_known_species(sim, species)
        parallel_simulator.add_worker(sim)
    simulators = parallel_simulator.get_simulators()
    
    amm = [] # ammissibli
    
    print("[INFO] start optimitation")
    
    attempts = optimizer.budget // len(simulators)
    best_result_global = math.inf
    for i in range(attempts):
        total_start_time = time.time()
        parameters = []
        for sim in simulators:
            parameter = optimizer.ask()
            parameters.append(parameter)
            set_parameters_to_model(sim, args,parameter)
                
        start_time = time.time()
        sols = parallel_simulator.simulate(error_handler)
        elapsed_time = time.time() - start_time

        errors = 0
        error_is_best = False
        not_errors = 0
        best_result = math.inf
        for (parameter, (value, error)) in zip(parameters, sols):
            if math.isnan(value):
                print(f"[FATAL ERROR] value is NaN, ignoring this result")
                exit(1)
            if value < best_result:
                if error: error_is_best = True
                else: error_is_best = False
                best_result = value
            if not error:
                not_errors += 1
                if value <= 0:
                    print(f"[INFO] ammissibile: {value}")
                    amm.append({k: v.value for k, v in parameter.items()})
                optimizer.tell(parameter, value)
            else:
                errors += 1
                optimizer.tell(parameter, value)
        if error_is_best:
            print("[WARNING] the best is an error")
        if best_result < best_result_global:
            best_result_global = best_result
        total_elapsed_time = time.time() - total_start_time
        print(f"[INFO] attempt {i+1}/{attempts}, best value: {best_result}, best result global: {best_result_global}, time: {elapsed_time}, total time: {total_elapsed_time}, errors: {errors}/{parallel_degree} not_errors: {not_errors}")
    amm.append(optimizer.provide_recommendation().value)
    return amm

def total_search(sbml: s2s.SBMLDoc, args, concentrations: dict[SpeciesId, float], seed: int):
    k_constants = sbml.get_kinetic_constants()
   
    workers = int(args.workers)
    parallel_degree = int(args.parallel_degree)
    parallel_simulator = s2s.ParallelSimulator(workers)
    
    error_handler = s2s.error_sum_create()
    s2s.error_sum_add_handler(error_handler, s2s.stability_error_create())
    s2s.error_sum_add_handler(error_handler, s2s.transitorial_error_create())
    
    to_search: list[tuple[ParameterId, list[float]]] = []
    combinations = 1
    for k in k_constants:
        choices = None
        if k.startswith("k_input") or k.startswith("k_output"):
            choices = [10**i for i in range(-24,-11)]
        else:
            choices = [10**i for i in range(-6,7)]
        combinations *= len(choices)
        to_search.append((k,choices))
    
    for _ in range(parallel_degree):
        sim = s2s.rr_simualtor(sbml)
        #for (species, _) in concentrations.items():
        #    s2s.rr_simulator_set_known_species(sim, species)
        parallel_simulator.add_worker(sim)
    simulators = parallel_simulator.get_simulators()
    amm = [] # ammissibli
     
    keys = k_constants
    choices_lists = [choices for (_, choices) in to_search]

    params_iter = (dict(zip(keys, combo)) for combo in itertools.product(*choices_lists))

    batch:list[dict[ParameterId, float]] = []
    best_value = math.inf
    best_combinations = None
    count = 0
    start = time.time()
    for i,param in enumerate(params_iter):
        errors = 0
        batch.append(param)
        if len(batch) == parallel_degree:
            local_best_value = math.inf
            local_best_combinations = None
            # assign each parameter set to a simulator and run
            for sim, p in zip(simulators, batch):
                for pid, val in p.items():
                    sim.set_parameter(pid, val)
            sols = parallel_simulator.simulate(error_handler)
            for p, (value, error) in zip(batch, sols):
                if value < local_best_value:
                    local_best_value = value
                    local_best_combinations = p.copy()
                if error:
                    errors += 1
                if not error and value < 1e-10:
                    amm.append(p.copy())
            if local_best_value < best_value:
                best_value = local_best_value
                best_combinations = local_best_combinations
            count += 1
            if count >= 1000:
                count = 0
                print(f"[INFO] {i+1}/{combinations}")
                print(f"[INFO]      best value: {best_value}")
                print(f"[INFO]      local best value: {local_best_value}")
                print(f"[INFO]      time: {time.time() - start}")
                print(f"[INFO]      errors: {errors}/{len(batch)}")
                print(f"[INFO]      len amm: {len(amm)}")
                start = time.time()
            batch = []
            

    # remaining partial batch
    if batch:
        for sim, p in zip(simulators, batch):
            for pid, val in p.items():
                sim.set_parameter(pid, val)
        sols = parallel_simulator.simulate(error_handler)
        for p, (value, error) in zip(batch, sols):
            if not error and value < 1e-12:
                amm.append(p.copy())
                if len(amm) >= 100:
                        return amm
    amm.append(best_combinations)
    return amm
            
    

def parse_args():
    parser = argparse.ArgumentParser(description="Optimize kinetic constants basing on the ordering")
    parser.add_argument("input_file", help="Path to the SBML input file")
    parser.add_argument("--proteomics", default="proteomics.json", help="path to proteomics")
    parser.add_argument("--workers",default="4", help="path to proteomics")
    parser.add_argument("--tissue",default="breast_cancer_cell", help="path to proteomics")
    parser.add_argument("--plot", action="store_true", help="Plot the simulation for each kinetic constants found")
    parser.add_argument("--output-file", default="kinetic_constants.json", help="Output file for the list of kinetic constants that satisfy the constraints")
    parser.add_argument("--budget", default="3000", help="Budget for the optimizator")
    parser.add_argument("--parallel-degree", default="16", help="Paralle degree of the optimizator")
    parser.add_argument("--scalar", action="store_true", help="Set this param if the domain is scalar")
    parser.add_argument("--unify", action="store_true", help="unify 2 set of kinetic constants")
    parser.add_argument("--file1", default=None ,help="file1 for unify")
    parser.add_argument("--file2", default=None, help="file2 to unify")
    parser.add_argument("--file3", default=None, help="file3 to unify (optional)")
    parser.add_argument("--serial", action="store_true", help="run serial")
    parser.add_argument("--optimizer", default=None, help="select algoritm")
    parser.add_argument("--all-combination", action="store_true", help="search for every singole combinations")
    return parser.parse_args()


def main():
    args = parse_args()
    if args.serial:
        opt2.main(args)
        return
    seed = int(time.time())
    # TODO: set seed
    s2s.set_seed(seed)
    

    sbml = s2s.SBMLDoc(args.input_file) # must be simulable, see simulable.sh
    proteomics: Proteomics = None
    with open(args.proteomics) as f:
        proteomics = json.load(f)
        
    concentrations = convert_ibaq_to_concentrations(sbml, proteomics, args.tissue)
    # sbml.random_start_concentration() # TODO: start random concentration every restart
    sols: dict[str, dict[Any]] = None
    if args.all_combination:
        sols = total_search(sbml, args, concentrations, seed)
    else:
        sols = optimize(sbml, args, concentrations, seed)
    
    result: list[dict[ParameterId, float]] = []
    if args.unify:
        # TODO: generalizza per un numero arbitrario di file
        for s in sols:
            s1 = dict()
            s1.update(s[args.file1])
            s1.update(s[args.file2])
            if args.file3:
                s1.update(s[args.file3])
            else:
                for (k, v) in s.items():
                    if k != args.file1 and k != args.file2:
                        s1[k] = v
            result.append(s1)
    else:
        result = sols
    with open(args.output_file, "w") as outfile:
        json.dump(result, outfile, indent=2)
    
if __name__ == "__main__":
    main()