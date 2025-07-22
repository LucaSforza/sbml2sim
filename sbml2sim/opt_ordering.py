import math
import s2s
import proteomic as ptc
from bioutils import Proteomics, convert_ibaq_to_concentrations, ParameterId, SpeciesId
import nevergrad as ng


import json
import argparse
import time
import sys
import random

def parameters_to_optimize(args) -> ng.p.Parameter:
    if args.scalar:
       return ng.p.Scalar(lower=-6,upper=6) 
    elif args.unify:
        pass # TODO
    else:
        return ng.p.Choice([10**(-3), 10**3])

def optimize(sbml: s2s.SBMLDoc, args, concentrations: dict[SpeciesId, float], seed: int) -> dict[ParameterId, float]:
    workers = int(args.workers)
    budget = int(args.budget)
    parallel_degree = int(args.parallel_degree)
    parallel_simulator = s2s.ParallelSimulator(workers)

    error_handler = s2s.ordering_error_create()
    
    for (species, conc) in concentrations.items():
        if sbml.is_input(species):
            sbml.set_initial_concentration(species, conc)
        elif not sbml.is_output(species):
            sbml.set_initial_concentration(species, conc)
            error_handler.add_real_concentration(species, conc)
    # TODO: choose error
    error_handler.order_real_concentration()
    
    kinetic_constants: list[ParameterId] = sbml.get_kinetic_constants()

    # parametri da ottimizzare    
    kinetic_param_dict = {}
    for kc in kinetic_constants:
        kinetic_param_dict[kc] = parameters_to_optimize(args)
        
        
    parametrization = ng.p.Dict(**kinetic_param_dict)
    optimizer = ng.optimizers.NGOpt(parametrization=parametrization, budget=budget, num_workers=parallel_degree)
    
    
    for _ in range(parallel_degree):
        sim = s2s.rr_simualtor(sbml)
        for (species, _) in concentrations.items():
            s2s.rr_simulator_set_known_species(sim, species)
        parallel_simulator.add_worker(sim)
    simulators = parallel_simulator.get_simulators()
    
    amm = [] # ammissibli
    
    print("[INFO] start optimitation")
    
    attempts = optimizer.budget // len(simulators)
    
    for i in range(attempts):
        total_start_time = time.time()
        parameters = []
        for sim in simulators:
            parameter = optimizer.ask()
            parameters.append(parameter)
            for (id,value) in parameter.items():
                if not args.scalar:
                    sim.set_parameter(id, value.value)
                else:
                    sim.set_parameter(id,10**value.value)
        start_time = time.time()
        sols = parallel_simulator.simulate(error_handler)
        elapsed_time = time.time() - start_time
        errors = 0
        not_errors = 0
        best_result = math.inf
        for (parameter, (value, error)) in zip(parameters, sols):
            if math.isnan(value):
                print(f"[FATAL ERROR] value is NaN, ignoring this result")
                exit(1)
            if value < best_result:
                best_result = value
            if not error:
                not_errors += 1
                if value < 1e-12:
                    print(f"[INFO] ammissibile: {value}")
                    amm.append({k: v.value for k, v in parameter.items()})
                optimizer.tell(parameter, value)
            else:
                errors += 1
                optimizer.tell(parameter, value)
        total_elapsed_time = time.time() - total_start_time
        print(f"[INFO] attempt {i}/{attempts}, best value: {best_result}, time: {elapsed_time}, total time: {total_elapsed_time}, errors: {errors}/{parallel_degree} not_errors: {not_errors}")

    return amm


def parse_args():
    parser = argparse.ArgumentParser(description="Optimize kinetic constants basing on the ordering")
    parser.add_argument("input_file", help="Path to the SBML input file")
    parser.add_argument("--proteomics", default="proteomics.json", help="path to proteomics")
    parser.add_argument("--workers",default="4", help="path to proteomics")
    parser.add_argument("--tissue",default="breast_cancer_cell", help="path to proteomics")
    parser.add_argument("--plot", action="store_true", help="Plot the simulation for each kinetic constants found")
    parser.add_argument("--output-file", default="kinetic_constants.json", help="Output file for the list of kinetic constants that satisfy the constraints")
    parser.add_argument("--budget", default="10000", help="Budget for the optimizator")
    parser.add_argument("--parallel-degree", default="40", help="Paralle degree of the optimizator")
    parser.add_argument("--scalar", action="store_true", help="Set this param if the domain is scalar")
    return parser.parse_args()


def main():
    seed = int(time.time())
    s2s.set_seed(seed)
    args = parse_args()

    sbml = s2s.SBMLDoc(args.input_file) # must be simulable, see simulable.sh
    proteomics: Proteomics = None
    with open(args.proteomics) as f:
        proteomics = json.load(f)
        
    concentrations = convert_ibaq_to_concentrations(sbml, proteomics, args.tissue)
    print("[INFO] random start concentrations")
    sys.stdout.flush()
    sbml.random_start_concentration() # TODO: start random concentration every restart
    sol = optimize(sbml, args, concentrations, seed)
    
    with open(args.output_file, "w") as outfile:
        json.dump(sol, outfile, indent=2)
    
if __name__ == "__main__":
    main()