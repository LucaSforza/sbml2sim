import math
import s2s
import proteomic as ptc
from bioutils import Proteomics, convert_ibaq_to_concentrations, ParameterId, SpeciesId
import nevergrad as ng

import json
import argparse
import time

def parse_args():
    parser = argparse.ArgumentParser(description="Optimize kinetic constants basing on the ordering")
    parser.add_argument("input_file", help="Path to the SBML input file")
    parser.add_argument("--proteomics", default="proteomics.json", help="path to proteomics")
    parser.add_argument("--workers",default=4, help="path to proteomics")
    parser.add_argument("--tissue",default="breast_cancer_cell", help="path to proteomics")
    parser.add_argument("--plot", action="store_true", help="Plot the simulation for each kinetic constants found")
    parser.add_argument("--output-file", default="kinetic_constants.json", help="Output file for the list of kinetic constants that satisfy the constraints")
    return parser.parse_args()


def optimize(sbml: s2s.SBMLDoc, args, concentrations: dict[SpeciesId, float]) -> dict[ParameterId, float]:
    capacity = 40
    
    workers = int(args.workers)
    parallel_simulator = s2s.ParallelSimulator(workers)
    
    for (species, conc) in concentrations.items():
        if sbml.is_input(species):
            sbml.set_initial_concentration(species, conc)
        elif not sbml.is_output(species):
            parallel_simulator.add_real_concentration(species, conc)
    
    # TODO: choose error
    parallel_simulator.order_real_concentration()
    
    kinetic_constants: list[ParameterId] = sbml.get_kinetic_constants()
    
    # Build parametrization with separate dictionaries
    # TODO: cerca quelli da migliorare
    kinetic_param_dict = {}
    for kc in kinetic_constants:
        # Parameter: kinetic constant (boolean parameter example)
        kinetic_param_dict[kc] = ng.p.Choice([1e-3, 1, 1e3])
        
    parametrization = ng.p.Dict(**kinetic_param_dict)
    optimizer = ng.optimizers.NGOpt21(parametrization=parametrization, budget=20_000, num_workers=capacity)
    
    for _ in range(capacity):
        parallel_simulator.add_worker(s2s.rr_simualtor(sbml))
    
    simulators = parallel_simulator.get_simulators()
    
    amm = [] # ammissibli
    
    print("[INFO] start optimitation")
    
    attempts = optimizer.budget // len(simulators)

    
    for i in range(attempts):
        start_time = time.time()
        parameters = []
        for sim in simulators:
            parameter = optimizer.ask()
            parameters.append(parameter)
            for (id,value) in parameter.items():
                sim.set_parameter(id, value.value)
        sols = parallel_simulator.simulate()
        errors = 0
        not_errors = 0
        for (parameter, (value, error)) in zip(parameters, sols):
            if not error:
                not_errors += 1
                if math.isnan(value):
                    print(f"[WARNING] value is NaN, ignoring this result")
                    optimizer.tell(parameter, 10**4)
                    continue
                if value < 1e-12:
                    print(f"[INFO] ammissibile: {value}")
                    amm.append({k: v.value for k, v in parameter.items()})
                optimizer.tell(parameter, value)
            else:
                errors += 1
                optimizer.tell(parameter, 10**4)
        elapsed_time = time.time() - start_time
        print(f"[INFO] attempt {i}/{attempts}, time: {elapsed_time}, errors: {errors}/{capacity} not_errors: {not_errors}")

    return amm

def main():
    args = parse_args()

    sbml = s2s.SBMLDoc(args.input_file) # must be simulable, see simulable.sh
    proteomics: Proteomics = None
    with open(args.proteomics) as f:
        proteomics = json.load(f)
        
    concentrations = convert_ibaq_to_concentrations(sbml, proteomics, args.tissue)
    
    sbml.random_start_concentration() # TODO: start random concentration every restart
    sol = optimize(sbml, args, concentrations)
    
    print(sol)
    
    with open(args.output_file, "w") as outfile:
        json.dump(sol, outfile, indent=2)
    
if __name__ == "__main__":
    main()