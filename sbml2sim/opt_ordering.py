import s2s
import proteomic as ptc
from bioutils import Proteomics, convert_ibaq_to_concentrations

import json
import argparse
import os
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

def main():
    args = parse_args()

    sbml = s2s.SBMLDoc(args.input_file) # must be simulable, see simulable.sh
    proteomics: Proteomics = None
    with open(args.proteomics) as f:
        proteomics = json.load(f)
        
    concentrations = convert_ibaq_to_concentrations(sbml, proteomics, args.tissue)
    
    sbml.input_start_random_concentration() # TODO: start random concentration every restart
    capacity = 40
    
    workers = int(args.workers)
    parallel_simulator = s2s.ParallelSimulator(workers, capacity)
    
    for (species, conc) in concentrations.items():
        if sbml.is_input(species):
            sbml.set_initial_concentration(species, conc)
        elif not sbml.is_output(species):
            parallel_simulator.add_real_concentration(species, conc)
            
    parallel_simulator.order_real_concentration()
    
    for _ in range(capacity):
        parallel_simulator.add_worker(s2s.rr_simualtor(sbml))
    total_start_time = time.time()
    for _ in range(10):
        start_time = time.time()
        results = parallel_simulator.simulate()
        end_time = time.time()
        elapsed_time = end_time - start_time
        print(f"[INFO] results: {results}")
        print(f"[INFO] time: {elapsed_time:.2f}")
    elapsed_time = time.time() - total_start_time
    print(f"[INFO] total time: {elapsed_time}")
    
if __name__ == "__main__":
    main()