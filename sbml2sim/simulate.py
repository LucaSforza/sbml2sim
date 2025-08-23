import math
from typing import Any
import numpy as np
import roadrunner as rr

from s2s import SBMLDoc, set_seed
from bioutils import SpeciesId, ParameterId, convert_ibaq_to_concentrations
import pandas as pd
import argparse
import json 
import time

SimulationResult=Any


def steady_state_residual(sbml, r: rr.RoadRunner, sim) -> float:
    # Prendi tutti i parametri che iniziano per 'avg_'
    t1_time = 0.8 * sim[-1, 0]  # tempo all'80%
    t2_time = sim[-1, 0]        # tempo finale

    # Trova l'indice del tempo più vicino a t1 e t2
    t_idx_1 = np.argmin(np.abs(sim[:, 0] - t1_time))
    t_idx_2 = len(sim) - 1  # ultimo punto

    avg_indices = [i for i, name in enumerate(sim.colnames)
                   if name.startswith("avg_")]

    if not avg_indices:
        print("[FATAL ERROR] not found any avg parameter")
        exit(1)

    # Calcola le differenze assolute tra t2 e t1
    diffs = [abs(sim[t_idx_2, i] - sim[t_idx_1, i]) for i in avg_indices]
    result = np.linalg.norm(diffs, ord=1)
    return result

def check_for_inconsistency(sim) -> float:
    names = [name for name in sim.colnames if not name.startswith("avg_") and name not in ("time", "get_time")]
    indices = [i for i, name in enumerate(sim.colnames) if name in names]
    values = sim[-1, indices]
    penalty_value = np.linalg.norm(values[values < -1e-1], ord=1)
    return penalty_value

def penalty(sbml, r: rr.RoadRunner,concentration: dict[SpeciesId, float],tissue: str, plot= False, output_file_name="simulation"):
    output_file_csv = None
    output_file_png = None
    if plot:
        output_file_csv = output_file_name+".csv"
        output_file_png = output_file_name+".png"
    try:
        # set_params_to_model(r, params)
        sim = r.simulate(0, 100, output_file=output_file_csv)  # o più lungo
        if plot:
            import matplotlib.pyplot as plt

            df = pd.read_csv(output_file_csv)
            avg_cols = [col for col in df.columns if col.startswith('avg_')]
            df[avg_cols].plot()
            plt.xlabel('Time')
            plt.ylabel('Value')
            plt.title('Simulation Results (avg_*)')
            plt.savefig(output_file_png)
    except Exception as e:
        # print(f"[WARNING] RoadRunner failed (integration error is normal): {e}")
        return (math.inf, math.inf)  # grande penalità se il modello va in crash (errori numerici)
    if not plot:
        penalty = steady_state_residual(sbml, r, sim)
        bio_pen = 0.0
        for specie_id, target_value in concentration.items():
            sim_id = f"{tissue}_{specie_id}"
            if sim_id in sim.colnames:
                idx = sim.colnames.index(sim_id)
                final_value = sim[-1, idx]
                # errore quadratico della quantità
                bio_pen += (final_value - target_value) ** 2
                # Somma anche l'errore quadratico dell'ordine di grandezza
                if final_value > 0 and target_value > 0:
                    log_err = (np.log10(final_value) - np.log10(target_value)) ** 2
                    bio_pen += log_err
                else:
                    bio_pen += 10**3
            # else:
            #     print(f"[FATAL ERROR] species {sim_id} not found in simulation result")
            #     exit(1)
        # penalty += check_for_inconsistency(sim)
    else:
        penalty = 0.9
        bio_pen = 0.0
    return (penalty, bio_pen)


# TODO: da rifare completamente
# returns the penalty
def simulate(sbml: SBMLDoc,concentration: dict[SpeciesId, float],tissue: str,  plot=False, output_file_name="simulation") -> tuple[float,float]:
    file_sbml = sbml.convert_to_string()
    r = rr.RoadRunner(file_sbml)
    return penalty(sbml, r,concentration, tissue, plot, output_file_name)

# returns the result of the simulation
# TODO: plot the real concentrations
def run(sbml: SBMLDoc, plot=False,only_avg=False, output_file="simulation", concentrations: dict[SpeciesId, float]=None) -> SimulationResult:
    file_sbml = sbml.convert_to_string()
    r = rr.RoadRunner(file_sbml)
    output_file_csv = output_file+".csv"
    output_file_png = output_file+".png"
    sim = r.simulate(0, 100, output_file=output_file_csv)
    if plot:
        import matplotlib.pyplot as plt
        if concentrations:
            df = pd.read_csv(output_file_csv)
            cols = [
                col for col in df.columns if col.find("time") == -1 
                and concentrations.get(col.strip("[]")) is not None and (not only_avg or col.startswith('avg_')) 
                and (only_avg or not col.startswith('avg_')) and not sbml.is_output(col.replace("avg_","",1))
            ]
            for (i,species) in enumerate(cols, start=1):
                print("[INFO] plot species: ", species)
                real_conc = concentrations[species.strip("[]")]
                sim_conc = float(df[species].iloc[-1])
                plt.plot(df["time"],df[species], label=species, linewidth=2.5)
                plt.axhline(y=real_conc, color='r', linestyle='--', label=f"{species} (real value)")
                plt.title(f"{species} concentration over time")
                plt.xlabel("Time")
                plt.ylabel("mol/L")
                plt.yscale('log')
                min_conc = real_conc if real_conc < sim_conc else sim_conc
                max_conc = real_conc if real_conc > sim_conc else sim_conc
                plt.ylim(bottom=10**(math.log(min_conc, 10) - 5), top=10**(math.log(max_conc, 10)+ 5))
                plt.legend([f"{species} mol/L"])
                plt.legend()
                plt.savefig(output_file+"_"+str(i)+".png")
                plt.clf()
            
        else:
            df = pd.read_csv(output_file_csv)
            cols = [col for col in df.columns if col.find("time") == -1 
                    and (not only_avg or col.startswith('avg_')) and (only_avg or not col.startswith('avg_')) 
                    and not sbml.is_output(col.replace("avg_","",1))
                ]
            df[cols].plot()
            plt.xlabel('Secondi')
            plt.ylabel('Concentrazione mol/L')
            plt.title('Simulazione')
            # plt.yscale('log')
            plt.savefig(output_file_png)
    return sim
        
        
def parse_args():
    parser = argparse.ArgumentParser(description="Simulate SBML model and plot results.")
    parser.add_argument("input_file", help="Path to the SBML input file")
    parser.add_argument("--plot", default=True, action="store_true", help="Plot the simulation results")
    parser.add_argument("--only-avg", action="store_true", help="Plot only avg_* columns")
    parser.add_argument("--output-file", default="simulation", help="Base name for output files (csv/png)")
    parser.add_argument("--kinetic-constants", default="ord_params.json", help="file path to kientic constants")
    parser.add_argument("--proteomics", default=None, help="file path to proteomics")
    parser.add_argument("--tissue", default="breast_cancer_cell", help="tissue name")
    parser.add_argument("--plot-constrained", action="store_true", help="plot only constrained species")
    parser.add_argument("--log-parameters", action="store_true", help="plot only constrained species")
    parser.add_argument("--random", action="store_true", help="random start concentrations")
    return parser.parse_args()

def main2():
    args = parse_args()
    sbml = SBMLDoc(args.input_file)
    if args.random:
        sbml.random_start_concentration()
    if args.kinetic_constants is not None:
        kinetic_constants: dict[str, dict[ParameterId, float]] = None
        with open(args.kinetic_constants) as f:
            kinetic_constants = json.load(f)
        k_constants = kinetic_constants["kinetic_constants"]
        output_constants = kinetic_constants["output_constants"]
        input_constants = kinetic_constants["input_constants"]
        
        for (param, value) in k_constants.items():
            if args.log_parameters:
                sbml.set_parameter(param, 10**value)
            else:
                sbml.set_parameter(param, value)
            
        for (param, value) in output_constants.items():
            sbml.set_parameter(param, 10**value)
            
        for (species, value) in input_constants.items():
            sbml.set_initial_concentration(species, 10**value)
    if args.proteomics is not None:
        proteomics = None
        with open(args.proteomics) as f:
            proteomics = json.load(f)
        concentrations = convert_ibaq_to_concentrations(sbml, proteomics, args.tissue)
        for (species, conc) in concentrations.items():
            if sbml.is_input(species) or not sbml.is_output(species):
                sbml.set_initial_concentration(species, conc)
    plot_concentrations = None
    if args.plot_constrained:
        print("[INFO] plot concentrations")
        plot_concentrations = concentrations
        
    run(sbml, plot=args.plot, only_avg=args.only_avg, output_file=args.output_file, concentrations=plot_concentrations)

def main():
    args = parse_args()
    seed = int(time.time())
    set_seed(seed)
    sbml = SBMLDoc(args.input_file)
    if args.random:
        sbml.random_start_concentration()
    if args.kinetic_constants is not None:
        kinetic_constants: list[dict[ParameterId, float]] = None
        with open(args.kinetic_constants) as f:
            kinetic_constants = json.load(f)
        if not isinstance(kinetic_constants, list):
            kinetic_constants = [kinetic_constants]
        # TODO: choose the element in the list
        print(kinetic_constants)
        for param_id, value in kinetic_constants[0].items():
            if args.log_parameters:
                sbml.set_parameter(param_id, 10**value)
            else:
                sbml.set_parameter(param_id, value)
    if args.proteomics is not None:
        proteomics = None
        with open(args.proteomics) as f:
            proteomics = json.load(f)
        concentrations = convert_ibaq_to_concentrations(sbml, proteomics, args.tissue)
        for (species, conc) in concentrations.items():
            if sbml.is_input(species) or not sbml.is_output(species):
                sbml.set_initial_concentration(species, conc)
    plot_concentrations = None
    if args.plot_constrained:
        print("[INFO] plot concentrations")
        plot_concentrations = concentrations
        
    run(sbml, plot=args.plot, only_avg=args.only_avg, output_file=args.output_file, concentrations=plot_concentrations)

if __name__ == "__main__":
    main2()