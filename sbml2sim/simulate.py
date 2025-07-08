import math
import numpy as np
import roadrunner as rr

from s2s import SBMLDoc
from bioutils import SpeciesId
import pandas as pd
import argparse
import sys


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
    # ritorna 10**30 moltiplicato per il numero di specie che hanno concentrazione negativa
    names = [name for name in sim.colnames if not name.startswith("avg_") and name not in ("time", "get_time")]
    # Trova gli indici delle specie rilevanti
    indices = [i for i, name in enumerate(sim.colnames) if name in names]
    # Prendi i valori all'ultimo tempo simulato
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

def run(sbml: SBMLDoc, plot=False,only_avg=False, output_file="simulation") -> None:
    file_sbml = sbml.convert_to_string()
    r = rr.RoadRunner(file_sbml)
    output_file_csv = output_file+".csv"
    output_file_png = output_file+".png"
    sim = r.simulate(0, 100, output_file=output_file_csv)
    if plot:
        import matplotlib.pyplot as plt

        df = pd.read_csv(output_file_csv)
        cols = [col for col in df.columns if col.find("time") == -1 and (not only_avg or col.startswith('avg_')) and (only_avg or not col.startswith('avg_')) and not sbml.is_output(col.replace("avg_","",1))]
        df[cols].plot()
        plt.xlabel('Time')
        plt.ylabel('Value')
        plt.title('Simulation Results')
        plt.savefig(output_file_png)
        
def parse_args():
    parser = argparse.ArgumentParser(description="Simulate SBML model and plot results.")
    parser.add_argument("input_file", help="Path to the SBML input file")
    parser.add_argument("--plot", action="store_true", help="Plot the simulation results")
    parser.add_argument("--only-avg", action="store_true", help="Plot only avg_* columns")
    parser.add_argument("--output-file", default="simulation", help="Base name for output files (csv/png)")
    return parser.parse_args()

if __name__ == "__main__":
    args = parse_args()
    sbml = SBMLDoc(args.input_file)
    run(sbml, plot=args.plot, only_avg=args.only_avg, output_file=args.output_file)