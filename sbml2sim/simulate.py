import math
import numpy as np
import roadrunner as rr

from s2s import SBMLDoc
import pandas as pd


def steady_state_residual(r: rr.RoadRunner) -> float:
    # Prendi tutti i parametri che iniziano per 'avg_'
    rates: rr._roadrunner.NamedArray = r.getRatesOfChangeNamedArray()
    avg_params = [rates[k] for k in rates.colnames if k.startswith('avg_')]
    if not avg_params:
        # nessun parametro avg
        print("[FATAL ERROR] not found any avg parameter")
        exit(1)
    return np.linalg.norm(avg_params, ord=1)

def penalty(r: rr.RoadRunner, plot= False):
    output_file = None
    if plot:
        output_file = "simulation.csv"
    try:
        # set_params_to_model(r, params)
        r.simulate(0, 100, output_file=output_file)  # o più lungo
        if plot:
            import matplotlib.pyplot as plt

            df = pd.read_csv(output_file)
            avg_cols = [col for col in df.columns if col.startswith('avg_')]
            df[avg_cols].plot()
            plt.xlabel('Time')
            plt.ylabel('Value')
            plt.title('Simulation Results (avg_*)')
            plt.legend()
            plt.savefig("simulation.png")
    except Exception as e:
        # print(f"[WARNING] RoadRunner failed (integration error is normal): {e}")
        return math.inf  # grande penalità se il modello va in crash (errori numerici)
    penalty = steady_state_residual(r)
    return penalty


# returns the penalty
def simulate(sbml: SBMLDoc, f: float, k_1: float, k_2: float, plot=False) -> float:
    # TODO: setta f, k_1 e k_2
    sbml.set_parameter("input_constant_f", f)
    sbml.set_parameter("input_constant_k_1", k_1)
    sbml.set_parameter("input_constant_k_2", k_2)
    file_sbml = sbml.convert_to_string()
    r = rr.RoadRunner(file_sbml)
    return penalty(r, plot)