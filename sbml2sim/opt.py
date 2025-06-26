import sys
import time
from typing import Any

import nevergrad as ng
import pandas as pd
import matplotlib.pyplot as plt

from sbml2sim.bioutils import choose_tissue_for_replication, get_proteomics, save_proteomics
import sbml2sim.s2s as s2s
from sbml2sim.proteomic import proteomic
import sbml2sim.proteomic as ptc

PROGRAM_NAME = sys.argv[0]

def usage():
    print(f"{PROGRAM_NAME} <SBML file path>")

# @returns SBML path
def parse_args() -> str:
    if len(sys.argv) <= 1:
        print("[ERROR] must specify the SBML file path")
        usage()
        exit(1)
    return sys.argv[1]

def init_model(sbml: s2s.SBMLDoc, file_path: str, tissue: str, proteomics: Any) -> s2s.SBMLDoc:
    random_seed = int(time.time() * 1000)
    s2s.set_seed(random_seed)
    sbml.add_kinetic_laws_if_not_exists()
    sbml.save_converted_file(file_path.replace(".","-modified."))
    new_sbml = sbml.replicate_model_per_tissue([tissue])
    new_sbml.random_kinetic_costant_value()
    new_sbml.random_start_concentration()
    new_sbml.add_time_to_model()
    new_sbml.add_avg_calculation_for_all_proteins()
    new_sbml.save_converted_file(file_path.replace(".","-tissues-modified."))
    # TODO: setta tutti i valori che conosciamo

def run_model(sbml: s2s.SBMLDoc, f: float, k_1: float, k_2: float) -> Any:
    pass

def utility_function(ng_params: dict[str, float], hidden_params: dict[str, Any]) -> float:
    """
    Funzione di utilità che prende in input:
    - ng_params: dizionario dei parametri che Nevergrad deve ottimizzare (str -> float)
    - hidden_params: dizionario dei parametri che Nevergrad non vede (fissi o di contesto)
    """
    # Le chiavi di ng_params sono sempre 'f', 'k_1', 'k_2'
    f = float(ng_params['f'])
    k_1 = float(ng_params['k_1'])
    k_2 = float(ng_params['k_2'])
    result = run_model(
        hidden_params['sbml'],
        f,
        k_1,
        k_2
    )
    return float(result)

def main():
    file_path = parse_args()
    sim_output = "result.csv"

    sbml: s2s.SBMLDoc = s2s.SBMLDoc(file_path)
    
    proteins: dict[str,str] = sbml.get_proteins_data()
    
    (proteomics, all_tissue_names) = get_proteomics(proteins)
    tissue_names = choose_tissue_for_replication(all_tissue_names, proteomics)
    print(tissue_names)
    save_proteomics(proteomics)

    # Setup per Nevergrad: le chiavi sono sempre 'f', 'k_1', 'k_2'
    
    sbml = init_model(sbml, file_path, "breast_cancer_cell", proteomics)
    
    kinetic_constants: list[str] = sbml.get_kinetic_constants()
    output_constants: list[str] = sbml.get_output_constants()
    
    # Build parametrization with separate dictionaries
    kinetic_param_dict = {}
    for kc in kinetic_constants:
        kinetic_param_dict[kc] = ng.p.Scalar(lower=-6.0, upper=6.0)
    output_param_dict = {}
    for oc in output_constants:
        output_param_dict[oc] = ng.p.Scalar(lower=-4.0, upper=1.0)
    param_dict = {
        "kinetic_constants": ng.p.Dict(**kinetic_param_dict),
        "output_constants": ng.p.Dict(**output_param_dict)
    }
    parametrization = ng.p.Dict(**param_dict)
    optimizer = ng.optimizers.OnePlusOne(parametrization=parametrization, budget=10)
    hidden_params = {'sbml': sbml}

    def ng_objective(ng_params):
        return utility_function(ng_params, hidden_params)

    recommendation = optimizer.minimize(ng_objective)
    print("Best parameters found:", recommendation.value)
