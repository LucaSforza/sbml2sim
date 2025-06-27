import json
import sys
import math
import secrets
import time
from typing import Any

import nevergrad as ng

from bioutils import convert_ibaq_to_concentrations, get_proteomics, assign_concentrations, set_compartement_size, SpeciesId, ParameterId, UniprodId
import s2s as s2s
from proteomic import proteomic
import proteomic as ptc
import simulate as sim

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

def init_model(sbml: s2s.SBMLDoc, file_path: str, tissue: str, concentrations: dict[SpeciesId, float]) -> s2s.SBMLDoc:
    random_seed = secrets.randbits(64)
    s2s.set_seed(random_seed)
    sbml.add_kinetic_laws_if_not_exists()
    sbml.save_converted_file(file_path.replace(".","-modified."))
    new_sbml = sbml.replicate_model_per_tissue([tissue]) # TODO: Technical debt
    new_sbml.random_kinetic_costant_value()
    new_sbml.random_start_concentration()
    new_sbml.add_time_to_model()
    new_sbml.add_avg_calculations_for_all_species()
    new_sbml.assigment_rule_for_inputs()
    
    new_sbml.save_converted_file(file_path.replace(".","-tissues-modified."))
    assign_concentrations(new_sbml, tissue, concentrations)
    return new_sbml
   
def assign_parameters(
    sbml: s2s.SBMLDoc,
    parameters: dict[ParameterId, float]
):
    for (param_id, value) in parameters.items():
        sbml.set_parameter(param_id, 10**value)
 
def set_sbml_for_attempt(
    sbml: s2s.SBMLDoc,
    tissue: str,
    kinetic_constants: dict[ParameterId, float],
    output_constants: dict[ParameterId, float],
    concentrations: dict[SpeciesId, float],
):
    # genera casualmente le concentrazioni e riposiziona le concentrazioni
    
    sbml.random_start_concentration()
    assign_concentrations(sbml,tissue, concentrations)
    # TODO: reset initial value for avg non-constant parameters
    assign_parameters(sbml, kinetic_constants)
    assign_parameters(sbml, output_constants)

# Global variable to count attempts
attempts = 0
best_results = math.inf

def utility_function(
    ng_params: dict[str, dict[ParameterId, float]],
    sbml: s2s.SBMLDoc,
    concentrations: dict[SpeciesId, float],
    tissue_name: str
) -> float:
    """
    Funzione di utilità che prende in input:
    - ng_params: dizionario dei parametri che Nevergrad deve ottimizzare (str -> float)
    - hidden_params: dizionario dei parametri che Nevergrad non vede (fissi o di contesto)
    """
    global attempts, best_results
    start_time = time.time()
    set_sbml_for_attempt(sbml, tissue_name, ng_params["kinetic_constants"], ng_params["output_constants"], concentrations)
    attempts += 1
    print(f"[INFO] attempt: {attempts}")
    # print(f"[INFO] params:\n {ng_params}")
    # Le chiavi di ng_params sono sempre 'f', 'k_1', 'k_2'
    result = 0.0
    integration_errors = 0
    for f in [10**-1]: # [10**(-i) for i in range(1,6+1)]:
        for k_1 in [0]: # [k for k in range(0,7+1)]:
            for k_2 in [0]: # [k for k in range(0,7+1)]:
                sim_result = sim.simulate(sbml, f, k_1, k_2)
                if math.isnan(sim_result):
                    print("[FATAL ERORR] utility function returned NaN")
                    exit(1)
                if math.isinf(sim_result):
                    if sim_result < 0:
                        print("[FATAL ERROR] Negative infinity returned by simulate function. An unexpected error occurred.")
                        exit(1)
                    else:
                        end_time = time.time()
                        
                        integration_errors += 1
                result += sim_result
    end_time = time.time()
    print(f"[INFO] Simulation {attempts} ended with errors: {integration_errors}/{1} in {end_time - start_time:.2f} seconds")
    utility = result + integration_errors*10**9
    print(f"[INFO] utility: {utility}")
    if utility < best_results:
        best_results = utility
    return utility


TISSUE="breast_cancer_cell" # TODO: param

def main():
    file_path = parse_args()

    sbml: s2s.SBMLDoc = s2s.SBMLDoc(file_path)
    # ref: https://bionumbers.hms.harvard.edu/bionumber.aspx?id=115154&ver=1&trm=cell+size+breast+cancer+cell+human+&org=
    volume_cell_breast_cancer_cell = 1.76 * 10**12  # nanometers
    set_compartement_size(sbml, volume_cell_breast_cancer_cell)
    proteins: dict[SpeciesId,UniprodId] = sbml.get_proteins_data()
    
    # TODO: remove all_tissue_names
    (proteomics, _) = get_proteomics(proteins) # APIs

    concentrations: dict[SpeciesId, float] = convert_ibaq_to_concentrations(sbml, proteomics, TISSUE)
    print(concentrations)
    # replica il modello per il tessuto del cancro al seno
    print("[INFO] init model")
    sbml = init_model(sbml, file_path, TISSUE, concentrations)
   
    print("[INFO] get kinetic constants") 
    kinetic_constants: list[ParameterId] = sbml.get_kinetic_constants()
    output_constants: list[ParameterId] = sbml.get_output_constants()
    
    # Build parametrization with separate dictionaries
    kinetic_param_dict = {}
    for kc in kinetic_constants:
        # Parameter: kinetic constant
        kinetic_param_dict[kc] = ng.p.Scalar(lower=-6.0, upper=6.0)
    output_param_dict = {}
    for oc in output_constants:
        # Parameter: output constant
        output_param_dict[oc] = 0 # ng.p.Scalar(lower=-6.0, upper=1.0)
    param_dict = {
        "kinetic_constants": ng.p.Dict(**kinetic_param_dict),
        "output_constants": ng.p.Dict(**output_param_dict)
    }
    parametrization = ng.p.Dict(**param_dict)
    optimizer = ng.optimizers.CMA(parametrization=parametrization, budget=80_000)

    def ng_objective(ng_params):
        # hidden parameters sbml and concentrations
         return utility_function(ng_params, sbml, concentrations, TISSUE)
    print("[INFO] Start Opt")
    start_opt = time.time()
    recommendation = optimizer.minimize(ng_objective)
    end_opt = time.time()
    print("[INFO] Best parameters found:\n", recommendation.value)
    print("[INFO] Best utility found:", best_results)
    print(f"[INFO] time:{end_opt - start_opt:.2f}")
    with open("parameters.json", "w") as f:
        json.dump(recommendation.value, f)
    result: dict[str, dict[ParameterId, float]] = recommendation.value
    set_sbml_for_attempt(sbml, TISSUE, result["kinetic_constants"], result["output_constants"], concentrations)
    sbml.save_converted_file(file_path.replace(".", "-real-tissues-modified."))
    sim.simulate(sbml, 10**(-1), 0, 0, plot=True)