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
    parameters: dict[ParameterId, float] | None
):
    if parameters is None: return
    for (param_id, value) in parameters.items():
        sbml.set_parameter(param_id, 10**value)
 
def set_sbml_for_attempt(
    sbml: s2s.SBMLDoc,
    tissue: str,
    kinetic_constants: dict[ParameterId, float] | None,
    output_constants: dict[ParameterId, float] | None,
    concentrations: dict[SpeciesId, float],
):
    # genera casualmente le concentrazioni e riposiziona le concentrazioni
    if kinetic_constants is None and output_constants is None:
        print("[FATAL ERROR] kinetic constants AND constants are None")
        exit(1)
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
    tissue_name: str,
    p: float,
    S: float,
) -> float:
    """
    Funzione di utilità che prende in input:
    - ng_params: dizionario dei parametri che Nevergrad deve ottimizzare (str -> float)
    - hidden_params: dizionario dei parametri che Nevergrad non vede (fissi o di contesto)
    """
    global attempts, best_results
    start_time = time.time()
    set_sbml_for_attempt(sbml, tissue_name, ng_params.value.get("kinetic_constants"), ng_params.value.get("output_constants"), concentrations)
    attempts += 1
    print(f"[INFO] attempt: {attempts}")
    # print(f"[INFO] params:\n {ng_params}")
    # Le chiavi di ng_params sono sempre 'f', 'k_1', 'k_2'
    result = 0.0
    integration_errors = 0

    (steady_pen, bio_pen) = sim.simulate(sbml, concentrations, tissue_name)
    if math.isnan(steady_pen) or math.isnan(bio_pen):
        print("[FATAL ERORR] utility function returned NaN")
        exit(1)
    if math.isinf(steady_pen):
        if steady_pen < 0:
            print("[FATAL ERROR] Negative infinity returned by simulate function. An unexpected error occurred.")
            exit(1)
        else:
            end_time = time.time()
            print(f"[INFO] Simulation {attempts} ended WITH errors in {end_time - start_time:.2f} seconds")
            return math.inf

    result += bio_pen*(1 - p)*S + steady_pen*p*S
    end_time = time.time()
    print(f"[INFO] Simulation {attempts} ended WITHOUT errors in {end_time - start_time:.2f} seconds")
    utility = result
    print(f"[INFO] utility: {utility}")
    if utility < best_results:
        best_results = utility
    return utility


TISSUE="breast_cancer_cell" # TODO: param

def search_for_kinetic_constants(sbml: s2s.SBMLDoc,file_path: str,  tissue: str, concentrations: dict[SpeciesId, float]) -> dict[ParameterId, float]:
    print("[INFO] Opt kinetic constants") 
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
        output_param_dict[oc] = ng.p.Scalar(lower=-20.0, upper=0.0)
    parametrization = ng.p.Dict(**{ "kinetic_constants": ng.p.Dict(**kinetic_param_dict), "output_constants": ng.p.Dict(**output_param_dict) })
    optimizer = ng.optimizers.NGOpt(parametrization=parametrization, budget=3_000)

    def ng_objective(ng_params):
        # hidden parameters sbml and concentrations
         return utility_function(ng_params, sbml, concentrations, tissue, 0.1, 10**4)
    print("[INFO] Start Opt")
    start_opt = time.time()
    for _ in range(optimizer.budget):
        x = optimizer.ask()
        loss = ng_objective(x)
        optimizer.tell(x, loss)
    end_opt = time.time()
    recommendation = optimizer.provide_recommendation()
    print("[INFO] Best parameters found:\n", recommendation.value)
    print("[INFO] Best utility found:", best_results)
    print(f"[INFO] time:{end_opt - start_opt:.2f}")
    with open("parameters_kinetic_constants.json", "w") as f:
        json.dump(recommendation.value, f)
    result: dict[str, dict[ParameterId, float]] = recommendation.value
    set_sbml_for_attempt(sbml, TISSUE, result["kinetic_constants"], result["output_constants"], concentrations)
    sbml.save_converted_file(file_path.replace(".", "-kinetic-constants."))
    sim.simulate(sbml ,concentrations, tissue,  plot=True, output_file_name="kinetic")
    return result

def main():
    global best_results, attempts
    file_path = parse_args()

    sbml: s2s.SBMLDoc = s2s.SBMLDoc(file_path)
    # ref: https://bionumbers.hms.harvard.edu/bionumber.aspx?id=115154&ver=1&trm=cell+size+breast+cancer+cell+human+&org=
    volume_cell_breast_cancer_cell = 1.76 * 10**12  # nanometers
    set_compartement_size(sbml, volume_cell_breast_cancer_cell)
    proteins: dict[SpeciesId,UniprodId] = sbml.get_proteins_data()
    
    (proteomics, _) = get_proteomics(proteins) # APIs

    concentrations: dict[SpeciesId, float] = convert_ibaq_to_concentrations(sbml, proteomics, TISSUE)
    print(concentrations)
    # replica il modello per il tessuto del cancro al seno
    print("[INFO] init model")
    sbml = init_model(sbml, file_path, TISSUE, concentrations)
    # le costanti cinetiche rappresentato un esponente x. Per avere il valore calcolcare 10^x
    result = search_for_kinetic_constants(sbml, file_path, TISSUE, concentrations)
    
    print("[INFO] kinetic costants: ")
    print(result["kinetic_constants"])
    print("[INFO] output constants")
    print(result["output_constants"])