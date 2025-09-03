import json
import sys
import math
import secrets
import time
from typing import Any

import nevergrad as ng

from optimizer import *

from bioutils import convert_ibaq_to_concentrations, get_proteomics, assign_concentrations, set_compartement_size, SpeciesId, ParameterId, UniprodId
import s2s as s2s
from proteomic import proteomic
import proteomic as ptc
import simulate as sim
import random as _py_random
import numpy as _np
import os as _os

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
    sbml: s2s.Simulator,
    parameters: dict[ParameterId, float] | None
):
    if parameters is None: return
    for (param_id, value) in parameters.items():
        sbml.set_parameter(param_id, 10**value)
 
def set_sbml_for_attempt(
    sbml: s2s.Simulator,
    tissue: str,
    kinetic_constants: dict[ParameterId, float] | None,
    output_constants: dict[ParameterId, float] | None,
    input_constants: dict[ParameterId, float] | None,
    concentrations: dict[SpeciesId, float],
):
    # genera casualmente le concentrazioni e riposiziona le concentrazioni
    if kinetic_constants is None and output_constants is None and input_constants is None:
        print("[FATAL ERROR] kinetic constants AND constants are None")
        exit(1)
    # sbml.random_start_concentrations()
    # TODO: reset initial value for avg non-constant parameters
    assign_parameters(sbml, kinetic_constants)
    assign_parameters(sbml, output_constants)
    for (species_id, value) in input_constants.items():
        sbml.set_initial_concentration(species_id, 10**value)

# Global variable to count attempts
attempts = 0
best_results = math.inf

def utility_function(
    ng_params: dict[str, dict[ParameterId, float]],
    simulator: s2s.ParallelSimulator,
    model: s2s.Simulator,
    error_handler: s2s.ErrorHandler,
    concentrations: dict[SpeciesId, float],
    tissue_name: str,
) -> float:
    """
    Funzione di utilità che prende in input:
    - ng_params: dizionario dei parametri che Nevergrad deve ottimizzare (str -> float)
    - hidden_params: dizionario dei parametri che Nevergrad non vede (fissi o di contesto)
    """
    global attempts, best_results
    start_time = time.time()
    set_sbml_for_attempt(model, tissue_name, ng_params.value.get("kinetic_constants"), ng_params.value.get("output_constants"), ng_params.value.get("input_constants"), concentrations)
    attempts += 1
    print(f"[INFO] attempt: {attempts}")

    (loss, is_error) = simulator.simulate(error_handler)[0]

    if math.isnan(loss):
        print("[FATAL ERORR] utility function returned NaN")
        exit(1)
    stringa = "WITHOUT"
    if is_error: 
        stringa = "WITH"
        loss = math.inf
    end_time = time.time()
    if loss < best_results:
        best_results = loss
    print(f"[INFO] Simulation {attempts} ended {stringa} errors in {end_time - start_time:.2f} seconds")
    print(f"[INFO] utility: {loss}")
    print(f"[INFO] best utility: {best_results}")
    print(f"[INFO] params: {ng_params.value}")
    return loss


TISSUE="breast_cancer_cell" # TODO: param

def search_for_kinetic_constants(budget: str, optimizer: str,sbml: s2s.SBMLDoc, model: s2s.Simulator, simulator: s2s.ParallelSimulator, error_handler: s2s.ErrorHandler) -> dict[ParameterId, float]:
    print("[INFO] Opt kinetic constants") 
    kinetic_constants: list[ParameterId] = sbml.get_kinetic_constants()
    # output_constants: list[ParameterId] = sbml.get_output_constants()
    # input_species: list[ParameterId] = [species for species in sbml.get_input_species() if species not in known_input]
    
    # Build parametrization with separate dictionaries
    kinetic_param_dict = {}
    for kc in kinetic_constants:
        # Parameter: kinetic constant
        if kc.startswith("k_input"):
            kinetic_param_dict[kc] = ng.p.Scalar(lower=-40.0, upper=-0.0)
        else:
            kinetic_param_dict[kc] = ng.p.Scalar(lower=-20.0, upper=20.0)
        
    opt = NevergradOpt(model, simulator, error_handler)
    
    
    parametrization = ng.p.Dict(**{ "kinetic_constants": ng.p.Dict(**kinetic_param_dict)})
    
    args = {"params": parametrization}
    
    if optimizer:
        args["optimizer"] = optimizer
    
    if budget:
        args["budget"] = int(budget)
    
    
    start = time.time()
    result, best = opt.minimize(args)
    print(f"[INFO] ended optimitation in {time.time() - start}")
    print(f"[INFO] best value: {best}")
    print("[INFO] result (JSON):")
    print(json.dumps(result, indent=2, default=lambda o: str(o)))
    
    return result
    """optimizer = ng.optimizers.NGOpt(parametrization=parametrization, budget=3_000)

    def ng_objective(ng_params):
        # hidden parameters sbml and concentrations
         return utility_function(ng_params, simulator,model, error_handler, concentrations, tissue)
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
    return result"""

def main(args: Any):
    global best_results, attempts
    # Imposta un seed casuale per questa esecuzione e rendilo riproducibile dove possibile
    random_seed = secrets.randbits(64)
    print(f"[INFO] Random seed: {random_seed}")
    _py_random.seed(random_seed)
    try:
        _np.random.seed(int(random_seed % (2**32)))
    except Exception:
        pass
    _os.environ["PYTHONHASHSEED"] = str(random_seed)
    try:
        s2s.set_seed(random_seed)
    except Exception:
        pass
    sbml: s2s.SBMLDoc = s2s.SBMLDoc(args.input_file)
    model = s2s.rr_simualtor(sbml)
    simualtor = s2s.ParallelSimulator(1)
    simualtor.add_worker(model) # TODO: cambia nome funzione
    # TODO: error handler
    error_handler = s2s.error_sum_create()
    s2s.error_sum_add_handler(error_handler, s2s.stability_error_create(), 10)
    s2s.error_sum_add_handler(error_handler, s2s.transitorial_error_create(), 1)
    result = search_for_kinetic_constants(args.budget, args.optimizer, sbml, model, simualtor, error_handler)
    print("[INFO] kinetic costants: ")
    print(result["kinetic_constants"])
    
    with open("ord_params.json", "w") as f:
        json.dump(result["kinetic_constants"],f)