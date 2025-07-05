import json
import os
import random
import time
from typing import Any, Iterable
import sys
import math

import pandas as pd
import matplotlib.pyplot as plt

import s2s
import simulate as sim
import proteomic as ptc
from proteomic import proteomic
from bioutils import *

DURATION = 10.0
RESULTS = "results.csv"

def plot_results(csv_path: str, output_file_path: str):
    df = pd.read_csv(csv_path)
    for col in df.columns:
        if col != "time" and col != "get_time":
            plt.plot(df["time"], df[col], label=col)
    plt.legend()
    plt.savefig(output_file_path)
    plt.clf()
    plt.close()

def test_random_start_concentration(sbml: s2s.SBMLDoc):
    sbml.random_start_concentration()
    sbml.simulate(output_file="rand-start-"+RESULTS, duration=DURATION)
    plot_results("rand-start-"+RESULTS, "random_start_concentration.png")
    print("[INFO] random start concentration simulation completed")

def test_random_protein_compound_start_concentration(sbml: s2s.SBMLDoc):
    sbml.small_compound_start_random_concentration()
    sbml.random_protein_concentrations()
    sbml.simulate(output_file="rand-prot-com-"+RESULTS, duration=DURATION)
    plot_results("rand-prot-com-"+RESULTS, "random_protein_compound_start_concentration.png")
    print("[INFO] random protein and compound simulation completed")

def test_clone_model_per_tissue(sbml: s2s.SBMLDoc, tissue: str, path: str, concentrations: dict[str, float]):
    
    new_sbml = sbml.replicate_model_per_tissue([tissue])
    new_sbml.random_start_concentration()
    assign_concentrations_to_small_compound(new_sbml)
    for (species, value) in concentrations.items():
        id = tissue+"_"+species
        print(f"for species {id} the mol/L is {value}")
        new_sbml.set_initial_concentration(id, value)
    # TODO: add start concentration for input casual
    new_sbml.add_time_to_model()
    new_sbml.add_avg_calculations_for_all_species()
    new_sbml.save_converted_file(path.replace(".","-real-tissues-modified."))
    print(new_sbml.get_kinetic_constants())
    print(new_sbml.get_output_constants())
    sim.simulate(new_sbml, concentrations, tissue, plot=True)
    print("[INFO] clone model per tissue simulation completed")

def test_all(sbml: s2s.SBMLDoc, tissue: str, path: str, concentrations: dict[str, dict[str, float]]):
    print(f"[INFO] cloning SBML for tissue: {tissue}")
    sbml.save_converted_file(path.replace(".","-modified."))
    test_clone_model_per_tissue(sbml, tissue, path, concentrations)


def parse_args() -> tuple[str,str]:
    import sys
    if len(sys.argv) != 3:
        print(sys.argv[0], "<reactome SBML> <tissue>")
    sbml_path = sys.argv[1]
    tissue = sys.argv[2]
    return (sbml_path, tissue)

def main():
    random_seed = int(time.time() * 1000)
    random.seed(random_seed)
    s2s.set_seed(random_seed)
    (sbml_path, tissue) = parse_args()
    sbml = s2s.SBMLDoc(sbml_path)
    # ref: https://bionumbers.hms.harvard.edu/bionumber.aspx?id=115154&ver=1&trm=cell+size+breast+cancer+cell+human+&org=
    volume_cell_breast_cancer_cell = 1.76 * 10**12  # nanometers
    set_compartement_size(sbml, volume_cell_breast_cancer_cell)
    proteins: dict[str,str] = sbml.get_proteins_data()

    (proteomics, all_tissue_names) = get_proteomics(proteins)
    if len(proteomics) == 0:
        print("[FATAL ERROR] proteomics are void")
        exit(1)
    tissue_names = choose_tissue_for_replication(all_tissue_names, proteomics)
    if not tissue in tissue_names:
        print("[FATAL ERROR] breast cancer cell not avaible")
    sbml.add_kinetic_laws_if_not_exists()
    sbml.random_kinetic_costant_value()
    sbml.set_zero_output_costant()
    
    concentrations = convert_ibaq_to_concentrations(sbml, proteomics, tissue)
    test_all(sbml, tissue, sbml_path, concentrations)

if __name__ == "__main__":
    main() 