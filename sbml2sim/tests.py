import json
import os
import random
import time
from typing import Iterable
import sys

import pandas as pd
import matplotlib.pyplot as plt

import sbml2sim as s2s
import proteomic as ptc
from proteomic import proteomic
from opt import choose_tissue_for_replication, save_proteomics

DURATION = 200.0
RESULTS = "results.csv"

def plot_results(csv_path: str, output_file_path: str):
    df = pd.read_csv(csv_path)
    for col in df.columns:
        if col.startswith("avg_"):
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

def test_clone_model_per_tissue(sbml: s2s.SBMLDoc, tissues: Iterable[str], path: str):
    new_sbml = sbml.replicate_model_per_tissue(tissues)
    new_sbml.add_time_to_model()
    new_sbml.add_avg_calculation_for_all_proteins()
    new_sbml.random_start_concentration()
    new_sbml.simulate(output_file="tissue-"+RESULTS, duration=DURATION)
    new_sbml.save_converted_file(path.replace(".","-tissues-modified."))
    plot_results("tissue-"+RESULTS, "tissues.png")
    print("[INFO] clone model per tissue simulation completed")

def test_all(sbml: s2s.SBMLDoc, tissues: Iterable[str], path: str):
    print(f"[INFO] cloning SBML for each tissue: {tissues}")
    test_clone_model_per_tissue(sbml, tissues, path)
    sbml.add_time_to_model()
    sbml.add_avg_calculation_for_all_proteins()
    sbml.save_converted_file(path.replace(".","-modified."))
    print("[INFO] simulate SBML, every protein and compound as a random start concentration")
    test_random_protein_compound_start_concentration(sbml)
    print("[INFO] simulate SBML, every species as a random start concentration")
    test_random_start_concentration(sbml)
    print("[INFO] end tests")

def nanometers_to_liters(x: float) -> float:
    return x*(10**(-24))

import math

def volume(r: float) -> float:
    return (4.0/3.0)*math.pi*(r**3)

def set_compartement_size(sbml: s2s.SBMLDoc):
    diameter_plasma_membrane = 10.0
    diameter_cell = 10000.0 #nano meters
    volume_cell = volume(diameter_cell)
    volume_plasma_membrane = volume(diameter_plasma_membrane)
    # nucleo occupa 20% del volume interno
    volume_nucleoplasm = .2*(volume_cell - volume_plasma_membrane)
    volume_cytosol = (volume_cell - volume_plasma_membrane) - volume_nucleoplasm
    
    for i in range(sbml.get_num_compartements()):
        name: str = sbml.get_name_compartement(i)
        match name:
            case "plasma membrane":
                sbml.set_volume_compartement(i, volume_plasma_membrane)
                pass
            case "cytosol":
                sbml.set_volume_compartement(i, volume_cytosol)
                pass
            case "nucleoplasm":
                sbml.set_volume_compartement(i, volume_nucleoplasm)
            case _:
                print(f"[FATAL ERROR] compartement {name} doen't exists")
                exit(1)

def main():
    random_seed = int(time.time() * 1000)
    random.seed(random_seed)
    s2s.set_seed(random_seed)
    
    sbml_path = sys.argv[1]
    sbml = s2s.SBMLDoc(sbml_path)
    set_compartement_size(sbml)
    
    proteins: dict[str,str] = sbml.get_proteins_data()
    all_tissue_names = set()
    proteomics: dict[str,tuple[str, list[proteomic]]] = dict()
    if os.path.exists("proteomics.json"):
        with open("proteomics.json", "r") as f:
            proteomics = json.load(f)
            for species, (protein, tissues) in proteomics.items():
                all_tissue_names.update(ptc.get_all_tissue_names(tissues))
    else:
        for species, protein in proteins.items():
            tissues = ptc.get_tissue(protein)
            proteomics[species] = (protein,tissues)
            all_tissue_names.update(ptc.get_all_tissue_names(tissues))
        save_proteomics(proteomics) 
    tissue_names = choose_tissue_for_replication(all_tissue_names, proteomics)
    sbml.add_kinetic_laws_if_not_exists()
    sbml.random_kinetic_costant_value()
    test_all(sbml, tissue_names, sbml_path)

if __name__ == "__main__":
    main() 