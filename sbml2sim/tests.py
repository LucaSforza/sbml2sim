import json
import os
import random
import time
from typing import Any, Iterable
import sys
import math

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
    for (species, value) in concentrations.items():
        id = tissue+"_"+species
        print(f"for species {id} the mol/L is {value}")
        new_sbml.set_initial_concentration(id, value)
    # TODO: add start concentration for input casual
    new_sbml.add_time_to_model()
    new_sbml.add_avg_calculations_for_all_species()
    new_sbml.save_converted_file(path.replace(".","-real-tissues-modified."))
    new_sbml.simulate(output_file="real-tissue-"+RESULTS, duration=DURATION)
    plot_results("real-tissue-"+RESULTS, "real-tissues.png")
    print("[INFO] clone model per tissue simulation completed")

def test_all(sbml: s2s.SBMLDoc, tissue: str, path: str, concentrations: dict[str, dict[str, float]]):
    print(f"[INFO] cloning SBML for tissue: {tissue}")
    sbml.save_converted_file(path.replace(".","-modified."))
    test_clone_model_per_tissue(sbml, tissue, path, concentrations)
    # sbml.add_time_to_model()
    # sbml.add_avg_calculations_for_all_species()
    # print("[INFO] simulate SBML, every species as a random start concentration")
    # test_random_start_concentration(sbml)
    # print("[INFO] simulate SBML, every protein and compound as a random start concentration")
    # test_random_protein_compound_start_concentration(sbml)
    # print("[INFO] end tests")



# @returns map species, concentration mol/L
def convert_ibaq_to_concentrations(sbml: s2s.SBMLDoc, proteomics: dict[str,tuple[str, list[proteomic]]] ,tissue: str) -> dict[str, float]:
    # reference: https://book.bionumbers.org/how-many-proteins-are-in-a-cell/
    proteins_in_a_cell = 2.0*(10**(-10))
    result = {}
    
    # calcola l'intensità totale del tessuto 
    total_intensity = 0.0
    for species_id, (_, tissue_list) in proteomics.items():
        for prot in tissue_list:
            if ptc.get_tissue_name(prot) == tissue:
                total_intensity += ptc.get_intensity(prot)
    print(f"Total intensity for tissue '{tissue}': {total_intensity}")
    
    for species_id, (_, tissue_list) in proteomics.items():
        compartment_id = sbml.get_compartement(species_id)
        volume_liters = sbml.get_volume_compartement(compartment_id)
        if volume_liters < 1e-15:
            print(f"[WARNING] Volume sospetto per compartimento {compartment_id}: {volume_liters} L")
        atomic_weight = ptc.get_mol_weight(tissue_list[0])
        tissue_conc = None
        for prot in tissue_list:
            if ptc.get_tissue_name(prot) == tissue:
                intensity = ptc.get_intensity(prot)
                # calcola la percentuale di presenza nel tessuto
                f = intensity/total_intensity
                m = f*proteins_in_a_cell # calcola la mole della specie
                n = m/atomic_weight # dividilo per peso atomico,cosi ad avere la mole delle singole proteine
                tissue_conc = n/volume_liters # calcola la mole per litro
                result[species_id] = tissue_conc
                break
    return result

# TODO: questa cosa è stupida, semplicemente inserisci le unità di misura nel file SBML
def nanometers_to_liters(x: float) -> float:
    return x*(10**(-24))

def volume(r: float) -> float:
    return (4.0/3.0)*math.pi*(r**3)

def set_compartement_size(sbml: s2s.SBMLDoc):
    diameter_plasma_membrane = 10.0
    # ref: https://bionumbers.hms.harvard.edu/bionumber.aspx?id=115154&ver=1&trm=cell+size+breast+cancer+cell+human+&org=
    diameter_cell = 1.76 * 10**12  # nanometers
    volume_cell = 1.76 * 10**12  # nanometers
    volume_plasma_membrane = volume(diameter_plasma_membrane)
    # nucleo occupa 20% del volume interno
    volume_nucleoplasm = .2*(volume_cell - volume_plasma_membrane)
    volume_cytosol = (volume_cell - volume_plasma_membrane) - volume_nucleoplasm
    
    for i in range(sbml.get_num_compartements()):
        name: str = sbml.get_name_compartement(i)
        match name:
            case "plasma membrane":
                sbml.set_volume_compartement(i, nanometers_to_liters(volume_plasma_membrane))
                pass
            case "cytosol":
                sbml.set_volume_compartement(i, nanometers_to_liters(volume_cytosol))
                pass
            case "nucleoplasm":
                sbml.set_volume_compartement(i, nanometers_to_liters(volume_nucleoplasm))
            case "extracellular region":
                sbml.set_volume_compartement(i, 7.0* 10**12)
            case _:
                print(f"[FATAL ERROR] compartement {name} doen't exists")
                exit(1)

def get_proteomics(proteins: dict[str,str]) -> tuple[Any, Any]:
    proteomics: dict[str,tuple[str, list[proteomic]]] = dict()
    all_tissue_names = set()
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
    return (proteomics, all_tissue_names)

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
    set_compartement_size(sbml)
    
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