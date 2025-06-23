# import nevergrad as ng
# import pandas as pd
# import matplotlib.pyplot as plt

import sbml2sim as s2s
import uniprod
import proteomic as ptc
from proteomic import proteomic

import sys
from typing import Any
import json
from collections import Counter
import os
import random
import time

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

def get_all_ids(protein_data: dict[str, str]) -> set[str]:
    result = set()
    for _,ids in protein_data.items():
        result.add(ids)
    return result

def map_proteins_to_genes(protein_data: dict[str, str]) -> Any | None:
    genes_id = get_all_ids(protein_data)
    
    job_id = uniprod.submit_id_mapping("UniProtKB_AC-ID", "UniProtKB", genes_id)
    
    if uniprod.check_id_mapping_results_ready(job_id):
        link = uniprod.get_id_mapping_results_link(job_id)
        results = uniprod.get_id_mapping_results_search(link)
        if results.get('failedIds') is not None:
            print(f"[WARNING] failed ids: {results['failedIds']}")
        return results['results']
    else:
        print("[FATAL ERROR] failed request")
        return None
    
def get_map_protein_gene(protein_data: dict[str, str]) -> dict[str, str]:
    results = map_proteins_to_genes(protein_data)
    final_result = dict()
    if results is not None:
        for result in results:
            _from = result["from"]
            _to = result["to"]
            genes = _to["genes"]
            if len(genes) > 1:
                print("[WARNING] an ID is associated with more than one gene")
            for gene in genes:
                name = gene["geneName"]["value"]
                final_result[_from] = name
    return final_result

def save_proteomics(
    proteomics: dict[str, tuple[str, list[proteomic]]],
    file_name: str = "proteomics.json"
) -> None:
    with open(file_name, "w") as f:
        json.dump(proteomics, f)
        
def choose_tissue_for_replication(all_tissue_names: set[str], proteomics: dict[str,tuple[str, list[proteomic]]]) -> set[str]:
    tissue_names = set()
    
    for name in all_tissue_names:
        invariant = True
        for (_,(_,ps)) in proteomics.items():
            found = False
            for p in ps:
                if ptc.get_tissue_name(p) == name:
                    found = True
                    break
            if not found:
                invariant = False
                break
        if invariant:
            tissue_names.add(name)
    
    if not tissue_names:
        # Trova il/i tessuto/i più frequente/i tra tutte le proteine

        all_tissues = []
        for _, (_, ps) in proteomics.items():
            all_tissues.extend([ptc.get_tissue_name(p) for p in ps])

        if not all_tissues:
            print("[ERROR] Nessun tessuto trovato tra le proteine.")
            exit(1)

        counter = Counter(all_tissues)
        max_count = max(counter.values())
        most_common_tissues = {tissue for tissue, count in counter.items() if count == max_count}
        tissue_names = most_common_tissues
    return tissue_names
def main():
    file_path = parse_args()
    sim_output = "result.csv"

    sbml: s2s.SBMLDoc = s2s.SBMLDoc(file_path)
    
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
    
    tissue_names = choose_tissue_for_replication(all_tissue_names, proteomics)
    print(tissue_names)
    save_proteomics(proteomics)
    random_seed = int(time.time() * 1000)
    random.seed(random_seed)
    s2s.set_seed(random_seed)
    sbml.add_kinetic_laws_if_not_exists()
    for i in range(sbml.get_num_compartements()):
        name: str = sbml.get_name_compartement(i)
        match name:
            case "plasma membrane":
                # Calcola il volume di una sfera di raggio 10 nanometri (nm)
                # raggio_nm = 10
                # raggio_m = raggio_nm * 1e-9  # converte in metri
                # volume_m3 = (4/3) * 3.141592653589793 * (raggio_m ** 3)
                # sbml.set_volume_compartement(i,volume_m3)
                sbml.set_volume_compartement(i, 12.6)
                pass
            case "extracellular region":
                # volume_entire_cell = 11 * 1e24  # nanometri cubi
                # volume_extracellular = volume_entire_cell * 0.3
                # sbml.set_volume_compartement(i, volume_extracellular)
                sbml.set_volume_compartement(i, 2720)
            case "cytosol":
                volume_entire_cell = 2250
                volume_cytosol = volume_entire_cell * 0.7
                sbml.set_volume_compartement(i, volume_cytosol)
                sbml.set_volume_compartement(i, 1000)
                pass
            case _:
                print(f"[FATAL ERROR] compartement {name} doen't exists")
                exit(1)
    sbml.save_converted_file(file_path.replace(".","-modified."))
    new_sbml = sbml.replicate_model_per_tissue(tissue_names)
    new_sbml.random_kinetic_costant_value()
    new_sbml.small_compound_start_random_concentration()
    new_sbml.random_protein_concentrations()
    new_sbml.add_time_to_model()
    new_sbml.add_avg_calculation_for_all_proteins()
    new_sbml.simulate(sim_output, duration=1000.0)
    new_sbml.save_converted_file(file_path.replace(".","-tissues-modified."))
    
    