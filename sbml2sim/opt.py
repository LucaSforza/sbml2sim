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
    
def main():
    file_path = parse_args()
    sim_output = "result.csv"

    sbml: s2s.SBMLDoc = s2s.SBMLDoc(file_path)
    proteins: dict[str,str] = sbml.get_proteins_data()

    tissue_names = set()
    proteomics: dict[str,tuple[str, list[proteomic]]] = dict()
    for species, protein in proteins.items():
        tissues = ptc.get_tissue(protein)
        proteomics[species] = (protein,tissues)
        tissue_names.update(ptc.get_all_tissue_names(tissues))
    
    save_proteomics(proteomics)
    new_sbml = sbml.replicate_model_per_tissue(tissue_names)
    new_sbml.save_converted_file(file_path.replace(".","-modified."))
    
            
    
    # sbml.random_protein_concentrations()
    # sbml.save_converted_file("a.sbml")
    # sbml.simulate(output_file=sim_output)

    # df = pd.read_csv(sim_output)
    # for col in df.columns:
    #     if col != "time" and col.startswith("avg_"):
    #         plt.plot(df["time"], df[col], label=col)
    # plt.legend()
    # plt.savefig(options["png_outputpath"])
    
    