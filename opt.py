import nevergrad as ng
import pandas as pd
import matplotlib.pyplot as plt

import sbmlconverter as sc
import uniprod

import sys
import getopt
from typing import Any

PROGRAM_NAME = sys.argv[0]

def usage():
    print(f"{PROGRAM_NAME} <SBML file path> [OPTIONS]")
    print("    -d <horizon>                 Horizon for the simulation, default: 100")
    print("    -o <output path>             Output path for the new SBML file, default: a.out")
    print("    -s <simulation results path> Indicate the path for saving the simulation results, default: simulation_results.csv")
    print("    -r                           Random start concentrations for all the species")
    print("    -c                           All convinience law")
    print("    -h                           Print those informations")
    print("    -g                           Dump protein species and their corresponding genes")
    print("    -f <output path>             Path to save the simulations result plot")
    
def parse_args():
    try:
        opts, args = getopt.gnu_getopt(
            sys.argv[1:], 
            "d:o:s:rchgf:", 
            []
        )
    except getopt.GetoptError as err:
        print(err)
        usage()
        sys.exit(2)

    # Defaults
    options = {
        "horizon": 100,
        "output": "a.out",
        "simresults": "simulation_results.csv",
        "random": False,
        "convinience": False,
        "dumps": False,
        "png_outputpath": "data.png",
    }

    # First positional argument: SBML file path
    if not args:
        usage()
        sys.exit(2)
    options["sbml_path"] = args[0]
    
    for o, a in opts:
        match o:
            case "-d":
                options["horizon"] = int(a)
            case "-o":
                options["output"] = a
            case "-s":
                options["simresults"] = a
            case "-r":
                options["random"] = True
            case "-c":
                options["convinience"] = True
            case "-g":
                options["dumps"] = True
            case "-f":
                options["png_outputpath"] = a
            case "-h":
                usage()
                sys.exit(0)
            case _:
                usage()
                sys.exit(1)
    return options

def get_all_ids(genes_data: dict[str, list[str]]) -> set[str]:
    result = set()
    for _,ids in genes_data.items():
        result.update(ids)
    return result

def map_ids_to_genes(genes_data: dict[str, list[str]]) -> Any | None:
    genes_id = get_all_ids(genes_data)
    
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
    
    
def main():
    options = parse_args()

    sbml = sc.SBMLDoc(options["sbml_path"], options["convinience"])
    if options["random"]:
        # sbml.random_start_concentration()
        sbml.random_protein_concentrations()
    if options["dumps"]:
        sbml.dump_genes_data()
    sbml.save_converted_file(options["output"])
    sbml.simulate(output_file=options["simresults"], duration=options["horizon"])

    df = pd.read_csv(options["simresults"])
    for col in df.columns:
        if col != "time" and col.startswith("avg_"):
            plt.plot(df["time"], df[col], label=col)
    plt.legend()
    plt.savefig(options["png_outputpath"])
    
    results = map_ids_to_genes(sbml.get_genes_data())
    if results is not None:
        for result in results:
            _from = result["from"]
            _to = result["to"]
            genes = _to["genes"]
            if len(genes) > 1:
                print("[WARNING] a id is assosieted to more than one gene")
            for gene in genes:
                name = gene["geneName"]["value"]
                print(f"{_from} -> {name}")
    

if __name__ == "__main__":
    main()