import nevergrad as ng
import pandas as pd
import matplotlib.pyplot as plt

import sbmlconverter as sc

import sys
import getopt

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
        # Use gnu_getopt to allow options after positional arguments
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
        "png_outpath": "data.png",
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
                options["png_outpath"] = a
            case "-h":
                usage()
                sys.exit(0)
            case _:
                usage()
                sys.exit(1)
    return options
    

if __name__ == "__main__":
    options = parse_args()

    sbml = sc.SBMLDoc(options["sbml_path"], options["convinience"])
    if options["random"]:
        sbml.random_start_concentration()
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