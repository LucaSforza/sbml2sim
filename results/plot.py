import sys
import pandas as pd
import matplotlib.pyplot as plt

def plot_concentration(csv: pd.DataFrame, tissue_concentration: dict[str, float], name: str, stripped_name: str):
    print(f"[INFO] start plotting {name}")
    plt.plot(csv["time"], csv[name], label=name, linewidth=2.5)
    plt.axhline(y=tissue_concentration[stripped_name], color='r', linestyle='--', label=f"{stripped_name} (real value)")
    plt.title(f"{stripped_name} concentration over time")
    plt.xlabel("Time")
    plt.ylabel("mol/L")
    plt.yscale('log')
    plt.ylim(1e-5, 1)
    plt.legend([f"{stripped_name} mol/L"])
    plt.legend()
    plt.savefig(stripped_name+".png")
    plt.clf()
    
def plot_utility(log_file_path: str):
    attempt = 0
    with open(log_file_path) as f:
        for line in f.readlines():
            if line.startswith("[INFO] utility:"):
                attempt += 1
                utility = float(line.split(":")[-1])
                plt.scatter(attempt, utility, color='blue', s=5)
        plt.yscale("log")
        plt.xlabel("attempt")
        plt.ylabel("utility")
        plt.savefig("log.png")
        plt.clf()

def main():
    # TODO: riplotta kinetic.png
    # TODO: prendi tutte le informazioni dai log
    tissue = "breast_cancer_cell"
    concentration = {
        'species_6813727': 0.0003225823634385977, 
        'species_391267': 0.00028469824008714093,
        'species_379540': 0.00010945937352407899,
        'species_201861': 0.0005121705189427184,
        'species_5216230': 0.0007163125346026875,
        'species_379539': 0.00025669755152323703,
        'species_379537': 0.00037514093686477336,
        'species_379538': 0.001604830448462479,
        'species_379546': 0.0006385251942726382
    }
    
    tissue_concentration: dict[str, float] = dict()
    
    for (key, value) in concentration.items():
        tissue_concentration[tissue+"_"+key] = value
    print(tissue_concentration)
    
    csv = pd.read_csv("kinetic.csv")
    for col in csv.columns:
        stripped_col = col.strip("[]")
        if stripped_col in tissue_concentration.keys():
            plot_concentration(csv, tissue_concentration, col, stripped_col)

    print(f"[INFO] start plotting all concentration")
    for col in csv.columns:
        if col.startswith("avg_"):
            plt.plot(csv["time"], csv[col])
            plt.title(f"Concentration over time")
            plt.xlabel("Time")
            plt.ylabel("mol/L")
            # plt.yscale('log')
    plt.savefig("kinetic2.png")
    plt.clf()       
    
# TODO: ristampa kinetic.png con meno cose nella leggenda e in base logaritmica
if "__main__" == __name__:
    main()
    # plot_utility("log")