import matplotlib

import sys
import pandas as pd
import matplotlib.pyplot as plt


PROGRAM_NAME = sys.argv[0]

def usage():
    print(PROGRAM_NAME+" <file path to simulation_results.csv>")

# return path to the SBML file and output file
def parse_args() -> str:
    if len(sys.argv) < 2:
        usage()
        exit(1)
    results_path = sys.argv[1]
    return results_path

def main():
    results_path = parse_args()

    df = pd.read_csv(results_path)
    for col in df.columns:
        if col.startswith("avg_"):
            plt.plot(df["time"], df[col], label=col)
    plt.legend()
    plt.savefig("data.png")
    

if __name__ == "__main__":
    main()