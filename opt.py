import nevergrad as ng
import pandas as pd
import matplotlib.pyplot as plt

import sbmlconverter as sc

sbml = sc.SBMLDoc("sbmls/R-HSA-1643713.sbml")
sbml.random_start_concentration()
sbml.save_converted_file("R-HSA-1643713-modified.sbml")
sbml.simulate(duration=1000)

df = pd.read_csv("simulation_results.csv")
for col in df.columns:
    if col != "time" and col.startswith("avg_"):
        plt.plot(df["time"], df[col], label=col)
plt.legend()
plt.savefig("data.png")