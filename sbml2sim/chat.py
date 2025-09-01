import roadrunner as te
import json
import matplotlib.pyplot as plt

# === 1. Caricamento del modello ===
sbml_file = "converted_R-HSA-391251.sbml"
rr = te.loadSBMLModel(sbml_file)

# === 2. Simulazione ===
result = rr.simulate(0, 100, 1000)  # da 0 a 100 con 1000 punti

# === 3. Salvataggio traiettorie in grafico ===
plt.figure(figsize=(12, 6))
for i, species in enumerate(rr.model.getFloatingSpeciesIds()):
    plt.plot(result[:, 0], result[:, i+1], label=species)

plt.xlabel("Time")
plt.ylabel("Concentration")
plt.legend(loc="best", fontsize=6)
plt.tight_layout()
plt.savefig("trajectories.png", dpi=300)
plt.close()

# === 4. Riassunto delle concentrazioni iniziali e finali ===
species_ids = rr.model.getFloatingSpeciesIds()
initial_conc = {s: rr.model.getFloatingSpeciesInitialConcentration(i) 
                for i, s in enumerate(species_ids)}
final_conc = {s: result[-1, i+1] for i, s in enumerate(species_ids)}

summary = {
    "initial_concentrations": initial_conc,
    "final_concentrations": final_conc
}

with open("simulation_summary.json", "w") as f:
    json.dump(summary, f, indent=4)

print("✅ Simulazione completata.")
print("• Grafico salvato in trajectories.png")
print("• Riassunto salvato in simulation_summary.json")