from collections import Counter
import json
import math
import os
from typing import Any

from proteomic import proteomic
import s2s
import proteomic as ptc
import uniprod

AVOGRADO = 6.022 * 10**23

def assign_concentrations_to_small_compound(sbml: s2s.SBMLDoc):
    # compounds è una mappa species_id all'id CHEBI
    # ATP https://hmdb.ca/metabolites/HMDB0000538
    # ATP concentrazione uguale in tutti i compartimenti
    # 0.00154 mol/L
    
    # H2O https://hmdb.ca/metabolites/HMDB0002111
    # 55 mol/L
    
    # ADP https://hmdb.ca/metabolites/HMDB0001341
    # 2.7*(10**-4)
    
    # PI(4,5)P2 https://pubmed.ncbi.nlm.nih.gov/33441034/ boh forse 0.005 mol/L
    
    compounds: dict[str, str] = sbml.get_compounds_data()
    
    for (species_id, chebi_id) in compounds.items():
        chebi_id = int(chebi_id)
        if chebi_id == 30616:
            # ATP
            sbml.set_initial_concentration(species_id, 0.00154)
        elif chebi_id == 15377:
            # Water
            sbml.set_initial_concentration(species_id, 55)
            pass
    
    pass

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
        id = sbml.get_compartement(species_id)
        volume_liters = sbml.get_volume_compartement(id)
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

def volume(d: float) -> float:
    """
    Calcola il volume di una sfera dato il diametro d.
    d è il diametro.
    """
    r = d / 2.0  # r ora rappresenta il diametro, quindi lo divido per 2 per ottenere il raggio
    return (4.0/3.0)*math.pi*(r**3)

def set_compartement_size(sbml: s2s.SBMLDoc, volume_cell_nanometers: float):
    diameter_plasma_membrane = 10.0
    # ref: https://bionumbers.hms.harvard.edu/bionumber.aspx?id=115154&ver=1&trm=cell+size+breast+cancer+cell+human+&org=
    volume_cell = 1.76 * 10**12  # nanometers
    volume_plasma_membrane = volume(diameter_plasma_membrane)
    # nucleo occupa 20% del volume interno
    volume_nucleoplasm = .2*(volume_cell - volume_plasma_membrane)
    volume_cytosol = (volume_cell - volume_plasma_membrane) - volume_nucleoplasm
    
    for i in range(sbml.get_num_compartements()):
        name: str = sbml.get_name_compartement(i)
        if name == "plasma membrane":
            sbml.set_volume_compartement(i, nanometers_to_liters(volume_plasma_membrane))
        elif name == "cytosol":
            sbml.set_volume_compartement(i, nanometers_to_liters(volume_cytosol))
        elif name == "nucleoplasm":
            sbml.set_volume_compartement(i, nanometers_to_liters(volume_nucleoplasm))
        elif name == "extracellular region":
            sbml.set_volume_compartement(i, nanometers_to_liters(7*10**12))
        else:
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