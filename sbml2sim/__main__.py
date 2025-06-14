# from opt import main

# main()

from proteomic import *

proteomics = get_tissue("P00533")
# print("Proteomics for protein: P00533")
# for proteomic in proteomics:
#     print(f"    tissue id: {get_tissue_id(proteomic)}")
#     print(f"    tissue name: {get_tissue_name(proteomic)}")
#     print(f"    intensity: {get_intensity(proteomic)}")
#     print("-----------------------------------------------")

for tissue in get_all_tissue_names(proteomics):
    print(tissue)