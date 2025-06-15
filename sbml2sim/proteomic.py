from typing import Any, Iterable
import requests
from requests.adapters import HTTPAdapter, Retry

proteomic = dict[str, Any]

POLLING_INTERVAL = 3
API_URL = "https://www.proteomicsdb.org/proteomicsdb/logic/api"

retries = Retry(total=5, backoff_factor=0.25, status_forcelist=[500, 502, 503, 504])
session = requests.Session()
session.mount("https://", HTTPAdapter(max_retries=retries))

def get_proteomic_string_request(uniprod_id: str) -> str:
    return f"{API_URL}/proteinexpression.xsodata/InputParams(PROTEINFILTER='{uniprod_id}',MS_LEVEL=1,TISSUE_ID_SELECTION='',TISSUE_CATEGORY_SELECTION='tissue;fluid',SCOPE_SELECTION=1,GROUP_BY_TISSUE=1,CALCULATION_METHOD=0,EXP_ID=-1)/Results?$select=TISSUE_ID,TISSUE_NAME,UNNORMALIZED_INTENSITY&$format=json"

def get_intensity(proteomic: proteomic) -> float:
    return float(proteomic['UNNORMALIZED_INTENSITY'])

def get_tissue_id(proteomic: proteomic) -> str:
    return proteomic['TISSUE_ID']

def get_tissue_name(proteomic: proteomic) -> str:
    name: str = proteomic['TISSUE_NAME']
    return name.replace(' ','_')

def get_all_tissue_names(proteomics: list[proteomic]) -> set[str]:
    result = set()
    for proteomic in proteomics:
        result.add(get_tissue_name(proteomic))
    return result

def get_all_tissue_id(proteomics: list[proteomic]) -> set[str]:
    result = set()
    for proteomic in proteomics:
        result.add(get_tissue_id(proteomic))
    return result

def get_tissue_names_from_bto(proteomics: list[proteomic], btos: set[str]) -> set[str]:
    result = set()
    for proteomic in proteomics:
        if get_tissue_id(proteomic) in btos:
            result.add(get_tissue_name(proteomic))
    return result

def print_proteomic(proteomic: proteomic, prefix: str = "") -> None:
    print(prefix, "Tissue name: ", get_tissue_name(proteomic))
    print(prefix, "Tissue id:   ", get_tissue_id(proteomic))
    print(prefix, "Value:       ", get_intensity(proteomic))

# @returns list of proteomics
def get_tissue(uniprod_id: str) -> list[proteomic]:
    print(f"[INFO] requesting proteomics for protein: {uniprod_id}")
    response = requests.get(get_proteomic_string_request(uniprod_id))
    if response.status_code != 200:
        print(f"[ERROR] status code of the request: {response.status_code}")
        exit(1)
    print("[INFO] completed the request")
    return response.json()['d']['results']