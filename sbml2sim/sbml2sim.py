from ctypes import cdll
from ctypes import c_void_p, c_char_p, c_bool, c_int, c_double
from ctypes import POINTER, c_size_t
from typing import Iterable

lib = cdll.LoadLibrary("build/libsbmlconverter.so")

lib.Proteins_iterator.restype = c_void_p
lib.Proteins_iterator.argtypes = [c_void_p]

lib.Proteins_delete_iterator.restype = None
lib.Proteins_delete_iterator.argtypes = [c_void_p]

lib.Proteins_iterator_next.restype = c_void_p
lib.Proteins_iterator_next.argtypes = [c_void_p]

lib.Pair_delete.restype = None
lib.Pair_delete.argtypes = [c_void_p]

lib.Pair_first_c_str.restype = c_char_p
lib.Pair_first_c_str.argtypes = [c_void_p]

lib.Pair_second_c_str.restype = c_char_p
lib.Pair_second_c_str.argtypes = [c_void_p]

lib.SBMLDoc_replicate_model_per_tissue.restype = c_void_p
lib.SBMLDoc_replicate_model_per_tissue.argtypes = [c_void_p, POINTER(c_char_p), c_size_t]

def _list_to_pointer(string_list: Iterable[str]):
    array_type = c_char_p*len(string_list)
    pointer = array_type(*(s.encode('utf-8') for s in string_list))
    return POINTER(c_char_p)(pointer)

def _iterate_proteins(proteins_ptr):
    it = lib.Proteins_iterator(proteins_ptr)
    try:
        while True:
            
            pair_ptr = lib.Proteins_iterator_next(it)
            
            if not pair_ptr:
                break
            key   = lib.Pair_first_c_str(pair_ptr).decode('utf-8')
            
            value = lib.Pair_second_c_str(pair_ptr).decode('utf-8')
            
            lib.Pair_delete(pair_ptr)
            
            yield key, value
    finally:
        lib.Proteins_delete_iterator(it)

class SBMLDoc:
    lib.SBMLDoc_new.restype = c_void_p
    lib.SBMLDoc_new.argtypes = [c_char_p]

    lib.SBMLDoc_number_of_kinetic_costant.restype = c_int
    lib.SBMLDoc_number_of_kinetic_costant.argtypes = [c_void_p]

    lib.SBMLDoc_set_kinetic_constants.restype = None
    lib.SBMLDoc_set_kinetic_constants.argtypes = [c_void_p, c_int, c_double]

    lib.SBMLDoc_save_converted_file.restype = c_bool
    lib.SBMLDoc_save_converted_file.argtypes = [c_void_p, c_char_p]

    lib.SBMLDoc_delete.restype = None
    lib.SBMLDoc_delete.argtypes = [c_void_p]

    lib.SBMLDoc_simulate.restype = None
    lib.SBMLDoc_simulate.argtypes = [c_void_p, c_char_p, c_double]

    lib.SBMLDoc_random_start_concentration.restype = None
    lib.SBMLDoc_random_start_concentration.argtypes = [c_void_p]

    lib.SBMLDoc_delete.restype = None
    lib.SBMLDoc_delete.argtypes = [c_void_p]

    lib.SBMLDoc_dump_proteins_data.restype = None
    lib.SBMLDoc_dump_proteins_data.argtypes = [c_void_p]

    lib.SBMLDoc_get_proteins_data.restype = c_void_p
    lib.SBMLDoc_get_proteins_data.argtypes = [c_void_p]

    lib.SBMLDoc_is_protein.restype = c_bool
    lib.SBMLDoc_is_protein.argtypes = [c_void_p, c_char_p]

    lib.SBMLDoc_random_protein_concentrations.restype = None
    lib.SBMLDoc_random_protein_concentrations.argtypes = [c_void_p]

    def __init__(self, file_path: str = None):
        if file_path is not None:
            self.obj = lib.SBMLDoc_new(file_path.encode('utf-8'))

    def number_of_kinetic_constants(self) -> int:
        return lib.SBMLDoc_number_of_kinetic_costant(self.obj)

    def set_kinetic_constants(self, id: int, value: float):
        lib.SBMLDoc_set_kinetic_constants(self.obj, id, value)

    def save_converted_file(self, output_path: str) -> bool:
        return lib.SBMLDoc_save_converted_file(self.obj, output_path.encode('utf-8'))
    
    def random_start_concentration(self):
        return lib.SBMLDoc_random_start_concentration(self.obj)

    def simulate(self, output_file = "simulation_results.csv", duration = 10.0):
        lib.SBMLDoc_simulate(self.obj, output_file.encode('utf-8'), duration)
    
    def dump_proteins_data(self):
        lib.SBMLDoc_dump_proteins_data(self.obj)
    
    def get_proteins_data(self) -> dict[str,str]:
        
        result = dict()
        proteins_ptr = lib.SBMLDoc_get_proteins_data(self.obj)
        
        for species, protein in _iterate_proteins(proteins_ptr):
            result[species] = protein
            
        return result
    
    def is_protein(self, specie: str) -> bool:
        return lib.SBMLDoc_is_protein(self.obj, specie.encode('utf-8'))
    
    def random_protein_concentrations(self):
        lib.SBMLDoc_random_protein_concentrations(self.obj)
    
    def replicate_model_per_tissue(self, tissues: Iterable[str]):
        obj = lib.SBMLDoc_replicate_model_per_tissue(self.obj, _list_to_pointer(tissues), len(tissues))
        result = SBMLDoc()
        result.obj = obj
        return result

    def __del__(self):
        if hasattr(self, 'obj') and self.obj:
            lib.SBMLDoc_delete(self.obj)
            self.obj = None


def replicate_model_per_tissue(file_path: str, tissues: list[str]):
    obj = lib.replicate_model_per_tissue(file_path.encode('utf-8'),_list_to_pointer(tissues), len(tissues))
    result = SBMLDoc()
    result.obj = obj
    return result