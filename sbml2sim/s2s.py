from ctypes import cdll
from ctypes import c_void_p, c_char_p, c_bool, c_int, c_double
from ctypes import POINTER, c_size_t
import math
from typing import Generator, Iterable
from ctypes import c_uint

import sys

lib = cdll.LoadLibrary("build/libsbmlconverter.so")

lib.SpeciesToId_iterator.restype = c_void_p
lib.SpeciesToId_iterator.argtypes = [c_void_p]

lib.SpeciesToId_delete_iterator.restype = None
lib.SpeciesToId_delete_iterator.argtypes = [c_void_p]

lib.SpeciesToId_iterator_next.restype = c_void_p
lib.SpeciesToId_iterator_next.argtypes = [c_void_p]

lib.Pair_delete.restype = None
lib.Pair_delete.argtypes = [c_void_p]

lib.Pair_first_c_str.restype = c_char_p
lib.Pair_first_c_str.argtypes = [c_void_p]

lib.Pair_second_c_str.restype = c_char_p
lib.Pair_second_c_str.argtypes = [c_void_p]

lib.SBMLDoc_replicate_model_per_tissue.restype = c_void_p
lib.SBMLDoc_replicate_model_per_tissue.argtypes = [c_void_p, POINTER(c_char_p), c_size_t]

lib.SpeciesToId_iterator_end.restype = c_bool
lib.SpeciesToId_iterator_end.argtypes = [c_void_p, c_void_p]

lib.set_seed.restype = None
lib.set_seed.argtypes = [c_uint]

def set_seed(seed: int):
    lib.set_seed(seed)

def _list_to_pointer(string_list: Iterable[str]):
    array_type = c_char_p*len(string_list)
    pointer = array_type(*(s.encode('utf-8') for s in string_list))
    return POINTER(c_char_p)(pointer)

def _iterate_ids(SpeciesToId_ptr) -> Generator[tuple[str, str], None, None]:
    it = lib.SpeciesToId_iterator(SpeciesToId_ptr)
    try:
        while True:
            if lib.SpeciesToId_iterator_end(SpeciesToId_ptr, it):
                break
            pair_ptr = lib.SpeciesToId_iterator_next(it)
            if not pair_ptr:
                break
            key = lib.Pair_first_c_str(pair_ptr).decode('utf-8')
            value = lib.Pair_second_c_str(pair_ptr).decode('utf-8')
            lib.Pair_delete(pair_ptr)
            yield key, value
    finally:
        lib.SpeciesToId_delete_iterator(it)

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

    lib.SBMLDoc_get_proteins_data.restype = c_void_p
    lib.SBMLDoc_get_proteins_data.argtypes = [c_void_p]

    lib.SBMLDoc_is_protein.restype = c_bool
    lib.SBMLDoc_is_protein.argtypes = [c_void_p, c_char_p]

    lib.SBMLDoc_random_protein_concentrations.restype = None
    lib.SBMLDoc_random_protein_concentrations.argtypes = [c_void_p]
    
    
    lib.SBMLDoc_add_kinetic_laws_if_not_exists.restype = None
    lib.SBMLDoc_add_kinetic_laws_if_not_exists.argtypes = [c_void_p]

    lib.SBMLDoc_add_time_to_model.restype = None
    lib.SBMLDoc_add_time_to_model.argtypes = [c_void_p]

    lib.SBMLDoc_add_avg_calculations_for_all_species.restype = None
    lib.SBMLDoc_add_avg_calculations_for_all_species.argtypes = [c_void_p]

    lib.SBMLDoc_add_avg_calculation_for_all_proteins.restype = None
    lib.SBMLDoc_add_avg_calculation_for_all_proteins.argtypes = [c_void_p]

    lib.SBMLDoc_random_kinetic_costant_value.restype = None
    lib.SBMLDoc_random_kinetic_costant_value.argtypes = [c_void_p]

    lib.SBMDoc_small_compound_start_random_concentration.restype = None
    lib.SBMDoc_small_compound_start_random_concentration.argtypes = [c_void_p]

    lib.SBMLDoc_get_num_compartements.restype = c_uint
    lib.SBMLDoc_get_num_compartements.argtypes = [c_void_p]

    lib.SBMLDoc_set_volume_compartement.restype = None
    lib.SBMLDoc_set_volume_compartement.argtypes = [c_void_p, c_uint, c_double]

    lib.SBMLDoc_get_name_compartement.restype = c_char_p
    lib.SBMLDoc_get_name_compartement.argtypes = [c_void_p, c_uint]

    lib.SBMLDoc_get_volume_compartement.restype = c_double
    lib.SBMLDoc_get_volume_compartement.argtypes = [c_void_p, c_char_p]

    lib.SBMLDoc_get_compartement.restype = c_char_p
    lib.SBMLDoc_get_compartement.argtypes = [c_void_p, c_char_p]
    
    lib.SBMLDoc_set_initial_concentration.restype = None
    lib.SBMLDoc_set_initial_concentration.argtypes = [c_void_p, c_char_p, c_double]

    lib.SBMLDoc_input_start_random_concentration.restype = None
    lib.SBMLDoc_input_start_random_concentration.argtypes = [c_void_p]

    lib.SBMLDoc_set_zero_output_costant.restype = None
    lib.SBMLDoc_set_zero_output_costant.argtypes = [c_void_p]
    
    lib.SBMLDoc_get_compounds_data.restype = c_void_p
    lib.SBMLDoc_get_compounds_data.argtypes = [c_void_p]
    
    lib.SBMLDoc_set_parameter.restype = None
    lib.SBMLDoc_set_parameter.argtypes = [c_void_p, c_char_p, c_double]

    lib.SBMLDoc_assigment_rule_for_inputs.restype = None
    lib.SBMLDoc_assigment_rule_for_inputs.argtypes = [c_void_p]
    
    lib.SBMLDoc_remove_all_assigment_rules.restype = None
    lib.SBMLDoc_remove_all_assigment_rules.argtypes = [c_void_p]

    lib.SBMLDoc_convert_to_sbml_string.restype = c_void_p
    lib.SBMLDoc_convert_to_sbml_string.argtypes = [c_void_p]

    lib.deallocate_string.restype = None
    lib.deallocate_string.argtypes = [c_void_p]

    lib.SBMLDoc_get_kinetic_constants.restype = c_void_p
    lib.SBMLDoc_get_kinetic_constants.argtypes = [c_void_p]

    lib.SBMLDoc_get_output_constants.restype = c_void_p
    lib.SBMLDoc_get_output_constants.argtypes = [c_void_p]

    lib.SBMLDoc_delete_string_vector.restype = None
    lib.SBMLDoc_delete_string_vector.argtypes = [c_void_p]

    lib.SBMLDoc_string_vector_size.restype = c_size_t
    lib.SBMLDoc_string_vector_size.argtypes = [c_void_p]

    lib.SBMLDoc_string_vector_get.restype = c_char_p
    lib.SBMLDoc_string_vector_get.argtypes = [c_void_p, c_size_t]
    
    lib.SBMLDoc_is_output.restype = c_bool
    lib.SBMLDoc_is_output.argtypes = [c_void_p, c_char_p]

    lib.SBMLDoc_set_outputs_constants.restype = None
    lib.SBMLDoc_set_outputs_constants.argtypes = [c_void_p]

    lib.SBMLDoc_set_outputs_variable.restype = None
    lib.SBMLDoc_set_outputs_variable.argtypes = [c_void_p]

    def set_outputs_constants(self):
        lib.SBMLDoc_set_outputs_constants(self.obj)

    def set_outputs_variable(self):
        lib.SBMLDoc_set_outputs_variable(self.obj)

    def is_output(self, species_id: str) -> bool:
        return lib.SBMLDoc_is_output(self.obj, species_id.encode('utf-8'))

    def convert_to_string(self) -> str:
        sys.stdout.flush()
        sbml_ptr = lib.SBMLDoc_convert_to_sbml_string(self.obj)
        sbml_c_char_p = c_char_p(sbml_ptr)
        sbml_python_str = sbml_c_char_p.value.decode('utf-8')
        sbml_python_str = str(sbml_python_str)  # clone the string before deallocation
        lib.deallocate_string(sbml_ptr)
        return sbml_python_str

    def remove_all_assigment_rules(self):
        lib.SBMLDoc_remove_all_assigment_rules(self.obj)
    
    def set_parameter(self, id_parameter: str, value: float):
        lib.SBMLDoc_set_parameter(self.obj, id_parameter.encode('utf-8'), c_double(value))

    def assigment_rule_for_inputs(self):
        lib.SBMLDoc_assigment_rule_for_inputs(self.obj)
    
    def get_compounds_data(self) -> dict[str, str]:
        result = dict()
        compounds_ptr = lib.SBMLDoc_get_compounds_data(self.obj)
        for species, compound in _iterate_ids(compounds_ptr):
            result[species] = compound
        return result

    def set_zero_output_costant(self):
        lib.SBMLDoc_set_zero_output_costant(self.obj)
        
    def input_start_random_concentration(self):
        lib.SBMLDoc_input_start_random_concentration(self.obj)

    def set_initial_concentration(self, species_id: str, value: float):
        lib.SBMLDoc_set_initial_concentration(self.obj, species_id.encode('utf-8'), c_double(value))

    def get_volume_compartement(self, id: str) -> float:
        return lib.SBMLDoc_get_volume_compartement(self.obj, id.encode('utf-8'))

    def get_compartement(self, species_id: str) -> str:
        ptr = lib.SBMLDoc_get_compartement(self.obj, species_id.encode('utf-8'))
        return ptr.decode('utf-8')

    def get_num_compartements(self) -> int:
        return lib.SBMLDoc_get_num_compartements(self.obj)
    
    # TODO: id must be a string in the future
    def set_volume_compartement(self, id: int, volume: float):
        lib.SBMLDoc_set_volume_compartement(self.obj, c_uint(id), c_double(volume))
    
    def get_name_compartement(self, id: int) -> str:
        name_ptr = lib.SBMLDoc_get_name_compartement(self.obj, c_uint(id))
        return name_ptr.decode('utf-8') if name_ptr else ""

    def __init__(self, file_path: str = None):
        if file_path is not None:
            self.obj = lib.SBMLDoc_new(file_path.encode('utf-8'))
            
    def add_time_to_model(self):
        lib.SBMLDoc_add_time_to_model(self.obj)

    def add_avg_calculations_for_all_species(self):
        lib.SBMLDoc_add_avg_calculations_for_all_species(self.obj)

    def add_avg_calculation_for_all_proteins(self):
        lib.SBMLDoc_add_avg_calculation_for_all_proteins(self.obj)
        
    def add_kinetic_laws_if_not_exists(self):
        lib.SBMLDoc_add_kinetic_laws_if_not_exists(self.obj)

    def random_kinetic_costant_value(self):
        lib.SBMLDoc_random_kinetic_costant_value(self.obj)

    def small_compound_start_random_concentration(self):
        lib.SBMDoc_small_compound_start_random_concentration(self.obj)

    def number_of_kinetic_constants(self) -> int:
        return lib.SBMLDoc_number_of_kinetic_costant(self.obj)

    def set_kinetic_constants(self, id: int, value: float):
        lib.SBMLDoc_set_kinetic_constants(self.obj, id, value)

    def save_converted_file(self, output_path: str) -> bool:
        return lib.SBMLDoc_save_converted_file(self.obj, output_path.encode('utf-8'))
    
    def random_start_concentration(self):
        return lib.SBMLDoc_random_start_concentration(self.obj)

    def simulate(self, output_file = "simulation_results.csv", duration = 10.0):
        print("[DEPRECATED] SBMLDoc.simulate")
        exit(1)
    
    def get_proteins_data(self) -> dict[str,str]:
        
        result = dict()
        proteins_ptr = lib.SBMLDoc_get_proteins_data(self.obj)
        
        for species, protein in _iterate_ids(proteins_ptr):
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
    
    def get_kinetic_constants(self) -> list[str]:
        result = []
        ptr = lib.SBMLDoc_get_kinetic_constants(self.obj)
        size = lib.SBMLDoc_string_vector_size(ptr)
        for i in range(size):
            s = lib.SBMLDoc_string_vector_get(ptr, i)
            if s:
                result.append(str(s.decode('utf-8')))
        lib.SBMLDoc_delete_string_vector(ptr)
        return result
    
    def get_output_constants(self) -> list[str]:
        result = []
        ptr = lib.SBMLDoc_get_output_constants(self.obj)
        size = lib.SBMLDoc_string_vector_size(ptr)
        for i in range(size):
            s = lib.SBMLDoc_string_vector_get(ptr, i)
            if s:
                result.append(str(s.decode('utf-8')))
        lib.SBMLDoc_delete_string_vector(ptr)
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

