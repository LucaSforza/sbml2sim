from ctypes import cdll
from ctypes import c_void_p, c_char_p, c_bool, c_int, c_double
from ctypes import POINTER, c_size_t

lib = cdll.LoadLibrary("build/libsbmlconverter.so")

lib.SBMLDoc_new.restype = c_void_p
lib.SBMLDoc_new.argtypes = [c_char_p, c_bool]

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

lib.SBMLDoc_dump_genes_data.restype = None
lib.SBMLDoc_dump_genes_data.argtypes = [c_void_p]

lib.SBMLDoc_get_genes_data.restype = c_void_p
lib.SBMLDoc_get_genes_data.argtypes = [c_void_p]

lib.Genes_proteins_iterator.restype = c_void_p
lib.Genes_proteins_iterator.argtypes = [c_void_p]

lib.Genes_delete_proteins_iterator.restype = None
lib.Genes_delete_proteins_iterator.argtypes = [c_void_p]

lib.Genes_proteins_iterator_next.restype = c_void_p
lib.Genes_proteins_iterator_next.argtypes = [c_void_p]

lib.Pair_delete.restype = None
lib.Pair_delete.argtypes = [c_void_p]

lib.Pair_first_c_str.restype = c_char_p
lib.Pair_first_c_str.argtypes = [c_void_p]

lib.Pair_second_as_cstr_array.restype = POINTER(c_char_p)
lib.Pair_second_as_cstr_array.argtypes = [c_void_p, POINTER(c_size_t)]

lib.Pair_delete_cstr_array.restype = None
lib.Pair_delete_cstr_array.argtypes = [POINTER(c_char_p)]

lib.SBMLDoc_is_protein.restype = c_bool
lib.SBMLDoc_is_protein.argtypes = [c_void_p, c_char_p]

lib.SBMLDoc_random_protein_concentrations.restype = None
lib.SBMLDoc_random_protein_concentrations.argtypes = [c_void_p]

def _iterate_genes(genes_ptr):
    it = lib.Genes_proteins_iterator(genes_ptr)
    try:
        while True:
            pair_ptr = lib.Genes_proteins_iterator_next(it)
            if not pair_ptr:
                break
            key = lib.Pair_first_c_str(pair_ptr).decode('utf-8')
            size = c_size_t()
            arr = lib.Pair_second_as_cstr_array(pair_ptr, size)
            values = [arr[i].decode('utf-8') for i in range(size.value)]
            lib.Pair_delete_cstr_array(arr)
            lib.Pair_delete(pair_ptr)
            yield key, values
    finally:
        lib.Genes_delete_proteins_iterator(it)

class SBMLDoc:
    def __init__(self, file_path: str, all_convenience_law: bool = False, avg_for_only_proteins: bool = False):
        flags = (1 if all_convenience_law else 0) | (2 if avg_for_only_proteins else 0)
        self.obj = lib.SBMLDoc_new(file_path.encode('utf-8'), flags)

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
    
    def dump_genes_data(self):
        lib.SBMLDoc_dump_genes_data(self.obj)
    
    def get_genes_data(self) -> dict[str,list[str]]:
        result = dict()
        genes_ptr = lib.SBMLDoc_get_genes_data(self.obj)
        
        for protein, genes in _iterate_genes(genes_ptr):
            result[protein] = genes
            
        return result
    
    def is_protein(self, specie: str) -> bool:
        return lib.SBMLDoc_is_protein(self.obj, specie.encode('utf-8'))
    
    def random_protein_concentrations(self):
        lib.SBMLDoc_random_protein_concentrations(self.obj)

    def __del__(self):
        if hasattr(self, 'obj') and self.obj:
            lib.SBMLDoc_delete(self.obj)
            self.obj = None