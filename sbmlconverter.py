from ctypes import cdll
from ctypes import c_void_p, c_char_p, c_bool, c_int, c_double

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

class SBMLDoc:
    def __init__(self, file_path: str, all_convenience_law: bool = False):
        self.obj = lib.SBMLDoc_new(file_path.encode('utf-8'), all_convenience_law)

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

    def __del__(self):
        if hasattr(self, 'obj') and self.obj:
            lib.SBMLDoc_delete(self.obj)
            self.obj = None