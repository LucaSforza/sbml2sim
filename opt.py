import nevergrad as ng

import sbmlconverter as sc

sbml = sc.SBMLDoc("sbmls/R-HSA-1643713.sbml")
sbml.random_start_concentration()
sbml.simulate()