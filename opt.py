import nevergrad as ng

import sbmlconverter as sc

sbml = sc.SBMLDoc("sbmls/R-HSA-1643713.sbml")
sbml.random_start_concentration()
sbml.simulate()

# def square(x):
#     return sum((x - .5)**2)

# optimizer = ng.optimizers.Wiz(parametrization=2,budget=100)
# recommendation = optimizer.minimize(square)
# print(recommendation.value)  # recommended value