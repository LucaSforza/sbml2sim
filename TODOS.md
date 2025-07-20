# 21/07/2025 - 31/07/2025

- Use randomized initial concentrations for each restart; run the simulation multiple times during evaluation of loss function.
- Run the project on the cluster and parallelize the computations.
- Check if the system is stable, otherwise penalize
- Complete the class to execute all possible combinations.
- write an email to the prof with all the results

# 20/07/2025

- make the project presentable
- Improve the simulation script to allow plotting of individual constrained species, including the true value, as well as a plot showing all constrained species together.
- Support multiple kinetic constant files as input.

# 19/07/2025
- resolve TODOs.

# 13/06/2025

- Run the project on the cluster and parallelize the computations.

# 12/06/2025
- Check if the system is stable, otherwise penalize
- Enable simulation of a model using an external file for kinetic constants.

# 11/06/2025

## TODOs:
- Enable simulation of a model using an external file for kinetic constants.
- The constants for a given enzyme are the same across all reactions in which it participates.
- Check if the system is stable, otherwise penalize
- Run the project on the cluster and parallelize the computations.

## OPTIONALS_
- The constants for a given enzyme are the same across all reactions in which it participates.

# 9/06/2025

## TODOs:
- the same from yesterday
- create binding for simulator.
- run simulation in parallel with different costant
- run using Nevergrad

# 8/06/2025

## TODOs:
- Create a simulatable file that includes proteomics data.
- Enable simulation of a model using an external file for kinetic constants.
- Simulate a model where parameters can be specified via a boolean array. (reaction are only fast or slow)
- The constants for a given enzyme are the same across all reactions in which it participates.
- implement prof idea:
    - Given a parametric model and a boolean array (reactions are either fast or slow), find all arrays that satisfy the following constraints:
        - The concentrations of the species must be ordered according to the experimental data (if available).
        - The system is stable, meaning the concentrations remain within a certain range, specifically between 0 and the concentration of water.
- Run the project on the cluster and parallelize the computations.


# 7/06/2025

## TODOs:
- create scripts for simulating an SBML model
- implement prof idea: TODO: write prof idea

# 26/06/2025

## TODOs
- Divide the optimitation in 2 steps:
    1. Search of the kinetic constants
    2. Search of the output constants

## Optional TODOS
- create unit measure

# 26/06/2025

## TODOs
- Begin black-box optimization process
    - use roadrunner only on python, the C++ library now only convert the SBML into a simulable one
    - take the libSBML code and put in the repo, create a decent build system that build the library
    - simulate the model in python
    - create the utility function
        - if floating points errors give a very hight penalty
        - if the system is unstable give a hight penalty
        - minimize the quadratic error of the outputs and internal species

## Optional TODOS
- create unit measure

# 24/06/2025

## TODOs
- Begin black-box optimization process
- Retrieve data for small compounds.
- Implement prof idea
    - for each input x set the value of x: x(t) = A + B sin(omega t + phi)
    - omega = 2 pi f where f in {10^-6,10^-5,...,10^-1}
    - phi = k pi 4 where k in {0, 1, 2, ..., 7}
    - B = k A/2 where k in {0,1,...,4}

## Things Discovered

- If a protein has hasVersion it means that is a modified species

## Optional TODOS

- create unit measure

# 23/06/2025

## TODOs

- Set the initial concentration of proteins using data obtained from the APIs.
- Retrieve data for small compounds.
- Convert iBAQ into mol/L

## Optional TODOS

- create unit measure

# 22/06/2025

## TODOs

- Set the initial concentration of proteins using data obtained from the APIs.
- Retrieve data for small compounds.
- Convert iBAQ into mol/L
- Find volumes for the compartements

## Optional TODOS

- create unit measure

# 21/06/2025

## TODOs

## Things Discovered

- Protein complexes should not be expanded, but should all be set to 0 before the simulation

## TODOs

- Replicate the biological model only for tissues with complete available data. (easy)
- Set the initial concentration of proteins using data obtained from the APIs.
- Retrieve data for small compounds.
- Isolate the inputs and make them constant. (easy)
- Convert iBAQ into mol/L
- Find volumes for the compartements

# Optional TODOs

- create unit measure

# 19/06/2025

## TODOs

- Expand all DefinedSets and clone all reactions; eliminate all DefinedSets as species.
- Expand all Protein Complexes.
- Replicate the biological model only for tissues with complete available data.
- Delete all drugs and drug-related reactions.
- Set the initial concentration of proteins using data obtained from the APIs.
- Retrieve data for small compounds.
- Isolate the inputs and make them constant.
- Do not use the convenience rate law anymore.

## Questions

- The outputs are species that appear only as products and never as reactants, but the simulator accumulates them. Should I add a degradation reaction? If so, what rate should it have? A new constant to be determined, but should it be different for each tissue or the same for all? Should this constant be proportional to the concentration or not?

- If the small compounds have different units of measurement compared to proteins, how should I handle this? Should I normalize them? If so, how?
