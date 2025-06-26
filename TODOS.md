# 24/06/2025

## TODOs
- Begin black-box optimization process
- Implement prof idea
    - for each input x set the value of x: x(t) = A + B sin(omega t + phi)
    - omega = 2 pi f where f in {10^-6,10^-5,...,10^-1}
    - phi = k pi 4 where k in {0, 1, 2, ..., 7}
    - B = k A/2 where k in {0,1,...,4}

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
