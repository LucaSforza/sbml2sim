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
