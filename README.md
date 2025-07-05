# sbml2sim

sbml2sim è uno strumento che permette di convertire e simulare modelli in formato SBML (Systems Biology Markup Language) di reactome, facilitando l’identificazione delle costanti cinetiche delle reazioni tramite tecniche di ottimizzazione black box.

# utilizzo


Per fare la build tramite docker:

```bash
./build.sh
```

Per verificare che tutto funzioni:

```bash
./test.sh
```

Per convertire un file SBML di reactome e trovare le costanti cinetiche:

```bash
./opt.sh <SBML file path>
```

Per esempio:

```bash
./opt.sh sbmls/R-HSA-391251.sbml
```
