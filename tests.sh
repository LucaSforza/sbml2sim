#!/bin/sh

PATHWAY_ID=R-HSA-391251

docker run -it -e DISPLAY=$DISPLAY --net=host --name sbml2sim sbml2sim /venv-sbml2sim/bin/python3 sbml2sim/tests.py sbmls/$PATHWAY_ID.sbml breast_cancer_cell


./take_file.sh simulation.png
./take_file.sh simulation.csv
./take_file.sh sbmls/$PATHWAY_ID-real-tissues-modified.sbml
./take_file.sh proteomics.json
./take_file.sh sbmls/$PATHWAY_ID-modified.sbml

docker rm -f sbml2sim 2>/dev/null || true
