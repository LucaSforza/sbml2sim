
#!/bin/bash

if [ $# -eq 0 ]; then
    echo "[ERROR] missing argument"
    echo "[INFO]: $0 <sbml file path>"
    exit 1
fi

file_path="$1"

PYTHON=/venv-sbml2sim/bin/python3

docker run -it -e \
    DISPLAY=$DISPLAY --net=host \
    --name sbml2sim sbml2sim \
    $PYTHON sbml2sim $file_path breast_cancer_cell

PATHWAY_ID=$(basename "$file_path" .sbml)
./take_file.sh sbmls/$PATHWAY_ID-real-tissues-modified.sbml
./take_file.sh parameters.json
./take_file.sh simulation.csv
./take_file.sh simulation.png

docker rm -f sbml2sim 2>/dev/null || true