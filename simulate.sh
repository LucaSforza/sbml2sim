#!/bin/bash

if [ $# -eq 0 ]; then
    echo "[ERROR] missing argument"
    echo "[INFO]: $0 <sbml file path>"
    exit 1
fi

file_path="$1"
output_path=converted_$(basename $file_path)

PYTHON=/venv-sbml2sim/bin/python3

# TODO: - Enable simulation of a model using an external file for kinetic constants.

docker run -it -e \
    --net=host \
    --name sbml2sim sbml2sim \
    $PYTHON sbml2sim/simulate.py $file_path --plot

./take_file.sh simulation.csv
./take_file.sh simulation.png

docker rm -f sbml2sim 2>/dev/null || true