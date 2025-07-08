#!/bin/bash

if [ $# -eq 0 ]; then
    echo "[ERROR] missing argument"
    echo "[INFO]: $0 <sbml file path>"
    exit 1
fi

file_path="$1"
output_path=converted_$(basename $file_path)

PYTHON=/venv-sbml2sim/bin/python3

docker run -it -e \
    --net=host \
    --name sbml2sim sbml2sim \
    $PYTHON sbml2sim/convert.py $file_path $output_path

./take_file.sh $output_path

docker rm -f sbml2sim 2>/dev/null || true