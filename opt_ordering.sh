#!/bin/bash

if [ $# -eq 0 ]; then
    echo "[ERROR] missing argument"
    echo "[INFO]: $0 <sbml file path> [--proteomics <proteomics.json>] [--workers N] [--tissue TISSUE] [--plot] [--output-file FILE]"
    exit 1
fi

file_path="$1"
shift

PYTHON=/venv-sbml2sim/bin/python3

docker run -it \
    -e DISPLAY=$DISPLAY --net=host \
    --name sbml2sim sbml2sim \
    $PYTHON sbml2sim/opt_ordering.py "$file_path" "$@"

docker rm -f sbml2sim 2>/dev/null || true