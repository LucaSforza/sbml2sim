#/bin/bash

if [ -z "$1" ]; then
  echo "Error: No filename provided."
  exit 1
fi

echo "[INFO] coping $1"

docker cp sbml2sim:/app/$1 .