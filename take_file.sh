#/bin/bash

if [ -z "$1" ]; then
  echo "Error: No filename provided."
  exit 1
fi

docker cp sbml2sim:/app/$1 .