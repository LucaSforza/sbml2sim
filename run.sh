#/bin/sh

xhost +local:docker

docker run -it -e DISPLAY=$DISPLAY --net=host --name sbml2sim sbml2sim

xhost -local:docker