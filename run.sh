#/bin/sh

xhost +local:docker

docker run -it --rm -e DISPLAY=$DISPLAY --net=host --cap-add=SYS_PTRACE --security-opt seccomp=unconfined --ulimit core=-1 --name sbml2sim sbml2sim bash

xhost -local:docker