all: lib

build:
	mkdir -p build

lib: build src/sbml_converter.cpp
	g++ -shared -I./include/ -I/usr/include/raptor2 -I/usr/include/rasqal -L./lib64 src/sbml_converter.cpp -o build/libsbmlconverter.so -lroadrunner -lsbml -lrdf -fPIC