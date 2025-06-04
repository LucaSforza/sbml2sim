all: libsbmlconverter.so

INCLUDE_DIRS= -I./include/
LIB_DIRS= -L./lib64

CXX =g++
CXXFLAGS= $(INCLUDE_DIRS) -fPIC
LDFLAGS= $(LIB_DIRS) -lroadrunner -lsbml


folder:
	mkdir -p build

libsbmlconverter.so: folder
	$(CXX) $(CXXFLAGS) -shared -o build/$@ src/bindings.cpp $(LDFLAGS)