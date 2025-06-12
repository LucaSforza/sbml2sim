all: libsbmlconverter.so

INCLUDE_DIRS= -I./include/
LIB_DIRS= -L./lib64

CXX =g++
CXXFLAGS= $(INCLUDE_DIRS) -fPIC
LDFLAGS= $(LIB_DIRS) -lroadrunner -lsbml


folder:
	mkdir -p build

libsbmlconverter.so: folder
	$(CXX) $(CXXFLAGS) -ggdb -shared -o build/$@ src/bindings.cpp $(LDFLAGS)

main: libsbmlconverter.so
	$(CXX) $(CXXFLAGS) -ggdb -o build/$@ src/main.cpp build/libsbmlconverter.so $(LDFLAGS)