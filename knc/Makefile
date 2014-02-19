
ifndef TAU_MAKEFILE
CXX = g++
else
CXX = tau_cxx.sh -optTauSelectFile=select.tau
endif
CXXFLAGS = -O3 -g -fopenmp
CPPFLAGS = -I.

LD = $(CXX)
LDFLAGS = -fopenmp
LIBS = -lm

RM = /bin/rm -f

EXE = fixedgrid

SRC = $(wildcard *.cpp)
HDR = $(wildcard *.hpp)
OBJ = $(SRC:.cpp=.o)

.PHONY: all clean

all: $(EXE)

$(EXE): $(OBJ)
	$(LD) $(LDFLAGS) $(OBJ) -o $(EXE) $(LIBS)

%.o: %.cpp $(HDR)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -c $< -o $@

clean: 
	$(RM) $(EXE) $(OBJ)

