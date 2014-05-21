
TIMERS_ENABLED = 1

ifndef TAU_MAKEFILE
  CXX = g++
else
  CXX = tau_cxx.sh -optTauSelectFile=select.tau
endif

REPORT_CXXFLAGS = -Wall
OPT_CXXFLAGS = -O3 -g -fopenmp

LD = $(CXX)
REPORT_LDFLAGS = -Wall
OPT_LDFLAGS = -O3 -g -fopenmp

LIBS = -lm

################################################################################

CXXFLAGS = $(OPT_CXXFLAGS) $(MIC_CXXFLAGS) $(REPORT_CXXFLAGS)
LDFLAGS = $(OPT_LDFLAGS) $(MIC_LDFLAGS) $(REPORT_LDFLAGS)

CPPFLAGS = -DTIMERS_ENABLED=$(TIMERS_ENABLED) -I.

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
	$(RM) $(EXE) $(OBJ) *MIC.o

