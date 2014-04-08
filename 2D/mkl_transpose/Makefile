
TIMERS_ENABLED = 1

ifndef TAU_MAKEFILE
  CXX = icpc
else
  CXX = tau_cxx.sh -optTauSelectFile=select.tau
endif

REPORT_CXXFLAGS = -vec-report2
OPT_CXXFLAGS = -O3 -g -openmp -mkl -xhost

LD = $(CXX)
REPORT_LDFLAGS = -vec-report2
OPT_LDFLAGS = -O3 -g -openmp -mkl

LIBS = -lm

################################################################################

CXXFLAGS = $(OPT_CXXFLAGS) $(REPORT_CXXFLAGS)
LDFLAGS = $(OPT_LDFLAGS) $(REPORT_LDFLAGS)

CPPFLAGS = -DTIMERS_ENABLED=$(TIMERS_ENABLED) -I.

RM = /bin/rm -f

EXE = fixedgrid

SRC = $(wildcard *.cpp)
HDR = $(wildcard *.hpp)
OBJ = $(SRC:.cpp=.o)
MAKEFILE = Makefile

.PHONY: all clean

all: $(EXE)

$(EXE): $(OBJ)
	$(LD) $(LDFLAGS) $(OBJ) -o $(EXE) $(LIBS)

%.o: %.cpp $(HDR) $(MAKEFILE)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -c $< -o $@

clean: 
	$(RM) $(EXE) $(OBJ) 

