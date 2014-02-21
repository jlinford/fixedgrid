
MIC_OFFLOAD = 1
TIMERS_ENABLED = 0

ifndef TAU_MAKEFILE
  CXX = icpc
else
  CXX = tau_cxx.sh -optTauSelectFile=select.tau
endif

ifeq ($(MIC_OFFLOAD),1)
  MIC_CXXFLAGS = -xhost
  MIC_LDFLAGS =
else
  MIC_CXXFLAGS = -mmic
  MIC_LDFLAGS = -mmic
endif

REPORT_CXXFLAGS = -opt-report-phase=offload -vec-report2
OPT_CXXFLAGS = -O3 -g -openmp

LD = $(CXX)
REPORT_LDFLAGS = -opt-report-phase=offload -vec-report2
OPT_LDFLAGS = -O3 -g -openmp

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

