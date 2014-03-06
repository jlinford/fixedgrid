TARGET = offload
TIMERS_ENABLED = 1
USE_BLOCKED_DISCRETIZE = 1

ifndef TAU_MAKEFILE
  CXX = icpc
else
  CXX = tau_cxx.sh -optTauSelectFile=select.tau
endif

ifeq ($(TARGET),host)
  MIC_CXXFLAGS = -xhost -no-offload
  MIC_LDFLAGS = 
else ifeq ($(TARGET),device)
  MIC_CXXFLAGS = -mmic
  MIC_LDFLAGS = -mmic
else ifeq ($(TARGET),offload)
  MIC_CXXFLAGS = -xhost
  MIC_LDFLAGS = 
else
  $(error Unknown TARGET)
endif

REPORT_CXXFLAGS = -opt-report-phase=offload -vec-report5
OPT_CXXFLAGS = -O3 -g -openmp

LD = $(CXX)
REPORT_LDFLAGS = -opt-report-phase=offload -vec-report5
OPT_LDFLAGS = -O3 -g -openmp

LIBS = -lm

################################################################################

CXXFLAGS = $(OPT_CXXFLAGS) $(MIC_CXXFLAGS) $(REPORT_CXXFLAGS)
LDFLAGS = $(OPT_LDFLAGS) $(MIC_LDFLAGS) $(REPORT_LDFLAGS)

CPPFLAGS = -DTIMERS_ENABLED=$(TIMERS_ENABLED) \
           -DUSE_BLOCKED_DISCRETIZE=$(USE_BLOCKED_DISCRETIZE) \
	   -I.

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
	$(RM) $(EXE) $(OBJ) *MIC.o

