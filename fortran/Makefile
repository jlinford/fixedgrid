
TIMERS_ENABLED = 0

ifndef TAU_MAKEFILE
  F90 = ifort
else
  F90 = tau_f90.sh
endif

ifdef MIC_OFFLOAD
  MIC_F90FLAGS = -xhost
  MIC_LDFLAGS =
else
  MIC_F90FLAGS = -mmic
  MIC_LDFLAGS = -mmic
endif

REPORT_F90FLAGS = -opt-report-phase=offload -vec-report5
OPT_F90FLAGS = -O3 -g -openmp

LD = $(F90)
REPORT_LDFLAGS = -opt-report-phase=offload -vec-report5
OPT_LDFLAGS = -O3 -g -openmp

LIBS = -lm

################################################################################

F90FLAGS = $(OPT_F90FLAGS) $(MIC_F90FLAGS) $(REPORT_F90FLAGS)
LDFLAGS = $(OPT_LDFLAGS) $(MIC_LDFLAGS) $(REPORT_LDFLAGS)

CPPFLAGS = -DTIMERS_ENABLED=$(TIMERS_ENABLED) -I.

RM = /bin/rm -f

EXE = fixedgrid

SRC = $(wildcard *.F90)
OBJ = $(SRC:.F90=.o)

.PHONY: all clean

all: $(EXE)

$(EXE): $(OBJ)
	$(LD) $(LDFLAGS) $(OBJ) -o $(EXE) $(LIBS)

%.o: %.F90
	$(F90) $(CPPFLAGS) $(F90FLAGS) -c $< -o $@

clean: 
	$(RM) $(EXE) $(OBJ) *MIC.o *.mod

