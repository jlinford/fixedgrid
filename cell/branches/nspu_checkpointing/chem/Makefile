CC_OPT_LEVEL := -O0
#######################################################################
# Target
########################################################################

# created embedded library
LIBRARY := lib_fixedgrid_chem.a

IMPORTS := saprc99_Function.o saprc99_Initialize.o saprc99_Integrator.o \
           saprc99_Jacobian.o saprc99_JacobianSP.o saprc99_LinearAlgebra.o \
           saprc99_Main.o saprc99_Rates.o saprc99_Util.o saprc99_Monitor.o \
           -lm

INCLUDE := -I../config

########################################################################
# make.footer
########################################################################
# make.footer is in the top of the SDK

ifdef CELL_TOP
 include $(CELL_TOP)/make.footer
else
 include /opt/ibm/cell-sdk/prototype/make.footer
endif

