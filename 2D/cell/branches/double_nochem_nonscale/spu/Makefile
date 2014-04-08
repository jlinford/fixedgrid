SPU_COMPILER := spuxlc
PPU_COMPLIER := ppuxlc

#######################################################################
# Target
########################################################################

PROGRAMS_spu    := fixedgrid_spu 

# created embedded library
LIBRARY_embed := lib_fixedgrid_spu.a

IMPORTS := discretize.o timer.o

INCLUDE := -I../config -I../util -I../chem -I..

########################################################################
# make.footer
########################################################################
# make.footer is in the top of the SDK

ifdef CELL_TOP
 include $(CELL_TOP)/make.footer
else
 include /opt/ibm/cell-sdk/prototype/make.footer
endif

