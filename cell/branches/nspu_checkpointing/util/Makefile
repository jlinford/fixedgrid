PPU_COMPILER := ppuxlc
CC_OPT_LEVEL := -O0

########################################################################
#                       Target
########################################################################

LIBRARY_ppu := lib_fixedgrid_util.a

IMPORTS := -lm

########################################################################
#                       Local Defines
########################################################################

INCLUDE := -I.. -I../config -I../chem

# imports the embedded simple_spu library
# allows consolidation of spu program into ppe binary
########################################################################
# make.footer
########################################################################
# make.footer is in the top of the SDK

ifdef CELL_TOP
 include $(CELL_TOP)/make.footer
else
 include /opt/ibm/cell-sdk/prototype/make.footer
endif

