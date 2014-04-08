SPU_COMPLIER := spuxlc
PPU_COMPILER := ppuxlc

########################################################################
# Subdirectories
########################################################################

DIRS := spu 

########################################################################
#                       Target
########################################################################

PROGRAM_ppu := fixedgrid

########################################################################
#                       Local Defines
########################################################################

IMPORTS := spu/lib_fixedgrid_spu.a -lspe2 -lpthread

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

