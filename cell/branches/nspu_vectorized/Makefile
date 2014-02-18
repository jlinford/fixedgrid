SPU_COMPLIER := spuxlc
PPU_COMPILER := ppuxlc
CC_OPT_LEVEL := -O0

########################################################################
# Subdirectories
########################################################################

DIRS := chem util spu

########################################################################
#                       Target
########################################################################

PROGRAM_ppu := fixedgrid

########################################################################
#                       Local Defines
########################################################################

IMPORTS := chem/lib_fixedgrid_chem.a util/lib_fixedgrid_util.a \
           spu/lib_fixedgrid_spu.a -lspe2 -lpthread -lm -lmisc

INCLUDE := -Iconfig -Ichem -Iutil

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

