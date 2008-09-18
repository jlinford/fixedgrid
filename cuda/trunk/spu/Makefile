#######################################################################
# Target
########################################################################

PROGRAM_spu		:= fixedgrid_spu 

LIBRARY_embed		:= lib_fixedgrid_spu.a

VPATH			:= ../common saprc99

LDFLAGS			:= -lsimdmath -lm

INCLUDE			:= -I../config

########################################################################
# make.footer
########################################################################
# make.footer is in the top of the SDK

# Override default opt level for SPU code
CC_OPT_LEVEL := -O5

include ../Makefile.defs

