#######################################################################
# Target
########################################################################

PROGRAM_spu		:= fixedgrid_spu 

LIBRARY_embed		:= lib_fixedgrid_spu.a

LINKERSCRIPT	= linker.script
H_BUILT_HERE	= $(LINKERSCRIPT)

VPATH			:= ../common saprc99

CFLAGS_gcc		:= 
CFLAGS_xlc		:= -qcompact
LDFLAGS			:= -Wl,-T,linker.script -lsimdmath -lm
#LDFLAGS			:= -lsimdmath -lm

INCLUDE			:= -I../config

########################################################################
# make.footer
########################################################################
# make.footer is in the top of the SDK

# Override default opt level for SPU code
CC_OPT_LEVEL := -O5

include ../Makefile.defs

ifdef GNU_ROOT
$(LINKERSCRIPT): linker.script.ed
	$(CP) $(GNU_ROOT)/spu/lib/ldscripts/elf32_spu.x linker.script
	ed linker.script < linker.script.ed > /dev/null 2>&1
else
$(LINKERSCRIPT): linker.script.ed
	$(CP) /usr/spu/lib/ldscripts/elf32_spu.x linker.script
	ed linker.script < linker.script.ed > /dev/null 2>&1
endif

cleanlocal:
	@$(RM) $(LINKERSCRIPT)
