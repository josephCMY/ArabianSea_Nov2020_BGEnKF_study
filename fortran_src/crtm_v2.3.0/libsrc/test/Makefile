# libsrc/test/Makefile.  Generated from Makefile.in by configure.

# check target makefile template

# Package-specific substitution variables
tarname         = crtm
endian_platform = little
endian_io       = little

# The file definitions.
EXE_FILE=check_crtm
SRC_FILE=$(EXE_FILE).f90
OBJ_FILE=${SRC_FILE:.f90=.o}

# Tool-specific substitution variables
FC       = ifort
FCFLAGS  = -I.. -assume byterecl -O3 -fp-model source -e08 -free -assume byterecl,realloc_lhs -c
FPPFLAGS = -DLITTLE_ENDIAN 
LDFLAGS  = -L.. 
LIBS     = -l$(tarname) 

# Colours for test result output
GREEN=\033[1;32m
RED=\033[1;31m
NC=\033[0m # No Color

# The targets
all: $(EXE_FILE)

check: $(EXE_FILE)
	@-rm -f .signal >/dev/null 2>&1
	@echo "========================================"; \
	echo "  Running test: $(endian_io)-endian I/O on $(endian_platform)-endian platform..."; \
        echo; \
	./$(EXE_FILE); \
	echo; \
	if [ -f .signal ]; then \
	  echo -e "$(GREEN)"; echo "  TEST SUCCESSFUL!"; echo -e "$(NC)"; \
	  rm .signal; \
	else \
	  echo -e "$(RED)"  ; echo "  TEST FAILED!"    ; echo -e "$(NC)"; \
	fi; \
	echo  "========================================"

clean:
	-rm -fr $(OBJ_FILE) $(EXE_FILE) gmon.out *.output *.bin >/dev/null 2>&1

distclean:
	-rm Makefile >/dev/null 2>&1

update:
	@update() \
	{ files=`find . -maxdepth 1 -name "$$1" -print`; \
	  if [ -n "$$files" ]; then \
	    mv $$files results; \
	  else \
	    echo "No $$1 files to update."; \
	  fi \
	}; \
	echo; echo; \
	update "*.output"; update "*.bin"

$(EXE_FILE): $(OBJ_FILE)
	$(FC) $(LDFLAGS) $(OBJ_FILE) -o $(EXE_FILE) $(LIBS)

# Specify targets that do not generate filesystem objects
.PHONY: all clean check update distclean

# Compilation rule include files
include ../make.rules
