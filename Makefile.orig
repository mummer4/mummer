#-------------------------------------------------------------------------------
#  Top level makefile for MUMmer 3.0
#
#  Dependencies: '/bin/sh', 'g++', 'gcc', 'csh', 'perl', 'sed'
#
#  'make all' builds all MUMmer code and scripts in the base directory
#
#  'make check' checks for the existance of the MUMmer dependencies
#
#  'make clean' removes *.o *~ core* and executable files
#
#  'make dist' creates a gzipped tarfile of the MUMmer directories
#
#  'make install' alias for 'make all' (for compatibility)
#
#  'make kurtz' builds Stefan's mummer program in the base directory
#
#  'make scripts' builds the MUMmer scripts in the base directory
#
#  'make tigr' builds TIGR's code in the base directory
#
#  'make uninstall' alias for 'make clean' (for compatibility)
#
#-------------------------------------------------------------------------------
SHELL = /bin/sh
VERSION := 3.23


TOP_DIR     := $(CURDIR)
BIN_DIR     := $(TOP_DIR)
AUX_BIN_DIR := $(TOP_DIR)/aux_bin

DOC_DIR       := $(TOP_DIR)/docs
SCRIPT_DIR    := $(TOP_DIR)/scripts
TIGR_SRC_DIR  := $(TOP_DIR)/src/tigr
KURTZ_SRC_DIR := $(TOP_DIR)/src/kurtz

CC   := $(filter /%,$(shell /bin/sh -c 'type gcc'))
CXX  := $(filter /%,$(shell /bin/sh -c 'type g++'))
SED  := $(filter /%,$(shell /bin/sh -c 'type sed'))
CSH  := $(filter /%,$(shell /bin/sh -c 'type csh'))
PERL := $(filter /%,$(shell /bin/sh -c 'type perl'))
AR   := $(filter /%,$(shell /bin/sh -c 'type ar'))

CXXFLAGS = -O3
CFLAGS = -O3
LDFLAGS  =

FLATS = ACKNOWLEDGEMENTS COPYRIGHT INSTALL LICENSE Makefile README ChangeLog



#-- EXPORT THESE VARIABLES TO OTHER MAKEFILES
export BIN_DIR AUX_BIN_DIR CXX CC CFLAGS CXXFLAGS LDFLAGS




#-- PHONY rules --#
.PHONY: all check clean dist scripts uninstall


all: kurtz tigr scripts


check:
ifndef TOP_DIR
	@echo "ERROR: could not find working directory"
endif
ifndef CC
	@echo "ERROR: 'gcc' GNU C compiler not found"
endif
ifndef CXX
	@echo "ERROR: 'g++' GNU C++ compiler not found"
endif
ifndef SED
	@echo "ERROR: 'sed' StreamEDitor not found"
endif
ifndef CSH
	@echo "ERROR: 'csh' C-shell not found"
endif
ifndef PERL
	@echo "ERROR: 'perl' PERL not found"
endif
ifndef AR
	@echo "ERROR: 'ar' GNU archiver not found"
endif
	@echo "check complete"


clean:
	rm -f *~ core*
	cd $(KURTZ_SRC_DIR); $(MAKE) clean
	cd $(TIGR_SRC_DIR); $(MAKE) clean
	cd $(SCRIPT_DIR); $(MAKE) clean
	cd $(DOC_DIR); $(MAKE) clean


dist: DISTDIR = MUMmer$(VERSION)
dist:
	mkdir $(DISTDIR)
	cp -r aux_bin $(DISTDIR)
	cp -r docs $(DISTDIR)
	cp -r scripts $(DISTDIR)
	cp -r src $(DISTDIR)
	cp $(FLATS) $(DISTDIR)
	rm -rf `find $(DISTDIR) -name CVS`
	tar -cvf $(DISTDIR).tar $(DISTDIR)
	gzip $(DISTDIR).tar
	rm -rf $(DISTDIR)


install: all


kurtz:
	cd $(KURTZ_SRC_DIR); $(MAKE) mummer


scripts:
	cd $(SCRIPT_DIR); $(MAKE) all


tigr:
	cd $(TIGR_SRC_DIR); $(MAKE) all


uninstall: clean


#-- END OF MAKEFILE --#
