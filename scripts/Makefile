#-- Imported variables from top level makefile
# BIN_DIR AUX_BIN_DIR CXX CC CFLAGS CXXFLAGS LDFLAGS

ifndef BIN_DIR
BIN_DIR := $(CURDIR)
endif
ifndef AUX_BIN_DIR
AUX_BIN_DIR := $(CURDIR)
endif
ifndef SCRIPT_DIR
SCRIPT_DIR := $(CURDIR)
endif

SCRIPT_DIR := $(CURDIR)
SED := $(filter /%,$(shell /bin/sh -c 'type sed'))
CSH := $(filter /%,$(shell /bin/sh -c 'type csh'))
PERL := $(filter /%,$(shell /bin/sh -c 'type perl'))
VPATH := $(BIN_DIR)

ALL := exact-tandems mapview mummerplot nucmer promer \
       run-mummer1 run-mummer3 nucmer2xfig dnadiff

#-- PHONY rules --#
.PHONY: all clean


all: $(ALL)


clean:
	rm -f *~
	cd $(BIN_DIR); rm -f $(ALL)




#-- not so PHONY rules --#
exact-tandems: exact-tandems.csh
	$(SED)  -e 's?__CSH_PATH?$(CSH)?g' \
		-e 's?__BIN_DIR?$(BIN_DIR)?g' \
		-e 's?__SCRIPT_DIR?$(SCRIPT_DIR)?g' \
		exact-tandems.csh > $(BIN_DIR)/exact-tandems
	chmod 755 $(BIN_DIR)/exact-tandems

mapview: mapview.pl
	$(SED)  -e 's?__PERL_PATH?$(PERL)?g' \
		-e 's?__SCRIPT_DIR?$(SCRIPT_DIR)?g' \
                mapview.pl > $(BIN_DIR)/mapview
	chmod 755 $(BIN_DIR)/mapview

mummerplot: mummerplot.pl Foundation.pm
	$(SED)  -e 's?__PERL_PATH?$(PERL)?g' \
                -e 's?__SCRIPT_DIR?$(SCRIPT_DIR)?g' \
                -e 's?__BIN_DIR?$(BIN_DIR)?g' \
                mummerplot.pl > $(BIN_DIR)/mummerplot
	chmod 755 $(BIN_DIR)/mummerplot

dnadiff: dnadiff.pl Foundation.pm
	$(SED) -e 's?__PERL_PATH?$(PERL)?g' \
               -e 's?__SCRIPT_DIR?$(SCRIPT_DIR)?g' \
               -e 's?__BIN_DIR?$(BIN_DIR)?g' \
                 dnadiff.pl > $(BIN_DIR)/dnadiff
	chmod 755 $(BIN_DIR)/dnadiff

nucmer: nucmer.pl Foundation.pm
	$(SED)  -e 's?__PERL_PATH?$(PERL)?g' \
                -e 's?__SCRIPT_DIR?$(SCRIPT_DIR)?g' \
                -e 's?__AUX_BIN_DIR?$(AUX_BIN_DIR)?g' \
                -e 's?__BIN_DIR?$(BIN_DIR)?g' \
                 nucmer.pl > $(BIN_DIR)/nucmer
	chmod 755 $(BIN_DIR)/nucmer

promer: promer.pl Foundation.pm
	$(SED)  -e 's?__PERL_PATH?$(PERL)?g' \
                -e 's?__SCRIPT_DIR?$(SCRIPT_DIR)?g' \
                -e 's?__AUX_BIN_DIR?$(AUX_BIN_DIR)?g' \
                -e 's?__BIN_DIR?$(BIN_DIR)?g' \
                promer.pl > $(BIN_DIR)/promer
	chmod 755 $(BIN_DIR)/promer

run-mummer1: run-mummer1.csh
	$(SED)  -e 's?__CSH_PATH?$(CSH)?g' \
		-e 's?__BIN_DIR?$(BIN_DIR)?g' \
		run-mummer1.csh > $(BIN_DIR)/run-mummer1
	chmod 755 $(BIN_DIR)/run-mummer1

run-mummer3: run-mummer3.csh
	$(SED)  -e 's?__CSH_PATH?$(CSH)?g' \
		-e 's?__BIN_DIR?$(BIN_DIR)?g' \
		run-mummer3.csh > $(BIN_DIR)/run-mummer3
	chmod 755 $(BIN_DIR)/run-mummer3

nucmer2xfig: nucmer2xfig.pl
	$(SED)  -e 's?__PERL_PATH?$(PERL)?g' \
		nucmer2xfig.pl > $(BIN_DIR)/nucmer2xfig
	chmod 755 $(BIN_DIR)/nucmer2xfig
