BUILD = build
PREFIX = $(shell pwd)/inst

INSTALL = install
BSCRIPTS = $(BUILD)/scripts
BTIGR = $(BUILD)/src/tigr
BMUMMER = $(BUILD)/src/kurtz/mm3src

all: tup

.PHONY: tup

tup:
	tup upd

install: tup
	@echo Installing into $(PREFIX)
	$(INSTALL) -d $(PREFIX)/aux_bin
	$(INSTALL) $(BSCRIPTS)/exact-tandems $(BSCRIPTS)/mapview \
	           $(BSCRIPTS)/mummerplot $(BSCRIPTS)/nucmer2xfig \
	           $(BSCRIPTS)/dnadiff $(BSCRIPTS)/nucmer $(BSCRIPTS)/promer \
	           $(BSCRIPTS)/run-mummer1 $(BSCRIPTS)/run-mummer3 \
		   $(BTIGR)/annotate $(BTIGR)/combineMUMs $(BTIGR)/delta-filter \
	           $(BTIGR)/gaps $(BTIGR)/repeat-match $(BTIGR)/show-aligns \
	           $(BTIGR)/show-coords $(BTIGR)/show-tiling \
	           $(BTIGR)/show-snps $(BTIGR)/show-diff \
	           $(BMUMMER)/mummer \
	           $(PREFIX)
	$(INSTALL) $(BTIGR)/postnuc $(BTIGR)/postpro $(BTIGR)/prenuc \
	           $(BTIGR)/prepro \
	           $(PREFIX)/aux_bin

