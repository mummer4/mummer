#-- Imported variables from top level makefile
# BIN_DIR AUX_BIN_DIR CXX CC CFLAGS CXXFLAGS LDFLAGS

ifndef BIN_DIR
BIN_DIR := $(CURDIR)
endif
ifndef AUX_BIN_DIR
AUX_BIN_DIR := $(CURDIR)
endif

OBJ_RULE = $(CXX) $(CXXFLAGS) $< -c -o $@
BIN_RULE = $(CXX) $(CXXFLAGS) $^ -o $(BIN_DIR)/$@; \
           chmod 755 $(BIN_DIR)/$@
AUX_BIN_RULE = $(CXX) $(CXXFLAGS) $^ -o $(AUX_BIN_DIR)/$@; \
           chmod 755 $(AUX_BIN_DIR)/$@
VPATH := $(AUX_BIN_DIR):$(BIN_DIR)

ALL := annotate combineMUMs delta-filter gaps mgaps \
     postnuc postpro prenuc prepro repeat-match \
     show-aligns show-coords show-tiling show-snps \
     show-diff


#-- PHONY rules --#
.PHONY: all clean


all: $(ALL)


clean:
	rm -f *.o *~
	cd $(BIN_DIR); rm -f $(ALL)
	cd $(AUX_BIN_DIR); rm -f $(ALL)




#-- not so PHONY rules --#
delta.o: delta.cc delta.hh
	$(OBJ_RULE)

tigrinc.o: tigrinc.cc tigrinc.hh
	$(OBJ_RULE)

sw_align.o: sw_align.cc sw_align.hh tigrinc.hh
	$(OBJ_RULE)

translate.o: translate.cc translate.hh
	$(OBJ_RULE)


annotate: annotate.cc tigrinc.o
	$(BIN_RULE)

combineMUMs: combineMUMs.cc tigrinc.o
	$(BIN_RULE)

delta-filter: delta-filter.cc tigrinc.o delta.o
	$(BIN_RULE)

gaps: gaps.cc tigrinc.o
	$(BIN_RULE)

mgaps: mgaps.cc tigrinc.o
	$(BIN_RULE)

postnuc: postnuc.cc tigrinc.o sw_align.o
	$(AUX_BIN_RULE)

postpro: postpro.cc tigrinc.o sw_align.o translate.o
	$(AUX_BIN_RULE)

prenuc: prenuc.cc tigrinc.o
	$(AUX_BIN_RULE)

prepro: prepro.cc tigrinc.o translate.o
	$(AUX_BIN_RULE)

repeat-match: repeat-match.cc tigrinc.o
	$(BIN_RULE)

show-aligns: show-aligns.cc tigrinc.o translate.o delta.o
	$(BIN_RULE)

show-coords: show-coords.cc tigrinc.o delta.o
	$(BIN_RULE)

show-tiling: show-tiling.cc tigrinc.o delta.o
	$(BIN_RULE)

show-snps: show-snps.cc tigrinc.o translate.o delta.o
	$(BIN_RULE)

show-diff: show-diff.cc tigrinc.o delta.o
	$(BIN_RULE)
