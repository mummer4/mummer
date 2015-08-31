SUBDIRS = . $(MAYBE_SWIG)
ACLOCAL_AMFLAGS = -I m4

AM_CPPFLAGS = -Wall -I$(srcdir)/include $(JELLYFISH2_0_CFLAGS)
AM_LDFLAGS = -lpthread $(JELLYFISH2_0_LIBS)

bin_PROGRAMS =
EXTRA_DIST =

# Yaggo automatic rules with silencing
V_YAGGO = $(V_YAGGO_$(V))
V_YAGGO_ = $(V_YAGGO_$(AM_DEFAULT_VERBOSITY))
V_YAGGO_0 = @echo "  YAGGO   " $@;
.yaggo.hpp:
	$(V_YAGGO)$(YAGGO) -o $@ $<

YAGGO_BUILT = # Append all file to be built by yaggo
BUILT_SOURCES = $(YAGGO_BUILT)

###########
# library #
###########
lib_LTLIBRARIES = libumdmummer.la
LDADD = libumdmummer.la
libumdmummer_la_SOURCES  = src/essaMEM/sparseSA.cpp src/essaMEM/qsufsort.c
libumdmummer_la_SOURCES += src/tigr/mgaps.cc src/tigr/postnuc.cc src/tigr/sw_align.cc src/tigr/tigrinc.cc
libumdmummer_la_SOURCES += src/umd/nucmer.cc
include_HEADERS = include/mummer/sparseSA.hpp				\
                  include/mummer/fasta.hpp include/mummer/postnuc.hh	\
                  include/mummer/translate.hh				\
                  include/mummer/sw_align.hh				\
                  include/mummer/tigrinc.hh include/mummer/nucmer.hpp	\
                  include/mummer/mgaps.hh include/mummer/delta.hh	\
                  include/mummer/sw_alignscore.hh



###########
# essaMEM #
###########
bin_PROGRAMS += essamummer
essamummer_SOURCES = src/essaMEM/fasta.cpp src/essaMEM/mummer.cpp

########
# tigr #
########
bin_PROGRAMS += mgaps postnuc
mgaps_SOURCES = src/tigr/mgaps_main.cc
postnuc_SOURCES = src/tigr/postnuc_main.cc

#######
# umd #
#######
bin_PROGRAMS += simplenucmer
YAGGO_BUILT += src/umd/simplenucmer_cmdline.hpp
simplenucmer_SOURCES = src/umd/simplenucmer.cc

#########
# Tests #
#########
TEST_EXTENSIONS =
TESTS =

EXTRA_DIST += $(TESTS)
clean-local: clean-local-check
.PHONY: clean-local-check
clean-local-check:
	-cd tests; rm -f *
