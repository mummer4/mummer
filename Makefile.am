SUBDIRS = . $(MAYBE_SWIG)
ACLOCAL_AMFLAGS = -I m4

AM_CPPFLAGS = -Wall -I$(srcdir)/include
AM_LDFLAGS = -lpthread

bin_PROGRAMS =
EXTRA_DIST =

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
