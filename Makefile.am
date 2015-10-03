include yaggo.mk
# SUBDIRS = . tests
ACLOCAL_AMFLAGS = -I m4

AM_CPPFLAGS = -Wall -I$(srcdir)/include $(JELLYFISH2_0_CFLAGS)
AM_LDFLAGS = -lpthread

######################
# Automake variables #
######################
bin_PROGRAMS =
EXTRA_DIST =
CLEANFILES =
DISTCLEANFILES =
check_PROGRAMS =
TESTS =
TEST_EXTENSIONS =

###########
# library #
###########
lib_LTLIBRARIES = libumdmummer.la
LDADD = libumdmummer.la
libumdmummer_la_SOURCES  = src/essaMEM/sparseSA.cpp src/essaMEM/qsufsort.c
libumdmummer_la_SOURCES += src/tigr/mgaps.cc src/tigr/postnuc.cc src/tigr/sw_align.cc src/tigr/tigrinc.cc
libumdmummer_la_SOURCES += src/umd/nucmer.cc
libumdmummer_la_LIBS = $(JELLYFISH2_0_LIBS)
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
bin_PROGRAMS += nucmer
YAGGO_BUILT += src/umd/nucmer_cmdline.hpp
nucmer_SOURCES = src/umd/nucmer_main.cc

#################
# SWIG bindings #
#################
include swig/Makefile.am

#########
# Tests #
#########
include unittests/Makefile.am
include tests/Makefile.am
