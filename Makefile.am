ACLOCAL_AMFLAGS = -I m4

AM_CPPFLAGS = -Wall -I$(srcdir)/src
AM_LDFLAGS = -lpthread

bin_PROGRAMS =

###########
# library #
###########
lib_LTLIBRARIES = libumdmummer.la
LDADD = libumdmummer.la
libumdmummer_la_SOURCES  = src/essaMEM/sparseSA.cpp src/essaMEM/qsufsort.c
libumdmummer_la_SOURCES += src/tigr/mgaps.cc src/tigr/postnuc.cc src/tigr/sw_align.cc src/tigr/tigrinc.cc
libumdmummer_la_SOURCES += src/umd/nucmer.cc



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
