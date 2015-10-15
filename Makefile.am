##########################
# Automake configuration #
##########################
ACLOCAL_AMFLAGS = -I m4

AM_CPPFLAGS = -Wall -I$(srcdir)/include $(JELLYFISH2_0_CFLAGS)
AM_LDFLAGS = -lpthread

# Pre-declare all used variables.
bin_PROGRAMS =
bin_SCRIPTS =
noinst_SCRIPTS =
noinst_HEADERS =
EXTRA_DIST =
CLEANFILES =
DISTCLEANFILES =
check_PROGRAMS =
check_SCRIPTS =
TESTS =
TEST_EXTENSIONS =
BUILT_SOURCES =

# Special cleaning rules. Like CLEANFILES and DISTCLEANFILES, but for
# directories.
CLEANDIRS =
DISTCLEANDIRS =
clean-local: clean-local-dirs
distclean-local: distclean-local-dirs
.PHONY: clean-local-dirs distclean-local-dirs
clean-local-dirs:
	-rm -rf $(CLEANDIRS)
distclean-local-dirs:
	-rm -rf $(DISTCLEANDIRS)

# yaggo rules
include yaggo.mk

##############################
# Build library libumdmummer #
##############################
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
bin_PROGRAMS += mgaps postnuc show-coords show-aligns show-snps
mgaps_SOURCES = src/tigr/mgaps_main.cc
postnuc_SOURCES = src/tigr/postnuc_main.cc
show_coords_SOURCES = src/tigr/show-coords.cc src/tigr/delta.cc
show_aligns_SOURCES = src/tigr/show-aligns.cc src/tigr/delta.cc src/tigr/translate.cc
show_snps_SOURCES = src/tigr/show-snps.cc src/tigr/delta.cc src/tigr/translate.cc

###########
# scripts #
###########
# For each script, generate 2 versions: one with local path, to use
# directly; and one (in .libs/) to use after installation. Similar to
# what LIBTOOLS does.

# List of scripts to install
perl_scripts = mummerplot mapview dnadiff
perl_libs_scripts = $(patsubst %,.libs/%,$(perl_scripts))
bin_SCRIPTS += $(perl_libs_scripts)
noinst_SCRIPTS += $(perl_scripts)
EXTRA_DIST += $(patsubst %,scripts/%.pl,$(perl_scripts))
CLEANFILES += $(perl_libs_scripts) $(perl_scripts)

# Rules to generate the scripts, in .libs and in build directory.
perl_inst_subst  = sed -e 's,[@]PERL[@],$(PERL),g' \
                       -e 's,[@]BIN_DIR[@],$(bindir),g' \
                       -e 's,[@]LIB_DIR[@],$(perl_libdir),g'
perl_local_subst = sed -e 's,[@]PERL[@],$(PERL),g' \
                       -e 's,[@]BIN_DIR[@],$(abs_builddir),g' \
                       -e 's,[@]LIB_DIR[@],$(abs_srcdir)/scripts,g'
.libs/%: scripts/%.pl
	@mkdir -p $(dir $@) && $(perl_inst_subst) < $< > $@
%: scripts/%.pl
	$(AM_V_GEN)$(perl_local_subst) < $< > $@ && chmod a+x $@

# Supporting Perl module
perl_lib_SCRIPTS = scripts/Foundation.pm
perl_libdir = $(libdir)/mummer
EXTRA_DIST += scripts/Foundation.pm

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
