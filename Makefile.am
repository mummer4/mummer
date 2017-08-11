##########################
# Automake configuration #
##########################
ACLOCAL_AMFLAGS = -I m4

AM_CPPFLAGS = -Wall -I$(srcdir)/include
AM_CXXFLAGS = -std=c++11 $(EXTRA_CXXFLAGS)
AM_LDFLAGS = -lpthread

# Pre-declare all used variables.
bin_PROGRAMS =
bin_SCRIPTS =
pkglibexec_PROGRAMS =
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

# Where to install
script_libdir = $(pkglibdir)
pkgconfigdir = $(libdir)/pkgconfig
pkgconfig_DATA = mummer.pc

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

# For distcheck
AM_DISTCHECK_CONFIGURE_FLAGS=--enable-perl-binding --enable-ruby-binding --enable-python-binding

##############################
# Build library libumdmummer #
##############################
lib_LTLIBRARIES = libumdmummer.la
LDADD = libumdmummer.la
libumdmummer_la_SOURCES  = src/essaMEM/sparseSA.cpp src/essaMEM/sssort_compact.cc
libumdmummer_la_SOURCES += src/tigr/mgaps.cc src/tigr/postnuc.cc src/tigr/sw_align.cc src/tigr/tigrinc.cc
libumdmummer_la_SOURCES += src/umd/nucmer.cc

library_includedir = $(includedir)/mummer-@PACKAGE_VERSION@

# TODO: remove duplicated header and functionality between jellyfish and thread_pipe
nobase_library_include_HEADERS = include/mummer/sparseSA.hpp			\
                                 include/mummer/fasta.hpp			\
                                 include/mummer/postnuc.hh			\
                                 include/mummer/translate.hh			\
                                 include/mummer/sw_align.hh			\
                                 include/mummer/tigrinc.hh			\
                                 include/mummer/nucmer.hpp			\
                                 include/mummer/mgaps.hh			\
                                 include/mummer/delta.hh			\
                                 include/mummer/sw_alignscore.hh		\
                                 include/mummer/sparseSA_imp.hpp		\
                                 include/jellyfish/circular_buffer.hpp		\
                                 include/jellyfish/cooperative_pool2.hpp	\
                                 include/jellyfish/cpp_array.hpp		\
                                 include/jellyfish/divisor.hpp			\
                                 include/jellyfish/stream_manager.hpp		\
                                 include/jellyfish/whole_sequence_parser.hpp	\
                                 include/mummer/timer.hpp
nobase_library_include_HEADERS += include/thread_pipe.hpp			\
                                  include/thread_pipe/circular_buffer.hpp	\
                                  include/thread_pipe/cooperative_pool2.hpp	\
                                  include/thread_pipe/divisor.hpp		\
                                  include/thread_pipe/traits.hpp
nobase_library_include_HEADERS +=							\
                                  include/compactsufsort/compactsufsort.hpp		\
                                  include/compactsufsort/const_iterator_traits.hpp	\
                                  include/compactsufsort/said_traits.hpp		\
                                  include/compactsufsort/trsort_imp.hpp			\
                                  include/compactsufsort/compactsufsort_imp.hpp		\
                                  include/compactsufsort/divsufsort_private.h		\
                                  include/compactsufsort/sssort_imp.hpp
nobase_library_include_HEADERS += include/mummer/48bit_index.hpp		\
                                  include/mummer/48bit_iterator.hpp		\
                                  include/mummer/const_iterator_traits.hpp
nobase_library_include_HEADERS += include/mummer/dset.hpp		\
                                  include/mummer/openmp_qsort.hpp	\
                                  include/mt_skip_list/common.hpp	\
                                  include/mt_skip_list/set.hpp		\
                                  include/mummer/redirect_to_pager.hpp

###########
# essaMEM #
###########
bin_PROGRAMS += mummer
mummer_SOURCES = src/essaMEM/fasta.cpp src/essaMEM/mummer.cpp

########
# tigr #
########
bin_PROGRAMS += show-coords show-aligns show-snps show-tiling	\
                show-diff repeat-match annotate combineMUMs	\
                delta-filter
pkglibexec_PROGRAMS += mgaps postnuc prepro postpro
mgaps_SOURCES = src/tigr/mgaps_main.cc
postnuc_SOURCES = src/tigr/postnuc_main.cc
show_coords_SOURCES = src/tigr/show-coords.cc src/tigr/delta.cc
show_aligns_SOURCES = src/tigr/show-aligns.cc src/tigr/delta.cc src/tigr/translate.cc
show_snps_SOURCES = src/tigr/show-snps.cc src/tigr/delta.cc src/tigr/translate.cc
show_tiling_SOURCES = src/tigr/show-tiling.cc src/tigr/delta.cc
show_diff_SOURCES = src/tigr/show-diff.cc src/tigr/delta.cc
repeat_match_SOURCES = src/tigr/repeat-match.cc
annotate_SOURCES = src/tigr/annotate.cc
combineMUMs_SOURCES = src/tigr/combineMUMs.cc
delta_filter_SOURCES = src/tigr/delta-filter.cc src/tigr/delta.cc
prepro_SOURCES = src/tigr/prepro.cc src/tigr/translate.cc
postpro_SOURCES = src/tigr/postpro.cc src/tigr/translate.cc

# Obsolete stuff
EXTRA_DIST += src/tigr/gaps.cc

###########
# scripts #
###########
# For each script, generate 2 versions: one with local path, to use
# directly; and one (in .libs/) to use after installation. Similar to
# what LIBTOOLS does.

# List of scripts to install
perl_scripts = mummerplot dnadiff promer
shell_scripts = exact-tandems
all_scripts = $(perl_scripts) $(shell_scripts)

# Obsolete stuff
EXTRA_DIST += scripts/run-mummer1.sh scripts/run-mummer3.sh	\
              scripts/mapview.pl scripts/nucmer2xfig.pl

libs_scripts = $(patsubst %,.libs/%,$(all_scripts))
bin_SCRIPTS += $(libs_scripts)
noinst_SCRIPTS += $(all_scripts)
EXTRA_DIST += $(patsubst %,scripts/%.pl,$(perl_scripts))
EXTRA_DIST += $(patsubst %,scripts/%.sh,$(shell_scripts))
CLEANFILES += $(libs_scripts) $(all_scripts)

# Rules to generate the scripts, in .libs and in build directory.
script_inst_subst  = sed -e 's,[@]PERL[@],$(PERL),g' \
                       -e 's,[@]BIN_DIR[@],$(bindir),g' \
                       -e 's,[@]LIB_DIR[@],$(script_libdir),g' \
	               -e 's,[@]LIBEXEC_DIR[@],$(pkglibexecdir),g' \
		       -e 's,[@]GNUPLOT_EXE[@],$(GNUPLOT),g'
script_local_subst = sed -e 's,[@]PERL[@],$(PERL),g' \
                       -e 's,[@]BIN_DIR[@],$(abs_builddir),g' \
                       -e 's,[@]LIB_DIR[@],$(abs_srcdir)/scripts,g' \
	               -e 's,[@]LIBEXEC_DIR[@],$(abs_builddir),g' \
                       -e 's,[@]GNUPLOT_EXE[@],$(GNUPLOT),g'

.libs/%: scripts/%.pl
	@mkdir -p $(dir $@) && $(script_inst_subst) < $< > $@
%: scripts/%.pl
	$(AM_V_GEN)$(script_local_subst) < $< > $@ && chmod a+x $@
.libs/%: scripts/%.sh
	@mkdir -p $(dir $@) && $(script_inst_subst) < $< > $@
%: scripts/%.sh
	$(AM_V_GEN)$(script_local_subst) < $< > $@ && chmod a+x $@

# Supporting modules
script_lib_SCRIPTS = scripts/Foundation.pm scripts/tandem-repeat.awk
EXTRA_DIST += $(script_lib_SCRIPTS)

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
