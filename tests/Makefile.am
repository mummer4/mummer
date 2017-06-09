######################################################
# Build helper programs and scripts used for testing #
######################################################
check_PROGRAMS += %D%/generate_sequences %D%/ufasta %D%/check_cigar %D%/check_LCP
check_SCRIPTS += %D%/testsh
CLEANFILES += $(check_SCRIPTS)

# Build testsh. A wrapper around /bin/shell that runs the tests in the
# proper environment
%D%/testsh: %D%/testsh.in
	$(AM_V_GEN)sed -e 's,[@]BUILDDIR[@],$(abs_builddir),g' \
	               -e 's,[@]SRCDIR[@],$(abs_srcdir),g' $< > $@
	@chmod a+x $@
EXTRA_DIST += %D%/testsh.in

# test_md5 is a wrapper around md5sum to check outputs
EXTRA_DIST += %D%/test_md5

# Build generate_sequences, used to generage pseudo-random sequences.
%C%_generate_sequences_SOURCES = %D%/generate_sequences.cc
%C%_generate_sequences_CPPFLAGS = $(CPPFLAGS) -I%D%
YAGGO_BUILT += %D%/generate_sequences_cmdline.hpp

# Build ufasta, used to sort fasta files for reliable comparisons.
%C%_ufasta_SOURCES = %D%/ufasta_src/sort.cc %D%/ufasta_src/dsort.cc	\
 %D%/ufasta_src/ufasta.cc %D%/ufasta_src/common.hpp
%C%_ufasta_CPPFLAGS = $(CPPFLAGS) -I%D%/ufasta_src -I$(srcdir)/tests/ufasta_src
YAGGO_BUILT += %D%/ufasta_src/dsort_cmdline.hpp %D%/ufasta_src/sort_cmdline.hpp

# Build check_cigar, to check validity of cigar strings.
%C%_check_cigar_SOURCES = %D%/check_cigar.cc
%C%_check_cigar_CPPFLAGS = $(CPPFLAGS) -I%D%
YAGGO_BUILT += %D%/check_cigar_cmdline.hpp

# Build check_LCP. Builds LCP from .sa file.
%C%_check_LCP_SOURCES = %D%/check_LCP.cc src/essaMEM/fasta.cpp
%C%_check_LCP_CPPFLAGS = $(AM_CPPFLAGS) -I%D%
YAGGO_BUILT += %D%/check_LCP_cmdline.hpp

####################################
# Generate pseudo random sequences #
####################################
# small_reads_*.fa are random sequences, which should change every time
# they are generated. Used in consistency tests.
%D%/data/small_reads%0.fa %D%/data/small_reads%1.fa: %D%/generate_sequences
	$(AM_V_GEN)mkdir -p $(dir $@); %D%/generate_sequences -p %D%/data/small -G 1M --numbers 1000 --numbers 1000

# seed_reads_*.fa are pseudo-random sequences. They are identical
# every time. Used for testing identical results with older version of
# MUMmer.
EXTRA_DIST += %D%/data_seed
%D%/data/seed_reads%0.fa %D%/data/seed_reads%1.fa %D%/data/seed_reads%2.fa: %D%/generate%sequences
	$(AM_V_GEN)mkdir -p $(dir $@); %D%/generate_sequences --seed-file $(srcdir)/tests/data_seed \
	  -p %D%/data/seed -G 1M \
	  --errors 1 --lengths 100 --numbers 1000 \
	  --errors 5 --lengths 1000 --numbers 1000 \
	  --errors 2 --lengths 50k --numbers 10

# Clean data directory
CLEANDIRS += %D%/data

#########################
# Tests scripts to run. #
#########################
# Shell script (.sh) tests are run by are wrapper testsh.
TEST_EXTENSIONS += .sh
SH_LOG_COMPILER = %D%/testsh

# List of tests to run
script_tests = %D%/save_load.sh %D%/batch.sh %D%/mummer.sh %D%/nucmer.sh %D%/sam.sh %D%/genome.sh %D%/delta-filter.sh
EXTRA_DIST += $(script_tests)
TESTS += $(script_tests)

# testsh creates a special directory for each test when it
# runs. Clean them.
CLEANDIRS += $(patsubst %D%/%.sh,%D%/%,$(script_tests))

# Dependencies of the tests on the generated data.
%D%/save_load.log: %D%/data/small_reads_0.fa %D%/data/small_reads_1.fa
%D%/batch.log: %D%/data/small_reads_0.fa %D%/data/small_reads_1.fa

%D%/mummer.log: %D%/data/seed_reads_0.fa %D%/data/seed_reads_1.fa
%D%/nucmer.log: %D%/data/seed_reads_0.fa %D%/data/seed_reads_1.fa
%D%/sam.log: %D%/data/seed_reads_0.fa %D%/data/seed_reads_1.fa
%D%/genome.log: %D%/data/seed_reads_2.fa
%D%/delta-filter.log: %D%/data/small_reads_0.fa %D%/data/small_reads_1.fa
