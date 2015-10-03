check_PROGRAMS += %D%/generate_sequences

%C%_generate_sequences_SOURCES = %D%/generate_sequences.cc
%C%_generate_sequences_CPPFLAGS = $(CPPFLAGS) -I%D%
YAGGO_BUILT += %D%/generate_sequences_cmdline.hpp

%D%/data/small_reads%0.fa %D%/data/small_reads%1.fa: %D%/generate_sequences
	$(AM_V_GEN)mkdir -p $(dir $@); %D%/generate_sequences -p %D%/data/small -G 1M --numbers 1000 --numbers 1000
CLEANFILES += %D%/small_reads_0.fa %D%/small_reads_1.fa

TEST_EXTENSIONS += .sh
SH_LOG_COMPILER = $(srcdir)/tests/testsh
TESTS += %D%/save_load.sh

%D%/save_load.log: %D%/data/small_reads_0.fa %D%/data/small_reads_1.fa

