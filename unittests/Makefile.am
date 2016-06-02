##########################
# Google testing library #
##########################
# Build rules for libraries.
check_LTLIBRARIES = %D%/libgtest.la %D%/libgtest_main.la

%C%_libgtest_la_SOURCES = %D%/gtest/src/gtest-all.cc
%C%_libgtest_main_la_SOURCES = %D%/gtest/src/gtest_main.cc
%C%_libgtest_main_la_LIBADD = %D%/libgtest.la
%C%_libgtest_la_CXXFLAGS = -I$(srcdir)/unittests
%C%_libgtest_main_la_CXXFLAGS = -I$(srcdir)/unittests

GTEST_SRC = %D%/gtest/src/gtest-all.cc	\
	    %D%/gtest/src/gtest_main.cc	\
	    %D%/gtest/gtest.h %D%/gtest/test.hpp

EXTRA_DIST += $(GTEST_SRC)

#####################
# Unittest programs #
#####################
unittests_programs = %D%/test_all
check_PROGRAMS += $(unittests_programs)
TESTS += $(unittests_programs)

%C%_test_all_SOURCES = %D%/test_nucmer.cc				\
 %D%/test_cooperative_pool2.cc %D%/test_whole_sequence_parser.cc	\
 %D%/test_sparse_sa.cc %D%/test_qsort.cc
%C%_test_all_LDADD = $(LDADD) %D%/libgtest_main.la
%C%_test_all_CXXFLAGS = $(AM_CXXFLAGS) -I$(srcdir)/unittests


