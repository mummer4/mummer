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
	    %D%/gtest/gtest.h

EXTRA_DIST += $(GTEST_SRC)

check_PROGRAMS += %D%/test_nucmer
unittests_test_nucmer_SOURCES = %D%/test_nucmer.cc
unittests_test_nucmer_LDADD = $(LDADD) unittests/libgtest_main.la
unittests_test_nucmer_CXXFLAGS = -I$(srcdir)/unittests
TESTS += %D%/test_nucmer
