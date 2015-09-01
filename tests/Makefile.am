AM_CPPFLAGS = -Wall -Werror -g -I$(top_srcdir)/include -I. $(JELLYFISH2_0_CFLAGS)
AM_LDFLAGS = -lpthread
LDADD = ../libumdmummer.la $(JELLYFISH2_0_LIBS) libgtest_main.la

##########################
# Google testing library #
##########################
# Build rules for libraries.
check_LTLIBRARIES = libgtest.la libgtest_main.la

libgtest_la_SOURCES = gtest/src/gtest-all.cc
libgtest_main_la_SOURCES = gtest/src/gtest_main.cc
libgtest_main_la_LIBADD = libgtest.la
libgtest_la_CXXFLAGS = -I$(srcdir)/tests
libgtest_main_la_CXXFLAGS = -I$(srcdir)/tests

GTEST_SRC = gtest/src/gtest-all.cc	\
	    gtest/src/gtest_main.cc	\
	    gtest/gtest.h

EXTRA_DIST = $(GTEST_SRC)


#########
# Tests #
#########
check_PROGRAMS = test_nucmer
test_nucmer_SOURCES = test_nucmer.cc

TEST_EXTENSIONS =
TESTS = test_nucmer

# EXTRA_DIST = $(TESTS)
clean-local: clean-local-check
.PHONY: clean-local-check
clean-local-check:
	-cd tests; rm -f *
