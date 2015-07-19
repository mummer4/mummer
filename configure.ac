AC_INIT([UMDMUMmer], [0.0.1], [gmarcais@umd.edu])
AC_CANONICAL_HOST
AC_CONFIG_MACRO_DIR([m4])
AM_INIT_AUTOMAKE([subdir-objects foreign parallel-tests color-tests])
AM_SILENT_RULES([yes])
AC_CONFIG_SRCDIR([.])
AC_CONFIG_HEADERS([config.h])
AC_PROG_LIBTOOL

# Change default compilation flags
AC_SUBST([ALL_CXXFLAGS], [-std=c++0x])
CXXFLAGS="-std=c++0x $CXXFLAGS"
AC_LANG(C++)
AC_PROG_CXX
AC_PROG_CC

AC_CONFIG_FILES([Makefile])
AC_OUTPUT
