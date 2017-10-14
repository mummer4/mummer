AC_INIT([MUMmer], [4.0.0beta2], [gmarcais@umd.edu])
AC_CANONICAL_HOST
AC_CONFIG_MACRO_DIR([m4])
AM_INIT_AUTOMAKE([subdir-objects foreign parallel-tests color-tests])
AM_SILENT_RULES([yes])
AC_CONFIG_SRCDIR([.])
AC_CONFIG_HEADERS([config.h])
AC_PROG_LIBTOOL

# Check for compiler and flags
AS_IF([test -z "$CXXFLAGS"], [cxxflags_given=no], [cxxflags_given=yes])
AC_LANG(C++)
AC_PROG_CXX
AC_PROG_CC
AC_PATH_PROG(ENV_PROG, [env])

# Add -std=c++0x always. Add -DNDEBUG only if CXXFLAGS is not given on command line
CXXFLAGS="-std=c++0x $CXXFLAGS"
AS_IF([test x$cxxflags_given = xno], [CXXFLAGS="-DNDEBUG $CXXFLAGS"])
EXTRA_CXXFLAGS=
AC_SUBST([EXTRA_CXXFLAGS])

# Check for openmp
AC_ARG_ENABLE([openmp], [AS_HELP_STRING([--disable-openmp], [Disable OpenMP, even if available])])
AS_IF([test "x$enable_openmp" != "xno"],
      # OpenMP enabled. Check if available
      [AX_OPENMP([AC_DEFINE(HAVE_OPENMP,1,[Define if OpenMP is enabled])])]
      [EXTRA_CXXFLAGS="$EXTRA_CXXFLAGS $OPENMP_CXXFLAGS"],
      # OpenMP disabled. Ignore pragma warnings
      [EXTRA_CXXFLAGS="$EXTRA_CXXFLAGS -Wno-unknown-pragmas"])

# Check for yaggo
AC_ARG_VAR([YAGGO], [Yaggo switch parser generator])
AS_IF([test "x$YAGGO" = "x"], [AC_PATH_PROG([YAGGO], [yaggo], [false])])

# Check for gnuplot
AC_ARG_VAR([GNUPLOT], [Gnuplot plotting program])
AS_IF([test "x$GNUPLOT" = "x"], [AC_PATH_PROGS([GNUPLOT], [gnuplot gnuplot5 gnuplot4], [false])])

# Check that type __int128 is supported and if the
AC_ARG_WITH([int128],
            [AS_HELP_STRING([--with-int128], [enable int128])],
            [], [with_int128=yes])
AS_IF([test "x$with_int128" != xno],
      [AC_CHECK_TYPE([__int128],
                     [AC_DEFINE([HAVE_INT128], [1], [Define if type __int128 is supported])])])


# Check if std::thread::sleep_for is supported
AC_MSG_CHECKING([if std::thread::sleep_for is supported])
try_nano_sleep=
AC_COMPILE_IFELSE([AC_LANG_PROGRAM([[#include <thread>]],
                                   [[std::this_thread::sleep_for(std::chrono::seconds(1))]])],
                  [AC_MSG_RESULT([yes])], [try_nano_sleep=yes])
AS_IF([test x$try_nano_sleep = xyes],
      [save_cxxflags="$CXXFLAGS"]
      [CXXFLAGS="-D_GLIBCXX_USE_NANOSLEEP $CXXFLAGS"]
      [AC_COMPILE_IFELSE([AC_LANG_PROGRAM([[#include <thread>]],
                                         [[std::this_thread::sleep_for(std::chrono::seconds(1))]])],
                         [AC_MSG_RESULT([yes])]
                         [EXTRA_CXXFLAGS="-D_GLIBCXX_USE_NANOSLEEP $EXTRA_CXXFLAGS"],
                         [AC_MSG_RESULT([no])]
                         [AC_MSG_FAILURE([Compiler does not support std::this_thread::sleep_for])])]
      [CXXFLAGS="$save_cxxflags"])

#
# SWIG and bindings
#
need_shared_libraries=no
# --enable-python-binding
AC_ARG_ENABLE([python-binding],
              [AC_HELP_STRING([--enable-python-binding@<:@=PATH@:>@], [create SWIG python module and install in PATH])])
# --enable-ruby-binding
AC_ARG_ENABLE([ruby-binding],
              [AC_HELP_STRING([--enable-ruby-binding@<:@=PATH@:>@], [create SWIG ruby module and install in PATH])])
# --enable-perl-binding
AC_ARG_ENABLE([perl-binding],
              [AC_HELP_STRING([--enable-perl-binding@<:@=PATH@:>@], [create SWIG perl module and install in PATH])])
# --enable-all-binding
AC_ARG_ENABLE([all-binding],
              [AC_HELP_STRING([--enable-all-binding@<:@=PATH@:>@], [create SWIG module for all languages and install in PATH])])


# --enable-swig
AC_ARG_ENABLE([swig],
              [AC_HELP_STRING([--enable-swig], [enable development of swig binding])])
AS_IF([test x$enable_swig = xyes],
      [AX_PKG_SWIG([3.0.0], [], [AC_MSG_ERROR([SWIG version 3 is required])])])
AS_IF([test -n "$SWIG"],
      [SWIG_ENABLE_CXX])
AM_CONDITIONAL([HAVE_SWIG], [test -n "$SWIG"])

# Python binding setup
AS_IF([test -z "$enable_python_binding"], [enable_python_binding="$enable_all_binding"])
AM_CONDITIONAL(PYTHON_BINDING, [test -n "$enable_python_binding" -a x$enable_python_binding != xno])
AM_COND_IF([PYTHON_BINDING],
           [AS_IF([test "x$enable_python_binding" != xyes], [PYTHON_SITE_PKG=$enable_python_binding],
                  [test "x$prefix" != xNONE], [PYTHON_SITE_PKG=$libdir/python])]
           [AX_PYTHON_DEVEL]
           [need_shared_libraries=yes])

# Ruby binding setup
AS_IF([test -z "$enable_ruby_binding"], [enable_ruby_binding="$enable_all_binding"])
AM_CONDITIONAL([RUBY_BINDING], [test -n "$enable_ruby_binding" -a x$enable_ruby_binding != xno])
AM_COND_IF([RUBY_BINDING],
           [AS_IF([test "x$enable_ruby_binding" != xyes], [RUBY_EXT_LIB=$enable_ruby_binding],
                  [test "x$prefix" != xNONE], [RUBY_EXT_LIB=$libdir/ruby])]
           [AX_RUBY_EXT]
           [need_shared_libraries=yes])

# Perl binding setup
AS_IF([test -z "$enable_perl_binding"], [enable_perl_binding="$enable_all_binding"])
AM_CONDITIONAL([PERL_BINDING], [test -n "$enable_perl_binding" -a x$enable_perl_binding != xno])
AM_COND_IF([PERL_BINDING],
           [AS_IF([test "x$enable_perl_binding" != xyes], [PERL_EXT_LIB=$enable_perl_binding],
                  [test "x$prefix" != xNONE], [PERL_EXT_LIB=$libdir/perl])]
           [need_shared_libraries=yes])
# Look for Perl even if no binding desired. It is required for some
# scripts (e.g. mummerplot)
AX_PERL_EXT
AS_IF([test x"$PERL" = x], [AC_MSG_ERROR([Perl5 not found])])

# Check that shared libraries are enabled if needed by SWIG bindings
AS_IF([test "x$need_shared_libraries" = xyes -a "x$enable_shared" = xno],
      [AC_MSG_ERROR([Shared libraries are required for the SWIG script binding but are disabled (probably by the switch --disable-shared)])])
AS_IF([test "x$need_shared_libraries" = xyes -a "x$can_build_shared" = xno],
      [AC_MSG_ERROR([Shared libraries (modules) are required for the SWIG script binding but cannot be build (probably not supported by the compiler)])])

AC_MSG_NOTICE([compiling with flags: CXXFLAGS=$CXXFLAGS])

AC_CONFIG_FILES([Makefile mummer.pc])
AC_OUTPUT
