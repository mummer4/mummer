# ===========================================================================
#        http://www.gnu.org/software/autoconf-archive/ax_ruby_ext.html
# ===========================================================================
#
# SYNOPSIS
#
#   AX_RUBY_EXT
#
# DESCRIPTION
#
#   Fetches the linker flags and C compiler flags for compiling and linking
#   Ruby binary extensions.  The macro substitutes RUBY_VERSION,
#   RUBY_EXT_INC, RUBY_EXT_LIB, RUBY_EXT_CPPFLAGS, RUBY_EXT_LDFLAGS and
#   RUBY_EXT_DLEXT variables if Ruby executable has been found.  It also
#   checks the same variables before trying to retrieve them from the Ruby
#   configuration.
#
#     RUBY_VERSION: version of the Ruby interpreter
#     RUBY_EXT_INC: Ruby include directory
#     RUBY_EXT_LIB: Ruby extensions destination directory
#     RUBY_EXT_CPPFLAGS: C preprocessor flags to compile extensions
#     RUBY_EXT_LDFLAGS: linker flags to build extensions
#     RUBY_EXT_DLEXT: extensions suffix for ruby modules (e.g. "so")
#
#   Examples:
#
#     AX_RUBY_EXT
#     if test x"$RUBY" = x; then
#        AC_ERROR(["cannot find Ruby"])
#     fi
#
# LICENSE
#
#   Copyright (c) 2011 Stanislav Sedov <stas@FreeBSD.org>
#
#   Redistribution and use in source and binary forms, with or without
#   modification, are permitted provided that the following conditions are
#   met:
#
#   1. Redistributions of source code must retain the above copyright
#   notice, this list of conditions and the following disclaimer.
#
#   2. Redistributions in binary form must reproduce the above copyright
#   notice, this list of conditions and the following disclaimer in the
#   documentation and/or other materials provided with the distribution.
#
#   THIS SOFTWARE IS PROVIDED BY THE AUTHOR AND CONTRIBUTORS ``AS IS'' AND
#   ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
#   IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
#   PURPOSE ARE DISCLAIMED.  IN NO EVENT SHALL THE AUTHOR OR CONTRIBUTORS BE
#   LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
#   CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
#   SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
#   INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
#   CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
#   ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF
#   THE POSSIBILITY OF SUCH DAMAGE.

#serial 2

AC_DEFUN([AX_RUBY_EXT],[
        #
        # Check if ruby executable exists.
        #
        AC_ARG_VAR([RUBY], [the Ruby interpreter])
        AC_PATH_PROGS(RUBY, ["${RUBY-ruby}"], [])

        if test -n "$RUBY" ; then
                #
                # Check Ruby version.
                #
                AC_MSG_CHECKING([for Ruby version])
                [RUBY_VERSION=`$RUBY -e 'puts RUBY_VERSION'`];
                AC_MSG_RESULT([$RUBY_VERSION])
                AC_SUBST(RUBY_VERSION)

                #
                # Check for the extensions target directory.
                #
                AC_MSG_CHECKING([for Ruby extensions target directory])
                AS_IF([test -z "$RUBY_EXT_LIB"],
                      [RUBY_EXT_LIB=`$RUBY -rrbconfig -e 'print RbConfig::expand(RbConfig::CONFIG.fetch("sitearchdir"))'`])
                AC_MSG_RESULT([$RUBY_EXT_LIB])
                AC_SUBST(RUBY_EXT_LIB)

                #
                # Check for include flags
                #
                AC_MSG_CHECKING([for Ruby include directory])
                AS_IF([test -z "$RUBY_EXT_CFLAGS"],
                      [RUBY_EXT_CFLAGS=`$RUBY -rrbconfig -e 'RbConfig::CONFIG.store("rubyarchhdrdir", File.join(RbConfig::CONFIG.fetch("rubyhdrdir"), RbConfig::CONFIG.fetch("arch"))) unless RbConfig::CONFIG.has_key?("rubyarchhdrdir"); %w(rubyhdrdir rubyarchhdrdir).each { |k| (v = RbConfig::CONFIG.fetch(k, nil)) and print(" -I", RbConfig::expand(v)) }'`])
                AC_MSG_RESULT([$RUBY_EXT_CFLAGS])
                AC_SUBST(RUBY_EXT_CFLAGS)

                #
                # Check for lib flags
                #
                AC_MSG_CHECKING([for Ruby libs])
                AS_IF([test -z "$RUBY_EXT_LIBS"],
                      [RUBY_EXT_LIBS=`$RUBY -rrbconfig -e '%w(LIBRUBYARG_SHARED LIBS).each { |k| (v = RbConfig::CONFIG.fetch(k, nil)) and print(" ", RbConfig::expand(v)) }'`])
                AC_MSG_RESULT([$RUBY_EXT_LIBS])
                AC_SUBST(RUBY_EXT_LIBS)


                # Fix LDFLAGS for OS X.  We don't want any -arch flags here, otherwise
                # linking might fail.  We also including the proper flags to create a bundle.
                AC_MSG_CHECKING([for Ruby extra ldflags])
                case "$host" in
                *darwin*)
                        RUBY_EXT_LDFLAGS=`echo ${RUBY_EXT_LDFLAGS} | sed -e "s,-arch [[^ ]]*,,g"`
                        RUBY_EXT_LDFLAGS="${RUBY_EXT_LDFLAGS} -bundle -undefined dynamic_lookup"
                        ;;
                esac
                AC_MSG_RESULT([$RUBY_EXT_LDFLAGS])
                AC_SUBST(RUBY_EXT_LDFLAGS)
        fi
])
