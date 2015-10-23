# SWIG
EXTRA_DIST += $(SWIG_TESTS)
SWIG_SRC = %D%/mummer.i %D%/nucmer.i
TESTS += $(SWIG_TESTS)
SWIG_TESTS =

# Definitions of test runners
TEST_EXTENSIONS += .py .rb .pl
PY_LOG_COMPILER = $(ENV_PROG)
AM_PY_LOG_FLAGS = PYTHONPATH=$(builddir)/%D%/python:$(builddir)/%D%/python/.libs $(PYTHON)
RB_LOG_COMPILER = $(RUBY)
AM_RB_LOG_FLAGS = -I$(builddir)/%D%/ruby/.libs
PL_LOG_COMPILER = $(PERL)
AM_PL_LOG_FLAGS = -I$(builddir)/%D%/perl5 -I$(builddir)/%D%/perl5/.libs -I$(srcdir)/%D%/perl5

SWIG_CPPFLAGS = -I$(srcdir)/include
SWIG_LIBS =

if HAVE_SWIG
SWIG_V_GEN = $(swig_v_GEN_$(V))
swig_v_GEN_ = $(swig_v_GEN_$(AM_DEFAULT_VERBOSITY))
swig_v_GEN_0 = @echo "  SWIG    " $@;
%/swig_wrap.cpp: $(SWIG_SRC)
	$(SWIG_V_GEN)$(SWIG) -$(notdir $*) $(SWIG_CPPFLAGS) -o $@ $<
else
%/swig_wrap.cc:
	@echo >&2 SWIG >= 3.x.x not found. Make sure it is install and rerun configure
	@false
endif

# Python support
if PYTHON_BINDING
PYTHON_BUILT = %D%/python/swig_wrap.cpp %D%/python/mummer.py
BUILT_SOURCES += $(PYTHON_BUILT)
EXTRA_DIST += $(PYTHON_BUILT)

pythonextdir = $(PYTHON_SITE_PKG)/mummer
pythonext_SCRIPTS = %D%/python/__init__.pyc
pythonext_LTLIBRARIES = %D%/python/_mummer.la
%C%_python__mummer_la_SOURCES = %D%/python/swig_wrap.cpp $(SWIG_SRC)
%C%_python__mummer_la_CPPFLAGS = $(PYTHON_CPPFLAGS) $(SWIG_CPPFLAGS)
%C%_python__mummer_la_LDFLAGS = -module
%C%_python__mummer_la_LIBS = $(SWIG_LIBS)
%C%_python__mummer_la_LIBADD = libumdmummer.la
DISTCLEANFILES += $(PYTHON_BUILT) $(pythonext_SCRIPTS)
PYTHONC_V_GEN = $(pythonc_v_GEN_$(V))
pythonc_v_GEN_ = $(pythonc_v_GEN_$(AM_DEFAULT_VERBOSITY))
pythonc_v_GEN_0 = @echo "  PYTHONC " $@;
%/__init__.pyc: %/mummer.py
	$(PYTHONC_V_GEN)$(PYTHON) -c 'import py_compile, sys; py_compile.compile(sys.argv[1], sys.argv[2])' $< $@
%D%/python/mummer.py: %D%/python/swig_wrap.cpp

# Tests
SWIG_TESTS += %D%/python/test_nucmer.py
endif

# Ruby support
if RUBY_BINDING
RUBY_BUILT = %D%/ruby/swig_wrap.cpp
BUILT_SOURCES += $(RUBY_BUILT)
rubyextdir = $(RUBY_EXT_LIB)
rubyext_LTLIBRARIES = %D%/ruby/mummer.la
%C%_ruby_mummer_la_SOURCES = %D%/ruby/swig_wrap.cpp $(SWIG_SRC)
%C%_ruby_mummer_la_CPPFLAGS = $(RUBY_EXT_CFLAGS) $(SWIG_CPPFLAGS)
%C%_ruby_mummer_la_LDFLAGS = -module
%C%_ruby_mummer_la_LIBS = $(SWIG_LIBS)
%C%_ruby_mummer_la_LIBADD = libumdmummer.la
DISTCLEANFILES += $(RUBY_BUILT)

# Tests
SWIG_TESTS += %D%/ruby/test_nucmer.rb
endif

# Perl5 support
if PERL_BINDING
PERL_BUILT = %D%/perl5/swig_wrap.cpp %D%/perl5/mummer.pm
BUILT_SOURCES += $(PERL_BUILT)
EXTRA_DIST += $(PERL_BUILT)

perlextdir = $(PERL_EXT_LIB)
perlext_SCRIPTS = %D%/perl5/mummer.pm
perlext_LTLIBRARIES = %D%/perl5/mummer.la
%C%_perl5_mummer_la_SOURCES = %D%/perl5/swig_wrap.cpp $(SWIG_SRC)
%C%_perl5_mummer_la_CPPFLAGS = $(PERL_EXT_CPPFLAGS) -I$(PERL_EXT_INC) $(SWIG_CPPFLAGS)
%C%_perl5_mummer_la_LDFLAGS = -module
%C%_perl5_mummer_la_LIBS = $(SWIG_LIBS)
%C%_perl5_mummer_la_LIBADD = libumdmummer.la
DISTCLEANFILES += $(PERL_BUILT)
%D%/perl5/mummer.pm: %D%/perl5/swig_wrap.cpp

# Tests
SWIG_TESTS += %D%/perl5/test_nucmer.pl
endif
