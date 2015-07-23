# SWIG
BUILT_SOURCES =
CLEANFILES =
EXTRA_DIST =
SWIG_SRC = mummer.i nucmer.i

LIBNAME = mummer

if HAVE_SWIG
SWIG_V_GEN = $(swig_v_GEN_$(V))
swig_v_GEN_ = $(swig_v_GEN_$(AM_DEFAULT_VERBOSITY))
swig_v_GEN_0 = @echo "  SWIG    " $@;
%/swig_wrap.cpp: $(SWIG_SRC)
	$(SWIG_V_GEN)$(SWIG) -$* -I$(srcdir)/../src -o $@ $<
else
%/swig_wrap.cc:
	@echo >&2 SWIG >= 3.x.x not found. Make sure it is install and rerun configure
	@false
endif

# Python support
if PYTHON_BINDING
PYTHON_BUILT = python/swig_wrap.cpp python/$(LIBNAME).py
BUILT_SOURCES += $(PYTHON_BUILT)

pythonextdir = $(PYTHON_SITE_PKG)/$(LIBNAME)
pythonext_SCRIPTS = python/__init__.pyc
pythonext_LTLIBRARIES = python/_mummer.la
python__mummer_la_SOURCES = python/swig_wrap.cpp $(SWIG_SRC)
python__mummer_la_CPPFLAGS = $(PYTHON_CPPFLAGS) -I$(srcdir)/../src
python__mummer_la_LDFLAGS = -module $(PYTHON_LDFLAGS)
python__mummer_la_LIBADD = ../libumdmummer.la
CLEANFILES += $(PYTHON_BUILT) $(dist_pythonext_DATA) python/$(LIBNAME).py
PYTHONC_V_GEN = $(pythonc_v_GEN_$(V))
pythonc_v_GEN_ = $(pythonc_v_GEN_$(AM_DEFAULT_VERBOSITY))
pythonc_v_GEN_0 = @echo "  PYTHONC " $@;
%/__init__.pyc: %/$(LIBNAME).py
	$(PYTHONC_V_GEN)$(PYTHON) -c 'import py_compile, sys; py_compile.compile(sys.argv[1], sys.argv[2])' $< $@
python/$(LIBNAME).py: python/swig_wrap.cpp
EXTRA_DIST += python/$(LIBNAME).py
endif

# Ruby support
if RUBY_BINDING
RUBY_BUILT = ruby/swig_wrap.cpp
BUILT_SOURCES += $(RUBY_BUILT)
rubyextdir = $(RUBY_EXT_LIB)
rubyext_LTLIBRARIES = ruby/mummer.la
ruby_mummer_la_SOURCES = ruby/swig_wrap.cpp $(SWIG_SRC)
ruby_mummer_la_CPPFLAGS = $(RUBY_EXT_CFLAGS) -I$(srcdir)/../src
ruby_mummer_la_LDFLAGS = -module $(RUBY_EXT_LIBS)
ruby_mummer_la_LIBADD = ../libumdmummer.la
CLEANFILES += $(RUBY_BUILT)
endif

# Perl5 support
if PERL_BINDING
PERL_BUILT = perl5/swig_wrap.cpp perl5/$(LIBNAME).pm
BUILT_SOURCES += $(PERL_BUILT)
perlextdir = $(PERL_EXT_LIB)
perlext_SCRIPTS = perl5/$(LIBNAME).pm
perlext_LTLIBRARIES = perl5/mummer.la
perl5_mummer_la_SOURCES = perl5/swig_wrap.cpp $(SWIG_SRC)
perl5_mummer_la_CPPFLAGS = $(PERL_EXT_CPPFLAGS) -I$(PERL_EXT_INC) -I$(srcdir)/../src
perl5_mummer_la_LDFLAGS = -module $(PERL5_EXT_LDFLAGS)
perl5_mummer_la_LIBADD = ../libumdmummer.la
CLEANFILES += $(PERL5_BUILT) perl5/$(LIBNAME).pm
perl5/$(LIBNAME).pm: perl5/swig_wrap.cpp
EXTRA_DIST += perl5/$(LIBNAME).pm
endif
