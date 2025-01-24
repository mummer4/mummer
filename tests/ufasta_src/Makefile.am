check_PROGRAMS += %D%/ufasta
%C%_SUB_CMDS = %D%/one.cc %D%/sizes.cc %D%/head.cc %D%/tail.cc %D%/rc.cc		\
               %D%/n50.cc %D%/extract.cc %D%/format.cc %D%/sort.cc %D%/dsort.cc \
               %D%/split.cc	%D%/rsort.cc
%C%_ufasta_SOURCES = %D%/ufasta.cc %D%/common.cc $(SUB_CMDS)
%C%_ufasta_CPPFLAGS = $(AM_CPPFLAGS) -I%D% -I$(srcdir)/%D%
YAGGO_BUILT += $(SUB_CMDS:.cc=_cmdline.hpp)
noinst_HEADERS += %D%/common.hpp %D%/seeded_prg.hpp
