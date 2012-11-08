#-- Imported variables from top level makefile
# BIN_DIR AUX_BIN_DIR CXX CC CFLAGS CXXFLAGS LDFLAGS

##include ../Makedef

SPLINTFLAGS=-f ../Splintoptions -DDEBUG

LD=$(CC)

##CFLAGS=${DEFINECFLAGS}

LIBBASE=libbase.a
LIBBASEDBG=libbase.dbg.a

.SUFFIXES: .dbg.o .c .o .inc

include Filelists.mf

all: $(LIBBASE) $(LIBBASEDBG)


protodef.h:
	Mkprotodef.sh `ls *.c` > protodef.h


$(LIBBASE): $(LIBOBJECTS)
	ar sruv $@ $(LIBOBJECTS)


include Filegoals.mf


$(LIBBASEDBG): $(LIBDEBUGOBJECTS)
	ar sruv $@ $(LIBDEBUGOBJECTS)


.PHONY:clean
clean:splintclean
	rm -f *.[oa] *.dbg.o *.inc *.aux *.dvi *.log *~


include Dependencies.mf
