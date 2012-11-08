#-- Imported variables from top level makefile
# BIN_DIR AUX_BIN_DIR CXX CC CFLAGS CXXFLAGS LDFLAGS

##include ../Makedef

LIBBASEDIR=../libbasedir
INCLUDEDIR=../libbasedir

# the following flags determines the maximal allowed input size of
# the suffix tree construction
# STREESMALL means that the maximal sequence length is 2 MB
# STREELARGE means that the maximal sequence length is 134 MB
# STREEHUGE  means that the maximal sequence length is 500 MB

#SIZEFLAG=-DSTREESMALL
#SIZEFLAG=-DSTREELARGE
SIZEFLAG=-DSTREEHUGE

override CFLAGS+=-I$(INCLUDEDIR) $(SIZEFLAG)
##CFLAGS=${DEFINECFLAGS} -I$(INCLUDEDIR) $(SIZEFLAG)
##LDFLAGS=${DEFINELDFLAGS}
SPLINTFLAGS=${SIZEFLAG} -I${INCLUDEDIR} -f ../Splintoptions -DDEBUG

#-DSTARTFACTOR=0.5

LD=${CC}

LIBBASE=${LIBBASEDIR}/libbase.a
LIBBASEDBG=${LIBBASEDIR}/libbase.dbg.a

PROTOFILES=access.c scanpref.c linkloc.c depthtab.c\
           ex2leav.c dfs.c overmax.c oversucc.c\
           addleafcount.c iterator.c streedbg.c

OBJECTS4=construct.4.o access.4.o scanpref.4.o linkloc.4.o depthtab.4.o\
          ex2leav.4.o dfs.4.o overmax.4.o oversucc.4.o\
          addleafcount.4.o iterator.4.o

OBJECTS=construct.o access.o scanpref.o linkloc.o depthtab.o\
          ex2leav.o dfs.o overmax.o oversucc.o\
          addleafcount.o iterator.o

DBGOBJECTS4=construct.dbg.4.o access.dbg.4.o streedbg.dbg.4.o\
             depthtab.dbg.4.o scanpref.dbg.4.o linkloc.dbg.4.o\
             ex2leav.dbg.4.o dfs.dbg.4.o overmax.dbg.4.o oversucc.dbg.4.o\
             addleafcount.dbg.4.o iterator.dbg.4.o

DBGOBJECTS=construct.dbg.o access.dbg.o streedbg.dbg.o depthtab.dbg.o\
             scanpref.dbg.o linkloc.dbg.o ex2leav.dbg.o dfs.dbg.o\
             overmax.dbg.o oversucc.dbg.o addleafcount.dbg.o iterator.dbg.o

SPLINTALL=access.splint scanpref.splint linkloc.splint depthtab.splint\
           ex2leav.splint dfs.splint overmax.splint oversucc.splint\
           addleafcount.splint iterator.splint streedbg.splint\
           loc.splint stree.splint

.SUFFIXES: .o .4.o .dbg.o .dbg.4.o .a .dbg.a .4.a .dbg.4.a .x .dbg.x

all: libs progs

libs: libstree.a libstree.dbg.a

progs: stree.x stree.dbg.x loc.x loc.dbg.x

include Filegoals.mf

libstree.4.a: $(OBJECTS4)
	ar sruv $@ $(OBJECTS4)

libstree.a: $(OBJECTS)
	ar sruv $@ $(OBJECTS)

libstree.dbg.4.a: $(DBGOBJECTS4)
	ar sruv $@ $(DBGOBJECTS4)

libstree.dbg.a: $(DBGOBJECTS)
	ar sruv $@ $(DBGOBJECTS)


stree.x: stree.o libstree.a
	$(LD) $(LDFLAGS) stree.o libstree.a $(LIBBASE) -o $@

loc.x: loc.o libstree.a
	$(LD) $(LDFLAGS) loc.o libstree.a $(LIBBASE) -o $@

stree.dbg.x: stree.dbg.o libstree.dbg.a 
	$(LD) $(LDFLAGS) stree.dbg.o libstree.dbg.a $(LIBBASEDBG) -o $@

loc.dbg.x: loc.dbg.o libstree.dbg.a 
	$(LD) $(LDFLAGS) loc.dbg.o libstree.dbg.a $(LIBBASEDBG) -o $@

streeproto.h: $(PROTOFILES) Mkstreeproto.sh
	@echo "make $@"
	@Mkstreeproto.sh ${PROTOFILES} > $@

.PHONY:clean
clean:splintclean
	rm -f *.[aox] *~

stree.dbg.o: streedbg.o construct.dbg.o

construct.dbg.o: construct.gen

include Dependencies.mf
