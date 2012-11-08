#-- Imported variables from top level makefile
# BIN_DIR AUX_BIN_DIR CXX CC CFLAGS CXXFLAGS LDFLAGS

##include ../Makedef

ifndef BIN_DIR
BIN_DIR := $(CURDIR)
endif

VPATH := $(BIN_DIR)

ALL := maxmat3.x maxmat3.dbg.x

LIBBASEDIR=../libbasedir
LIBSTREEDIR=../streesrc
INCLUDEDIR=-I${LIBBASEDIR} -I${LIBSTREEDIR}

override CFLAGS+=$(INCLUDEDIR)
##CFLAGS=${DEFINECFLAGS} $(INCLUDEDIR)
##LDFLAGS=${DEFINELDFLAGS}
SPLINTFLAGS=${INCLUDEDIR} -f ../Splintoptions -DDEBUG
LD=$(CC)

LIBBASE=$(LIBBASEDIR)/libbase.a
LIBBASEDBG=$(LIBBASEDIR)/libbase.dbg.a
LIBSTREE=$(LIBSTREEDIR)/libstree.a
LIBSTREEDBG=$(LIBSTREEDIR)/libstree.dbg.a

MUM3OBJECTS=maxmat3.o maxmatopt.o maxmatinp.o\
            findmumcand.o findmaxmat.o procmaxmat.o

MUM3DBGOBJECTS=maxmat3.dbg.o maxmatopt.dbg.o maxmatinp.dbg.o\
               findmumcand.dbg.o findmaxmat.dbg.o procmaxmat.dbg.o

SPLINTALL=maxmat3.splint maxmatopt.splint maxmatinp.splint\
          findmumcand.splint findmaxmat.splint procmaxmat.splint

.SUFFIXES: .dbg.o

all: $(ALL)

mummer: $(MUM3OBJECTS) $(LIBSTREE)
	$(LD) $(LDFLAGS) $(MUM3OBJECTS) $(LIBSTREE) $(LIBBASE) \
	-o $(BIN_DIR)/$@; chmod 755 $(BIN_DIR)/$@

maxmat3.x: $(MUM3OBJECTS) $(LIBSTREE)
	$(LD) $(LDFLAGS) $(MUM3OBJECTS) $(LIBSTREE) $(LIBBASE) \
	-o $(BIN_DIR)/$@; chmod 755 $(BIN_DIR)/$@

maxmat3.dbg.x: ${MUM3DBGOBJECTS} $(LIBSTREEDBG)
	$(LD) $(LDFLAGS) $(MUM3DBGOBJECTS) $(LIBSTREEDBG) $(LIBBASEDBG) \
	-lm -o $(BIN_DIR)/$@; chmod 755 $(BIN_DIR)/$@

include Filegoals.mf

.PHONY:clean
clean:splintclean
	rm -f *.[ox] shit[123] *~
	cd $(BIN_DIR); rm -f $(ALL) mummer

include Dependencies.mf
