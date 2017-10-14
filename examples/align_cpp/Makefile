CC = g++
# Uncomment to use openmp
# OPENMP = -fopenmp
CXXFLAGS = -Wall -Werror -Wno-unknown-pragmas $(shell pkg-config --cflags mummer) $(OPENMP)
LDFLAGS = $(shell pkg-config --libs-only-L mummer) $(shell pkg-config --libs-only-L mummer | sed 's/-L/-Wl,-rpath,/') $(OPENMP)
LDLIBS = $(shell pkg-config --libs-only-l mummer)

all: align
align: align.o

clean:
	rm -f align.o align
