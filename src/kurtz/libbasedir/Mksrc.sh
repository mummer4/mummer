#!/bin/sh
VSTREESRC=${DIRVSTREE}/src/vstree/src
vcopy.sh `cat Importfiles`
grep -v -f Excludemulti ${VSTREESRC}/include/multidef.h > multidef.h
gawk -f Cuthere.awk ${VSTREESRC}/kurtz/multiseq.c | grep -v -f Excludemulti > multiseq.c
Mkprotodef.sh `ls *.c` > protodef.h
for filename in `ls *.[ch]`
do
  cat ../Copyright > tmp
  Skipfirstcom.pl ${filename} >> tmp
  mv tmp ${filename}
done
