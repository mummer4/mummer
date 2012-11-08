#!/bin/sh
if test $# -lt 1
then
  echo "Usage: $0 <cfilenamelist>"
  exit 1
fi
cat << ENDOFINCLUDE
/* 
  This file is generated. Do not edit.

  A Library for the Efficient Construction and Application of Suffix Trees

  Copyright (c) 2003 by Stefan Kurtz and The Institute for
  Genomic Research.  This is OSI Certified Open Source Software.
  Please see the file LICENSE for licensing information and
  the file ACKNOWLEDGEMENTS for names of contributors to the
  code base.
*/

#ifndef STREEPROTO_H
#define STREEPROTO_H

#ifdef __cplusplus
extern "C" {
#endif

ENDOFINCLUDE

grep -h '^#define CONSTRUCT' construct.c |\
    sed -e 's/^#define CONSTRUCT \(.*\)/\1;/'

cat << ENDOFINCLUDE

void freestree(Suffixtree *stree);
#ifdef __cplusplus
}
#endif

ENDOFINCLUDE

cproto -I ../libbasedir "$@"

echo "#endif"
