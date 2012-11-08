/*
  Copyright (c) 2003 by Stefan Kurtz and The Institute for
  Genomic Research.  This is OSI Certified Open Source Software.
  Please see the file LICENSE for licensing information and
  the file ACKNOWLEDGEMENTS for names of contributors to the
  code base.
*/

#ifndef SYMBOLDEF_H
#define SYMBOLDEF_H

#ifndef SYMBOLBYTES
#define SYMBOLBYTES 1
#endif

#if SYMBOLBYTES == 1
typedef Uchar SYMBOL;
#else
#if SYMBOLBYTES == 2
typedef Ushort SYMBOL;
#else
#if SYMBOLBYTES == 4
typedef Uint SYMBOL;
#endif
#endif
#endif

#endif
