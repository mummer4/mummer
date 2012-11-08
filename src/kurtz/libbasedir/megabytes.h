/*
  Copyright (c) 2003 by Stefan Kurtz and The Institute for
  Genomic Research.  This is OSI Certified Open Source Software.
  Please see the file LICENSE for licensing information and
  the file ACKNOWLEDGEMENTS for names of contributors to the
  code base.
*/

//\Ignore{

#ifndef MEGABYTES_H
#define MEGABYTES_H
#include "types.h"

//}

/*
  The following macro transforms bytes into megabytes.
*/

#define MEGABYTES(V)  ((double) (V)/((UintConst(1) << 20) - 1))

//\Ignore{

#endif

//}
