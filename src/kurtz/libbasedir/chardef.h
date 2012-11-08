/*
  Copyright (c) 2003 by Stefan Kurtz and The Institute for
  Genomic Research.  This is OSI Certified Open Source Software.
  Please see the file LICENSE for licensing information and
  the file ACKNOWLEDGEMENTS for names of contributors to the
  code base.
*/

//\Ignore{

#ifndef CHARDEF_H
#define CHARDEF_H
#include <limits.h>

//}

/*
  This file defines some character values used when storing
  multiple sequences.
*/

#define SEPARATOR       UCHAR_MAX         // separator symbol in multiple seq
#define WILDCARD        (SEPARATOR-1)     // wildcard symbol in multiple seq
#define UNDEFCHAR       (SEPARATOR-2)     // undefined character in multiple seq
#define ISSPECIAL(C)    ((C) >= WILDCARD) // either WILDCARD or SEPARATOR
#define ISNOTSPECIAL(C) ((C) < WILDCARD)  // neither WILDCARD nor SEPARATOR

//\Ignore{

#endif

//}
