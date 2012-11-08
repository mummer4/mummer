/*
  Copyright (c) 2003 by Stefan Kurtz and The Institute for
  Genomic Research.  This is OSI Certified Open Source Software.
  Please see the file LICENSE for licensing information and
  the file ACKNOWLEDGEMENTS for names of contributors to the
  code base.
*/

//\Ignore{

#ifndef MUMCAND_H
#define MUMCAND_H
#include "types.h"
#include "arraydef.h"

//}

/*
  The following structure stores MUM candidates. That is, maximal matches
  which are unique in the subject-sequence but not necessarily in the
  query sequence.
*/

typedef struct
{
  Uint mumlength,    // length of the mum
       dbstart,      // start position in the subject-sequence
       queryseq,     // number of the query sequence
       querystart;   // start position in the query sequence      
} MUMcandidate;     // \Typedef{MUMcandidate}

/*
  We store MUM-candidates in a table and hence declare a corresponding
  array.
*/

DECLAREARRAYSTRUCT(MUMcandidate);

//\Ignore{

#endif

//}
