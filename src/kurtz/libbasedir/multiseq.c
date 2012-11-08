/*
  Copyright (c) 2003 by Stefan Kurtz and The Institute for
  Genomic Research.  This is OSI Certified Open Source Software.
  Please see the file LICENSE for licensing information and
  the file ACKNOWLEDGEMENTS for names of contributors to the
  code base.
*/

//\Ignore{

#include <stdio.h>
#include <stdlib.h>
#include <limits.h>
#include <string.h>
#include <ctype.h>
#include <unistd.h>
#include "types.h"
#include "debugdef.h"
#include "spacedef.h"
#include "minmax.h"
#include "protodef.h"
#include "errordef.h"
#include "chardef.h"
#include "multidef.h"

//}

/*EE 
  This file implements the data type \texttt{Multiseq}.
*/

/*
  The newline symbol
*/

#define NEWLINE  '\n'

/* 
  undefined number of database sequences
*/

#define UNDEFNUMOFDBSEQ 0

/*
  Copy a multiple component of a multiple sequence
*/

#define COPYMULTISEQ(COMP)              multiseq1->COMP = multiseq2->COMP

/*
  Store the next position where a description starts
*/

#define STORESTARTDESC\
        if(multiseq->numofsequences >= allocatedstartdesc)\
        {\
          allocatedstartdesc += 128;\
          multiseq->startdesc\
            = ALLOCSPACE(multiseq->startdesc,Uint,allocatedstartdesc);\
        }\
        multiseq->startdesc[multiseq->numofsequences]\
          = multiseq->descspace.nextfreeUchar

/*EE
  \texttt{initmultiseq} initializes a \texttt{Multiseq} record
  such that it contains no sequences.
*/

void initmultiseq(Multiseq *multiseq)
{
  multiseq->startdesc = NULL;
  INITARRAY(&multiseq->markpos,Uint);
  INITARRAY(&multiseq->descspace,Uchar);
  multiseq->sequence = NULL;
  multiseq->rcsequence = NULL;
  multiseq->numofsequences = 0;
  multiseq->totallength = 0;
}


/*EE
  The following function frees the space allocated for a \texttt{multiseq}.
*/

void freemultiseq(Multiseq *multiseq)
{
  DEBUG0(2,"# freemultiseq\n");

  if(DELETEMEMORYMAP(multiseq->markpos.spaceUint) != 0)
  {
    FREEARRAY(&multiseq->markpos,Uint);
  }
  if(DELETEMEMORYMAP(multiseq->descspace.spaceUchar) != 0)
  {
    FREEARRAY(&multiseq->descspace,Uchar);
  }
  if(DELETEMEMORYMAP(multiseq->startdesc) != 0)
  {
    FREESPACE(multiseq->startdesc);
  }
  if(multiseq->originalsequence != NULL && 
     multiseq->originalsequence != multiseq->sequence)
  {
    if(DELETEMEMORYMAP(multiseq->originalsequence) != 0)
    {
      FREESPACE(multiseq->originalsequence);
    }
  }
  if(DELETEMEMORYMAP(multiseq->sequence) != 0)
  {
    FREESPACE(multiseq->sequence);
  }
  FREESPACE(multiseq->rcsequence);
}

/*EE
  The following function applies a function \texttt{apply} to all 
  sequences in a \texttt{multiseq}. \texttt{rcmode} is \texttt{True} 
  iff the function is to be applied to the reverse complemented sequence.
  In each call, \texttt{apply} has the following arguments:
  \begin{itemize}
  \item
  the first argument is \texttt{applyinfo}
  \item
  the second argument is the number of the current sequence
  \item
  the third argument is a pointer to the start of the sequence
  \item
  the fourth argument is the length of the sequence
  \end{itemize}
*/

Sint overallsequences(BOOL rcmode,Multiseq *multiseq,void *applyinfo,
                      Sint(*apply)(void *,Uint,Uchar *,Uint))
{
  Uint i;
  Uchar *seq, *start, *end;

  if(rcmode)
  {
    seq = multiseq->rcsequence;
  } else
  {
    seq = multiseq->sequence;
  }
  for(i = 0; i < multiseq->numofsequences; i++)
  {
    if(i == 0)
    {
      start = seq;
    } else
    {
      start = seq + multiseq->markpos.spaceUint[i-1] + 1;
    }
    if(i == multiseq->numofsequences - 1)
    {
      end = seq + multiseq->totallength;
    } else
    {
      end = seq + multiseq->markpos.spaceUint[i];
    }
    DEBUG3(3,"overallsequences: i=%lu, start=%lu, end=%lu\n",
              (Showuint) i,
              (Showuint) (start-seq),(Showuint) (end-seq));
    if(apply(applyinfo,i,start,(Uint) (end - start)) != 0)
    {
      return -1;
    }
  }
  return 0;
}


/*EE
  Suppose we have a sequence of length \texttt{totalwidth}
  and a \texttt{position} in the range \([0,\texttt{totalwidth}-1]\).
  Given a sorted array of separators \texttt{recordspes} of length
  \texttt{numofrecords} such that each separator is a position in
  the range \([0,\texttt{totalwidth}-1]\). The record separator
  divides the sequence into records numbered from 0 to 
  \texttt{numofrecords}.
  The following function \texttt{getrecordnum} delivers the number of
  the record in which position occurs.
  If the position is not in the correct range, then a negative error code
  is returned. The running time of \texttt{getrecordnum} is 
  \(O(\log_{2}\texttt{numofrecords})\).
*/

Sint getrecordnum(Uint *recordseps,Uint numofrecords,Uint totalwidth,
                  Uint position)
{
  Uint *leftptr, *midptr = NULL, *rightptr, len;

  if(numofrecords == UintConst(1) || position < recordseps[0])
  {
    return 0;
  }
  if(position > recordseps[numofrecords-2])
  { 
    if(position < totalwidth)
    {
      return numofrecords - 1;
    }
    ERROR1("cannot find position %lu",(Showuint) position);
    return -1;
  }

  DEBUG2(3,"getrecordnum for pos %lu in [0..%lu]\n",
          (Showuint) position,
          (Showuint) (numofrecords-2));
  leftptr = recordseps;
  rightptr = recordseps + numofrecords - 2;
  while (leftptr<=rightptr)
  {
    len = (Uint) (rightptr-leftptr);
    midptr = leftptr + DIV2(len);
    DEBUG1(3,"left=%lu,",(Showuint) (leftptr - recordseps));
    DEBUG1(3,"mid=%lu,",(Showuint) (midptr - recordseps));
    DEBUG1(3,"right=%lu,",(Showuint) (rightptr - recordseps));
    if(*midptr < position)
    {
      if(position < *(midptr+1))
      {
        return (Sint) (midptr - recordseps + 1);
      } 
      leftptr = midptr + 1;
    } else
    {
      if(*(midptr-1) < position)
      {
        return (Sint) (midptr - recordseps);
      }
      rightptr = midptr-1;
    }
  }
  ERROR1("cannot find position %lu",(Showuint) position);
  return -1;
}

/*EE
  Given a \texttt{multiseq}, and a position in \texttt{multiseq->sequence},
  the function \texttt{getseqnum} delivers the sequence number for
  \texttt{position}. If this cannot be found, then a negative error code
  is returned. The running time of \texttt{getseqnum} is \(O(\log_{2}k)\),
  where \(k\) is the number of sequences in \texttt{multiseq}.
*/

Sint getseqnum(Multiseq *multiseq,Uint position)
{
  return getrecordnum(multiseq->markpos.spaceUint,
                      multiseq->numofsequences,
                      multiseq->totallength,
                      position);
}


/*EE
  The function \texttt{pos2pospair} computes a pair \((i,j)\) stored
  in \texttt{pos->uint0} and \texttt{pos->uint1} such that
  \(i\) is the sequence number for \texttt{position} in \texttt{multiseq},
  and \(j\) is the relative index for \texttt{position} in the \(i\)th
  sequence \(T_{i}\). If this cannot be found, then a negative error code
  is returned. In case of success, the function returns 0.
  The running time of \texttt{pos2pair} is \(O(\log_{2}k)\),
  if \(k\) is the number of sequences in \texttt{multiseq}. 
*/

Sint pos2pospair(Multiseq *multiseq,PairUint *pos,Uint position)
{
  Sint retcode;

  retcode = getseqnum(multiseq,position);
  if(retcode < 0)
  {
    return -1;
  }
  pos->uint0 = (Uint) retcode;
  if(pos->uint0 == 0)
  {
    pos->uint1 = position;
  } else
  {
    pos->uint1 = position - multiseq->markpos.spaceUint[pos->uint0-1] - 1;
  }
  return 0;
}

