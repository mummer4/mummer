/*
  Copyright (c) 2003 by Stefan Kurtz and The Institute for
  Genomic Research.  This is OSI Certified Open Source Software.
  Please see the file LICENSE for licensing information and
  the file ACKNOWLEDGEMENTS for names of contributors to the
  code base.
*/

//}

//\FILEINFO{streesmall.h}{Stefan Kurtz}{November 1999}

//\Ignore{

#ifndef STREESTDEF_H
#define STREESTDEF_H
#include <stdio.h>
#include <stdlib.h>
#include <limits.h>

//}

/*
  This header file defines the constants and macros used for the improved 
  linked list implementation technique for suffix trees as described
  in \cite{KUR:BAL:1999}. The implementation technique requires
  one integer for each leaf, two integers for each small node and
  3 integers for each large node, provided 
  \(n\leq 2^{21}-1=2~\emph{megabytes}\).
  For more details, see \cite{KUR:1998,KUR:BAL:1999} and the bit layout,
  as depicted in the file \texttt{bitlayout.ps}.
*/

#define SMALLINTS           2
#define LARGEINTS           3
#define MULTBYSMALLINTS(V)  ((V) << 1)
#define DIVBYSMALLINTS(V)   ((V) >> 1)

#define DISTBITS            16
#define MAXDISTANCE         ((1 << DISTBITS) - 1)

#define INDEXBITS           23
#define NILBIT              (1 << INDEXBITS)
#define LARGEBIT            (1 << INDEXBITS)
#define MAXINDEX            (NILBIT - 1)

#define ISLEAF(V)           ((V) < stree->branchnodeoffset)
#define ISLARGE(V)          ((V) & LARGEBIT)
#define MAKELEAF(V)         (V)
#define MAKELARGE(V)        (V)
#define MAKELARGELEAF(V)    (V)
#define MAKEBRANCHADDR(V)   ((V) + stree->branchnodeoffset)
#define GETLEAFINDEX(V)     (V)
#define GETBRANCHINDEX(V)   ((V) - stree->branchnodeoffset)

#define SMALLDEPTHBITS      11
#define SMALLDEPTHMARK      (1 << (INDEXBITS+1))
#define SMALLDEPTH          ((1 << SMALLDEPTHBITS) - 1)
#define ISSMALLDEPTH(V)     (((V) & ~SMALLDEPTH) == 0)

#define NILPTR(P)           ((P) & NILBIT)
#define UNDEFINEDREFERENCE  (~((Uint) 0))
#define MAXTEXTLEN          2097150

#define GETCHILD(B)         ((*(B)) & MAXINDEX)
#define SETCHILD(B,VAL)     SETVAL(B,((*(B)) & (511 << 23)) | (VAL))

#define GETDISTANCE(B)      (((*(B)) >> 24) | (((*((B)+1)) >> 16) & (255 << 8)))

#define SETDISTANCE(B,VAL)\
        SETVAL(B,((*(B)) & ((1 << 23) - 1)) | ((VAL) << 24));\
        SETVAL(B+1,((*((B)+1)) & ((1 << 24) - 1)) | (((VAL) & (255 << 8)) << 16))

#define GETBROTHER(B)        ((*((B)+1)) & (NILBIT | MAXINDEX))
#define SETBROTHER(B,VAL)    SETVAL((B)+1,((*((B)+1)) & (255 << 24)) | (VAL))

#define GETDEPTH(B)          getdepth(B)
#define GETHEADPOS(B)        ((*((B)+2)) & ((1 << 21) - 1))

#define SETDEPTHHEADPOS(DP,HD)  setdepthheadposition(stree,DP,HD)

#define GETSUFFIXLINK(B)     getlargelinkconstruction(stree)

#define SETNEWCHILD(B,CH)   SETVAL(B,(CH) | LARGEBIT)
#define SETNEWCHILDBROTHER(CH,BR)\
        SETVAL(stree->nextfreebranch,(CH) | LARGEBIT);\
        SETVAL(stree->nextfreebranch+1,BR)

#define SETSUFFIXLINK(SL)   setsuffixlink(stree,SL)

#define LEAFBROTHERVAL(V)      ((V) & (MAXINDEX | NILBIT))
#define SETLEAFBROTHER(B,VAL)  *(B) = (*(B) & (255 << 24)) | (VAL)

#define GETCHAINEND(C,B,D)     C = (B) + (MULTBYSMALLINTS(D)) 
#define SETBRANCHNODEOFFSET    stree->branchnodeoffset = stree->textlen + 1

//\Ignore{

#ifdef DEBUG
#define CHILDREFERSTOLEAF(B)  ISLEAF(*(B))
#endif

//}

/*
  This following function looks up the depth of a large node. 
  \emph{btptr} refers to the first integer of that node.
*/

static Uint getdepth(Uint *btptr)
{
  if(*(btptr+1) & SMALLDEPTHMARK)
  {
    return *(btptr+2) >> 21;
  }
  return ((*btptr &     (7 << 29)) >> 11) |
         ((*(btptr+1) & (127 << 25)) >> 14) |
         (*(btptr+2) >> 21);
}

//\Ignore{

#if SYMBOLBYTES == 1
#define LARGESTCHARINDEX          UCHAR_MAX
#else
#define LARGESTCHARINDEX          stree->lastcharindex
#endif


#endif

//}
