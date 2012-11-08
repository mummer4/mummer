/*
  Copyright (c) 2003 by Stefan Kurtz and The Institute for
  Genomic Research.  This is OSI Certified Open Source Software.
  Please see the file LICENSE for licensing information and
  the file ACKNOWLEDGEMENTS for names of contributors to the
  code base.
*/

//}

//\FILEINFO{streelarge.h}{Stefan Kurtz}{November 1999}

//\Ignore{

#ifndef STREELARGE_H
#define STREELARGE_H
#include <stdio.h>
#include <stdlib.h>
#include <limits.h>

/*
  This header file defines the constants and macros used for the improved
  linked list implementation technique for suffix trees as described
  in \cite{KUR:1998}. The implementation technique requires
  one integer for each leaf, two integers for each small node and
  four integers for each large node, provided
  \(n\leq 2^{27}-1=128~\emph{megabytes}\).
  For more details, see \cite{KUR:1998} and the bit layout,
  as depicted in the file \texttt{bitlayout.ps}.
*/

#define SMALLINTS           2                  // # of integers for small node
#define LARGEINTS           4                  // # of integers for large node
#define MULTBYSMALLINTS(V)  ((V) << 1)         // multiply by SMALLINTS
#define DIVBYSMALLINTS(V)   ((V) >> 1)         // div by SMALLINTS

#define DISTBITS            5                  // # bits for distance
#define DISTSHIFT           (32 - DISTBITS)    // shift for distance bits
#define MAXDISTANCE         ((1U << DISTBITS) - 1U)  // maximal distance
#define MASKDISTANCE        (((Uint) MAXDISTANCE) << DISTSHIFT) // mask dist bit

#define INDEXBITS           29     // # of bits for index value
#define LEAFBIT             1      // bit for masking a leaf address
#define NILBIT              (1U << INDEXBITS)  // bit for masking NIL value
#define MAXINDEX            (NILBIT - 1)   // MAXIMAL index value
#define EXTRAPATT           (3U << 30)     // mask extra value

#define SMALLDEPTHMARK      (1U << 31)  // mask for small depth value
#define SMALLDEPTHBITS      10          // number of bits for small depth
#define SHIFTMIDDLE         6           // shift middle part
#define SMALLDEPTH          ((1U << SMALLDEPTHBITS) - 1)  // maximal small depth
#define ISSMALLDEPTH(V)     (((V) & ~SMALLDEPTH) == 0)  // is depth small
#define MAXTLEN             ((1U << (INDEXBITS-2)) - 1)  // maximal text length
#define UNUSEDINLEAF        2                            // unused bits in leaf
#define HIGHERSIZE          (~((1U << (INDEXBITS-1-UNUSEDINLEAF)) - 1)) 
#define SHIFTHIGHER         4
#define LOWERLINKBITS       (32 - (SMALLDEPTHBITS+1))  // # bits for lower link
#define LOWERLINKSIZE       ((1U << LOWERLINKBITS) - 1)  // max value of llink
#define LOWERLINKPATT       (~(SMALLDEPTH | SMALLDEPTHMARK)) // mask llink
#define MIDDLELINKPATT      (~MAXTLEN)  // middle of llink

/*
  We use the least significant bit to discriminate references to leafs
  and branching nodes. Since base addresses are even, the unset least
  significant bit of a reference identifies a base address. For a leaf
  reference we shift the leaf number one to the right and
  set the least significant bit.
*/

#define ISLEAF(V)           ((V) & LEAFBIT)
#define ISLARGE(V)          (((V) & MASKDISTANCE) == 0)
#define MAKELEAF(V)         (((V) << 1) | LEAFBIT)
#define MAKELARGE(V)        (V)
#define MAKELARGELEAF(V)    MAKELEAF(V)
#define GETLEAFINDEX(V)     ((V) >> 1)
#define GETBRANCHINDEX(V)   (V)

#define NILPTR(P)           ((P) & NILBIT)
#define UNDEFINEDREFERENCE  (~((Uint) 0))
#define MAXTEXTLEN          ((MAXINDEX/LARGEINTS) - 3)

#define GETCHILD(B)          ((((*(B)) << 2) & MAXINDEX) | ((*((B)+1)) >> 30))
#define GETDISTANCE(B)       ((*(B)) >> DISTSHIFT)
#define GETBROTHER(B)        ((*((B)+1)) & (MAXINDEX | NILBIT))
#define GETDEPTH(B)          getdepth(B)
#define GETHEADPOS(B)        (*((B)+3) & MAXTLEN)
#define GETSUFFIXLINK(B)     getlargelinkconstruction(stree)
#define SETCHILD(B,VAL)\
        SETVAL(B,((*(B)) & MASKDISTANCE) | ((VAL) >> 2));\
        SETVAL((B)+1,((*((B)+1)) & (MAXINDEX | NILBIT)) | ((VAL) << 30))

#define SETBROTHER(B,VAL)    SETVAL((B) + 1,((*((B)+1)) & EXTRAPATT) | (VAL))
#define SETDISTANCE(B,VAL)\
        SETVAL(B,((*(B)) & (~MASKDISTANCE)) | ((VAL) << DISTSHIFT))
#define SETDEPTHHEADPOS(DP,HN) setdepthheadposition(stree,DP,HN)

#define SETNEWCHILD(B,VAL)     SETVAL(B,VAL)
#define SETNEWCHILDBROTHER(CH,BR)\
        SETVAL(stree->nextfreebranch,(CH) >> 2);\
        SETVAL(stree->nextfreebranch+1,((CH) << 30) | (BR))

#define SETSUFFIXLINK(SL)   setsuffixlink(stree,SL)

#define LEAFBROTHERVAL(V)      ((V) & (MAXINDEX | NILBIT))
#define SETLEAFBROTHER(B,VAL)  *(B) = (*(B) & EXTRAPATT) | (VAL)

#define GETCHAINEND(C,B,D)     C = (B) + (MULTBYSMALLINTS(D))
#define MAKEBRANCHADDR(V)      (V)
#define SETBRANCHNODEOFFSET    /* nothing */

//\Ignore{

#ifdef DEBUG
#define CHILDREFERSTOLEAF(B)   ISLEAF(*(B))
#endif

//}

/*
  This following function looks up the depth of a large node.
  \emph{btptr} refers to the first integer of that node.
*/

static Uint getdepth(Uint *btptr)
{
  Uint thirdval = *(btptr+2);

  if(thirdval & SMALLDEPTHMARK)
  {
    return thirdval & SMALLDEPTH;
  }
  return thirdval & MAXTLEN;
}

//\Ignore{

#if SYMBOLBYTES == 1
#define LARGESTCHARINDEX          UCHAR_MAX
#else
#define LARGESTCHARINDEX          stree->lastcharindex
#endif

#endif

//}
