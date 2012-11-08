/*
  Copyright (c) 2003 by Stefan Kurtz and The Institute for
  Genomic Research.  This is OSI Certified Open Source Software.
  Please see the file LICENSE for licensing information and
  the file ACKNOWLEDGEMENTS for names of contributors to the
  code base.
*/

#include "streedef.h"
#include "streeacc.h"

#define P3(F,A,B,C)   F(A,B,C)
#define PNO(F,A,B,C)  /* Nothing */

#define PROCESS(PL,PB)\
        while(!NILPTR(node))\
        {\
          if(ISLEAF(node))\
          {\
            leafindex = GETLEAFINDEX(node);\
            PL(processleaf,stree,leafindex,info);\
            node = LEAFBROTHERVAL(stree->leaftab[leafindex]);\
          } else\
          {\
            btptr = stree->branchtab + GETBRANCHINDEX(node);\
            PB(processbranch,stree,btptr,info);\
            node = GETBROTHER(btptr);\
          }\
        }

void oversuccsstree(Suffixtree *stree,Bref bnode,
                    void(*processleaf)(Suffixtree *,Uint,void *),
                    void(*processbranch)(Suffixtree *,Bref,void *),
                    void *info)
{
  Uint leafindex, node, *btptr;

  node = GETCHILD(bnode);
  if(processleaf != NULL)
  { 
    if(processbranch != NULL)
    {
      PROCESS(P3,P3);
    } else
    {
      PROCESS(P3,PNO);
    }
  } else
  {
    if(processbranch != NULL)
    {
      PROCESS(PNO,P3);
    } else
    {
      fprintf(stderr,"processleaf and processbranch are both undefined\n");
      exit(EXIT_FAILURE);
    }
  }
}
