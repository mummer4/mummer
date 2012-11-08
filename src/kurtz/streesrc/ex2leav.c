/*
  Copyright (c) 2003 by Stefan Kurtz and The Institute for
  Genomic Research.  This is OSI Certified Open Source Software.
  Please see the file LICENSE for licensing information and
  the file ACKNOWLEDGEMENTS for names of contributors to the
  code base.
*/

#include "debugdef.h"
#include "streedef.h"
#include "streeacc.h"

BOOL exactlytwoleavesstree(Suffixtree *stree,PairUint *twoleaves,Bref start)
{
  BOOL firstleaffound = False;
  Uint tmpval, node;

  node = GETCHILD(start);
  while(True)
  {
    if(ISLEAF(node))
    {
      if(firstleaffound)
      {
        twoleaves->uint1 = GETLEAFINDEX(node);
        if(twoleaves->uint0 > twoleaves->uint1)
        {
          tmpval = twoleaves->uint1;
          twoleaves->uint1 = twoleaves->uint0;
          twoleaves->uint0 = tmpval;
        }
        DEBUG2(3,"has two leafs: %lu %lu\n",(Showuint) twoleaves->uint0,
                                            (Showuint) twoleaves->uint1);
        if(NILPTR(LEAFBROTHERVAL(stree->leaftab[GETLEAFINDEX(node)])))
        {
          return True;
        }
        return False;
      } 
      twoleaves->uint0 = GETLEAFINDEX(node);
      DEBUG1(3,"first successor is leaf %lu\n",(Showuint) twoleaves->uint0);
      firstleaffound = True;
      node = LEAFBROTHERVAL(stree->leaftab[GETLEAFINDEX(node)]);
    } else
    {
      DEBUG0(3,"has branch successor\n");
      return False;
    }
  }
}
