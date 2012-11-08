/*
  Copyright (c) 2003 by Stefan Kurtz and The Institute for
  Genomic Research.  This is OSI Certified Open Source Software.
  Please see the file LICENSE for licensing information and
  the file ACKNOWLEDGEMENTS for names of contributors to the
  code base.
*/

#include "types.h"
#include "protodef.h"
#include "debugdef.h"
#include "streedef.h"
#include "spacedef.h"
#include "arraydef.h"

typedef struct
{
  Suffixtree *stree;                      // suffix tree info
  ArrayUint countstack;
} Countstate;

static Sint processleaf(/*@unused@*/ Uint leafnumber,
                        /*@unused@*/ Bref lcpnode,void *info)
{
  Countstate *cstate = (Countstate *) info;

  cstate->countstack.spaceUint[cstate->countstack.nextfreeUint - 1]++;
  return 0;
}

static BOOL processbranch1(/*@unused@*/ Bref nodeptr,void *info)
{
  Countstate *cstate = (Countstate *) info;

  CHECKARRAYSPACE(&(cstate->countstack),Uint,128);
  cstate->countstack.spaceUint[cstate->countstack.nextfreeUint++] = 0;
  return True;
}

#ifdef REPNUM
static Uint tmpleafcount, repnum = 0;
#endif

static Sint processbranch2(Bref nodeptr,void *info)
{
  Countstate *cstate = (Countstate *) info;
  Branchinfo branchinfo;
  Uint *father;

  cstate->countstack.nextfreeUint--;
  father = cstate->countstack.spaceUint + cstate->countstack.nextfreeUint - 1;
  *father += *(father+1);
  getbranchinfostree(cstate->stree,ACCESSDEPTH | ACCESSHEADPOS,
                     &branchinfo,nodeptr);
#ifdef REPNUM
  if(branchinfo.depth >= 12)
  {
    tmpleafcount = *(father+1);
    repnum += ((tmpleafcount * tmpleafcount) >> 1);
  }
#endif
  cstate->stree->leafcounts[branchinfo.headposition] = *(father+1);
  DEBUG2(3,"leafcount(%lu) = %lu\n",
            (Showuint) BRADDR2NUM(cstate->stree,nodeptr),
            (Showuint) *(father+1));
  *(father+1) = 0;
  return 0;
}

/* 
   depth first traversal of the tree.
   push branch node on stack when visited first.
   suppose there is an edge \(father -> son\) such that son is the root of
   a subtree which has been processed. 
   If son is leaf then call processleaf(father,son). 
   If son is not a leaf, then call processbranch(father,son)
*/

Uint getleafcountstree(Suffixtree *stree,Bref nodeptr)
{
  Branchinfo branchinfo;

#ifdef DEBUG
  if(stree->leafcounts == NULL)
  {
    fprintf(stderr,"leafcounts are not stored yet. First call \"addleafcountsstree\"\n");
    exit(EXIT_FAILURE);
  }
#endif
  getbranchinfostree(stree,ACCESSHEADPOS,&branchinfo,nodeptr);
  return stree->leafcounts[branchinfo.headposition];
}

#ifdef DEBUG

static Sint countleafs(/*@unused@*/ Uint leafindex,
                       /*@unused@*/ Bref lcanode,void *info)
{
  Uint *lc = (Uint *) info;

  (*lc)++;
  return 0;
}

static void checkleafcount(Suffixtree *stree,
                           Bref nodeptr,
                           /*@unused@*/ Uint depth,
                           Uint headposition,
                           /*@unused@*/ void *info)
{
  Reference ref;
  Uint lc = 0;
  
  DEBUG1(3,"checkleafcount of %lu\n",(Showuint) BRADDR2NUM(stree,nodeptr));
  ref.toleaf = False;
  ref.address = nodeptr;
  if(depthfirststree(stree,&ref,countleafs,NULL,NULL,
                     NULL,NULL,(void *) &lc) != 0)
  {
    fprintf(stderr,"depthfirststree failed\n");
    exit(EXIT_FAILURE);
  }
  if(stree->leafcounts[headposition] != lc)
  {
    fprintf(stderr,"leafcount[%lu] = %lu != %lu lc\n",
                    (Showuint) BRADDR2NUM(stree,nodeptr),
                    (Showuint) stree->leafcounts[headposition],
                    (Showuint) lc);
    exit(EXIT_FAILURE);
  }
}

static void checkleafcountall(Suffixtree *stree)
{
  overallstree(stree,True,checkleafcount,NULL);
}
#endif

Sint addleafcountsstree(Suffixtree *stree)
{
  Countstate cstate;
  Reference rootref;
  Uint i;

  cstate.stree = stree;
  INITARRAY(&(cstate.countstack),Uint);
  CHECKARRAYSPACE(&(cstate.countstack),Uint,128);
  stree->leafcounts = ALLOCSPACE(NULL,Uint,stree->nextfreeleafnum+1);
  for(i=0; i<=stree->nextfreeleafnum; i++)
  {
    stree->leafcounts[i] = 0;
  }
  cstate.countstack.spaceUint[cstate.countstack.nextfreeUint++] = 0;
  rootref.toleaf = False;
  rootref.address = ROOT(stree);
  if(depthfirststree(stree,&rootref,processleaf,processbranch1,
                     processbranch2,NULL,NULL,(void *) &cstate) != 0)
  {
    return -1;
  }
  DEBUGCODE(1,checkleafcountall(stree));
  FREEARRAY(&(cstate.countstack),Uint);
#ifdef REPNUM
  printf("repnum=%lu\n",(Showuint) repnum);
#endif
  return 0;
}
