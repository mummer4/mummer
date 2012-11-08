/*
  Copyright (c) 2003 by Stefan Kurtz and The Institute for
  Genomic Research.  This is OSI Certified Open Source Software.
  Please see the file LICENSE for licensing information and
  the file ACKNOWLEDGEMENTS for names of contributors to the
  code base.
*/

#include "arraydef.h"
#include "spacedef.h"
#include "debugdef.h"
#include "streedef.h"
#include "streeacc.h"
#include "protodef.h"

#define SETCURRENT(V)\
        if(ISLEAF(V))\
        {\
          current->address = stree->leaftab + GETLEAFINDEX(V);\
          current->toleaf = True;\
        } else\
        {\
          current->address = stree->branchtab + GETBRANCHINDEX(V);\
          current->toleaf = False;\
          dfsstate->secondtime = False;\
        }

/*
  compute the first branching node of the suffix tree.
*/

/*@null@*/ Bref firstbranchingnode(Suffixtree *stree)
{
  if(stree->branchtab >= stree->nextfreebranch)
  {
    return NULL;
  } 
  return stree->branchtab;
}

/*
  Given branching node \emph{bptr}, compute the branching node following
  \emph{bptr} in table \emph{branchtab}. If \emph{bptr} is 
  the last branching node, then return \texttt{NULL}.
*/

/*@null@*/ Bref nextbranchingnode(Suffixtree *stree,Bref bptr)
{
  Bref nodeptr;

  if(bptr < stree->branchtab || bptr >= stree->nextfreebranch)
  {
    return NULL;
  }
  if(ISLARGE(*bptr))
  {
    nodeptr = bptr + LARGEINTS;
  } else
  {
    nodeptr = bptr + SMALLINTS;
  }
  if(nodeptr >= stree->nextfreebranch)
  {
    return NULL;
  } 
  return nodeptr;
}

/*
  Compute the first leaf of the suffix tree.
*/

Lref firstleaf(Suffixtree *stree)
{
  return stree->leaftab;
}

/*
  Given leaf \emph{lptr}, compute the leaf following 
  \emph{lptr} in table \emph{branchtab}. If \emph{lptr} is the last
  leaf, then return \texttt{NULL}.
*/

/*@null@*/ Lref nextleaf(Suffixtree *stree,Lref lptr)
{
  if(lptr < stree->leaftab || lptr >= stree->leaftab + stree->textlen)
  {
    return NULL;
  }
  return lptr+1;
}

/*
  Compute the first node in the suffix tree, which is a branching node if
  it exists.
*/

/*@null@*/ Reference *firstnode(Suffixtree *stree,Reference *refspace)
{
  if(stree->branchtab >= stree->nextfreebranch)
  {
    return NULL;
  } 
  refspace->toleaf = False;
  refspace->address = stree->branchtab;
  return refspace;
}

/*
  Given a node referenced by \emph{nref}, compute a reference
  to the node following \emph{nref}. If \emph{nref} is the last
  leaf, then return \texttt{NULL}.
*/

/*@null@*/ Reference *nextnode(Suffixtree *stree,Reference *nref,
                               Reference *refspace)
{
  Lref lptr, bptr;

  if(nref->toleaf)
  {
    lptr = nextleaf(stree,nref->address);
    if(lptr == NULL)
    {
      return NULL;
    }
    refspace->toleaf = True;
    refspace->address = lptr;
  } else
  {
    bptr = nextbranchingnode(stree,nref->address);
    if(bptr == NULL)
    {
      refspace->toleaf = True;
      refspace->address = firstleaf(stree);
      return NULL;
    }
    refspace->toleaf = False;
    refspace->address = bptr;
  }
  return refspace;
}

static void int2ref(Suffixtree *stree,Reference *ref,Uint i)
{
  if(ISLEAF(i))
  {
    ref->toleaf = True;
    ref->address = stree->leaftab + GETLEAFINDEX(i);
  } else
  {
    ref->toleaf = False;
    ref->address = stree->branchtab + GETBRANCHINDEX(i);
  }
}

/*
  Compute the first successor of a given branching node. 
*/

/*@null@*/ Reference *firstsucc(Suffixtree *stree,Bref bptr,
                                Reference *refspace)
{
  if(bptr < stree->branchtab || bptr >= stree->nextfreebranch)
  {
    return NULL;
  }
  int2ref(stree,refspace,GETCHILD(bptr));
  return refspace;
}

/*
  Compute the right brother of a given branching node. If there is none,
  then return \texttt{NULL}.
*/

/*@null@*/ Reference *rightbrother(Suffixtree *stree,Reference *node)
{
  Uint brotherval;

  if(node->toleaf)
  {
    brotherval = LEAFBROTHERVAL(*(node->address));
  } else
  {
    brotherval = GETBROTHER(node->address);
  }
  if(NILPTR(brotherval))
  {
    return NULL;
  }  
  int2ref(stree,node,brotherval);
  return node;
}

/*
  Compute the first node in a depth first traversal of the suffix tree.
*/

Reference *firstnodedfs(Suffixtree *stree,DFSstate *dfsstate,
                        Reference *current)
{
  if(!current->toleaf)
  {
    dfsstate->secondtime = False;
    INITARRAY(&(dfsstate->stack),Bref);
    STOREINARRAY(&(dfsstate->stack),Bref,128,current->address);
    SETCURRENT(GETCHILD(current->address));
  } 
  return current;
}

/*
  Given a current reference and the current state of a depth first
  order traversal, compute the next node, and update the state accordingly.
*/

/*@null@*/ Reference *nextnodedfs(Suffixtree *stree,Reference *current,
                                  DFSstate *dfsstate)
{
  Uint child, brotherval;

  if(current->toleaf)
  {
    if(dfsstate->stack.nextfreeBref == 0)
    {
      return NULL;
    }
    brotherval = LEAFBROTHERVAL(*(current->address));
    if(NILPTR(brotherval))
    {
      if(dfsstate->stack.nextfreeBref == UintConst(1))
      {
        INITARRAY(&(dfsstate->stack),Bref);
        return NULL;
      }
      (dfsstate->stack.nextfreeBref)--;
      current->address = dfsstate->stack.spaceBref[dfsstate->stack.nextfreeBref];
      current->toleaf = False;
      dfsstate->secondtime = True;
    } else
    {
      SETCURRENT(brotherval);
    }
  } else
  {
    if(dfsstate->secondtime)
    {
      brotherval = GETBROTHER(current->address);
      if(NILPTR(brotherval))
      {
        if(dfsstate->stack.nextfreeBref == UintConst(1))
        {
          INITARRAY(&(dfsstate->stack),Bref);
          return NULL;
        }
        (dfsstate->stack.nextfreeBref)--;
        current->address = dfsstate->stack.spaceBref[dfsstate->stack.nextfreeBref];
        current->toleaf = False;
        dfsstate->secondtime = True;
      } else
      {
        SETCURRENT(brotherval);
      }
    } else
    {
      STOREINARRAY(&dfsstate->stack,Bref,128,current->address);
      DEBUG1(3,"#push[%lu]=",(Showuint) (dfsstate->stack.nextfreeBref-1));
      DEBUG1(3,"%lu\n",(Showuint) BRADDR2NUM(stree,current->address));
      child = GETCHILD(current->address);
      SETCURRENT(child);
    }
  }
  return current;
}
