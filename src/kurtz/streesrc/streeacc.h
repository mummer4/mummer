/*
  Copyright (c) 2003 by Stefan Kurtz and The Institute for
  Genomic Research.  This is OSI Certified Open Source Software.
  Please see the file LICENSE for licensing information and
  the file ACKNOWLEDGEMENTS for names of contributors to the
  code base.
*/

#ifndef STREEACC_H
#define STREEACC_H

#ifdef STREESMALL
#include "streesmall.h"
#endif

#ifdef STREELARGE
#include "streelarge.h"
#endif

#ifdef STREEHUGE
#include "streehuge.h"
#endif

#ifdef DEBUG
#define SHOWVAL(S)    fprintf(stderr,"#%s %lu\n",#S,(Showuint) S)
#define SETVAL(E,VAL) *(E) = VAL;\
                      if((E) > stree->maxset)\
                      {\
                        stree->maxset = E;\
                      }
#else

//}

/*
  This file contains some macros for retrieving depth, headpositions,
  and suffix links.
*/

#define SETVAL(E,VAL) *(E) = VAL

//\Ignore{

#endif

//}

/*
  \texttt{GETBOTH} retrieves the \emph{depth} and the \emph{headposition} of 
  a branching node referred to by \texttt{PT}. In case, we need these values
  for a node of the current chain, the distance is not set. So we compute 
  it as the difference between the next free base address, and the base 
  address of the node the chain starts with. Then we refer to the current depth
  and the number of the current leaf. In case, the node is large, we can
  directly look up the values. In case, the node is small, we determine
  the distance, and a pointer to the large node at the end of the chain.
  Then we can retrieve the depth and the head positions from this, as
  proved in \cite{KUR:1998}, Observation 7.
*/

#define GETBOTH(DP,HP,PT) \
        if(stree->chainstart != NULL && (PT) >= stree->chainstart)\
        {\
          distance = 1 + \
                     DIVBYSMALLINTS((Uint) (stree->nextfreebranch - (PT)));\
          DP = stree->currentdepth + distance;\
          HP = stree->nextfreeleafnum - distance;\
        } else\
        {\
          if(ISLARGE(*(PT)))\
          {\
            DP = GETDEPTH(PT);\
            HP = GETHEADPOS(PT);\
          } else\
          {\
            distance = GETDISTANCE(PT);\
            GETCHAINEND(largeptr,PT,distance);\
            DP = GETDEPTH(largeptr) + distance;\
            HP = GETHEADPOS(largeptr) - distance;\
          }\
        }

/*
  The macros \texttt{GETONLYHEADPOS}, \texttt{GETONLYDEPTH}, and
  \texttt{GETDEPTHAFTERHEADPOS} retrieve the depth or the head position.
  This is done as in the previous macro, and we omit it here.
*/

//\Ignore{

#define GETONLYHEADPOS(HP,PT) \
        if(stree->chainstart != NULL && (PT) >= stree->chainstart)\
        {\
          distance = 1 + DIVBYSMALLINTS((Uint) (stree->nextfreebranch - (PT)));\
          HP = stree->nextfreeleafnum - distance;\
        } else\
        {\
          if(ISLARGE(*(PT)))\
          {\
            HP = GETHEADPOS(PT);\
          } else\
          {\
            distance = GETDISTANCE(PT);\
            GETCHAINEND(largeptr,PT,distance);\
            HP = GETHEADPOS(largeptr) - distance;\
          }\
        }

#define GETONLYDEPTH(DP,PT) \
        if(stree->chainstart != NULL && (PT) >= stree->chainstart)\
        {\
          distance = 1 + DIVBYSMALLINTS((Uint) (stree->nextfreebranch - (PT)));\
          DP = stree->currentdepth  + distance;\
        } else\
        {\
          if(ISLARGE(*(PT)))\
          {\
            DP = GETDEPTH(PT);\
          } else\
          {\
            distance = GETDISTANCE(PT);\
            GETCHAINEND(largeptr,PT,distance);\
            DP = GETDEPTH(largeptr) + distance;\
          }\
        }

#define GETDEPTHAFTERHEADPOS(DP,PT) \
        if(stree->chainstart != NULL && (PT) >= stree->chainstart)\
        {\
          DP = stree->currentdepth + distance;\
        } else\
        {\
          if(ISLARGE(*(PT)))\
          {\
            DP = GETDEPTH(PT);\
          } else\
          {\
            DP = GETDEPTH(largeptr) + distance;\
          }\
        }

#define GETHEADPOSAFTERDEPTH(HP,PT) \
        if(stree->chainstart != NULL && (PT) >= stree->chainstart)\
        {\
          HP = stree->nextfreeleafnum - distance;\
        } else\
        {\
          if(ISLARGE(*(PT)))\
          {\
            HP = GETHEADPOS(PT);\
          } else\
          {\
            HP = GETHEADPOS(largeptr) - distance;\
          }\
        }

#define NEXTNODE(PT)\
        if(ISLARGE(*(PT)))\
        {\
          PT += LARGEINTS;\
        } else\
        {\
          PT += SMALLINTS;\
        }

//}

/*
  The suffix link is always determined for the \emph{headnode}. If this
  is large, the we have to retrieve it from that node. Otherwise, the
  suffix link refers to the next node. In both cases, the depth of the 
  \emph{headnode} is decremented.
*/

#define FOLLOWSUFFIXLINK\
        if(ISLARGE(*(stree->headnode)))\
        {\
          stree->headnode = stree->branchtab + GETSUFFIXLINK(stree->headnode);\
        } else\
        {\
          stree->headnode += SMALLINTS;\
        }\
        stree->headnodedepth--

/*
  Whenever \emph{insertleaf} is called, \emph{onsuccpath} stores the 
  address of the new leaf.
  Whenever \emph{insertbranch} is called, \emph{onsuccpath} stores the 
  address of the new branching node. Both nodes are a successor of the
  node, for which a suffix link is possible to be computed in the
  next step. In case linear retrieval of suffix links is required, 
  it is possible to start at the node referenced by \emph{onsuccpath}.
*/

#if defined(STREELARGE) || defined(STREESMALL)
#define RECALLSUCC(S)             stree->onsuccpath = S
#else
#define RECALLSUCC(S)             /* Nothing */
#endif

/*
  The following three macros handle the setting of the suffix link in a 
  nil reference. The \emph{RECALL}-macros store the address of the 
  reference. In case the reference is a new leaf, this is marked. 
*/

#define RECALLNEWLEAFADDRESS(A)   stree->setlink = A;\
                                  stree->setatnewleaf = True
#define RECALLLEAFADDRESS(A)      stree->setlink = A;\
                                  stree->setatnewleaf = False
#define RECALLBRANCHADDRESS(A)    stree->setlink = (A) + 1;\
                                  stree->setatnewleaf = False

#ifdef STREEHUGE
#define SETNILBIT                 *(stree->setlink) = NILBIT
#else
#define SETNILBIT                 if(stree->setatnewleaf)\
                                  {\
                                    *(stree->setlink) = NILBIT;\
                                  } else\
                                  {\
                                    *(stree->setlink) |= NILBIT;\
                                  }
#endif

#define SETMAXBRANCHDEPTH(D)      if((D) > stree->maxbranchdepth)\
                                  {\
                                    stree->maxbranchdepth = D;\
                                  }

//\Ignore{

/*
#ifdef SHOWLEAD
#define LEADLEVEL 3
#else
#define LEADLEVEL 4
#endif
*/

#define LEADLEVEL 2

#ifdef DEBUG
#define SHOWINDEX(NODE)\
        if((NODE) == UNDEFINEDREFERENCE)\
        {\
          DEBUG0(LEADLEVEL,"UNDEFINEDREFERENCE");\
        } else\
        {\
          if(NILPTR(NODE))\
          {\
            DEBUG0(LEADLEVEL,"NILPTR");\
          } else\
          {\
            if(ISLEAF(NODE))\
            {\
              DEBUG1(LEADLEVEL,"Leaf %lu",(Showuint) GETLEAFINDEX(NODE));\
            } else\
            {\
              if(ISLARGE(stree->branchtab[GETBRANCHINDEX(NODE)]))\
              {\
                DEBUG1(LEADLEVEL,"Large %lu",(Showuint) GETBRANCHINDEX(NODE));\
              } else\
              {\
                DEBUG1(LEADLEVEL,"Small %lu",(Showuint) NODE);\
              }\
            }\
          }\
        }
#else
#define SHOWINDEX(NODE) /* Nothing */
#endif

#ifdef DEBUG
void showtable(Suffixtree *stree,BOOL final);
void checkstree(Suffixtree *stree);
void showstate(Suffixtree *stree);
void showstree(Suffixtree *stree);
void enumlocations(Suffixtree *stree,void(*processloc)(Suffixtree *stree,Location *));
void checklocation(Suffixtree *stree,Location *loc);
#endif

#endif

//}
