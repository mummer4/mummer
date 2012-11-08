/*
  Copyright (c) 2003 by Stefan Kurtz and The Institute for
  Genomic Research.  This is OSI Certified Open Source Software.
  Please see the file LICENSE for licensing information and
  the file ACKNOWLEDGEMENTS for names of contributors to the
  code base.
*/
#ifdef DEBUG
#include <string.h>
#include "types.h"
#include "intbits.h"
#include "visible.h"
#include "debugdef.h"
#include "streedef.h"
#include "streeacc.h"
#include "protodef.h"
#include "streeproto.h"

#define SETLEAFUSED(V)   SETIBIT(leafused,V)
#define ISLEAFUSED(V)    ISIBITSET(leafused,V)

#define SETBRANCHUSED(V) SETIBIT(branchused,V)
#define ISBRANCHUSED(V)  ISIBITSET(branchused,V)

#define PROCESSALL(STOP,DP)\
        for(r = (Sint) loc.edgelen-1; r >= STOP; r--)\
        {\
          loc.remain = (Uint) r;\
          loc.locstring.length = DP + loc.edgelen - r;\
          DEBUGCODE(3,showlocation(stdout,stree,&loc));\
          DEBUG0(3,"\n");\
          processloc(stree,&loc);\
        }

static void showthesymbolstring(FILE *fp,SYMBOL *tlast,SYMBOL *left,
                                SYMBOL *right)
{
  SYMBOL *ptr;

  for(ptr=left; ptr<=right; ptr++)
  {
    if(ptr == tlast)
    {
      (void) putc('~',fp);
      return;
    } 
    if(ptr > left + 10)
    {
      fprintf(fp,"...");
      return;
    }
    SHOWCHARFP(fp,*ptr);
  }
}

static char *showsymbol(SYMBOL c)
{
  static char outbuf[100+1];

  if(INVISIBLE(c))
  {
    sprintf(outbuf,"\\%lu",(Showuint) c);
  } else
  {
    sprintf(outbuf,"%c",c);
  }
  return outbuf;
}
 
 Uint getlargelinkstree(Suffixtree *stree,Uint *btptr,Uint depth);

void showtable(Suffixtree *stree,BOOL final)
{
  Uint *largeptr, *btptr, *succptr, *rcptr, i,
       succdepth, distance, 
       nodeaddress, succ, depth, child, brother, 
       headposition, suffixlink;
  Uint leafindex, edgelen;
  SYMBOL *leftpointer;

  for(rcptr = stree->rootchildren; 
      rcptr <= stree->rootchildren + LARGESTCHARINDEX;
      rcptr++)
  {
    if(*rcptr != UNDEFINEDREFERENCE)
    {
      printf("rootchildren[%c]=",(char) (rcptr - stree->rootchildren));
      if(ISLEAF(*rcptr))
      {
        printf("Leaf %lu\n",(Showuint) GETLEAFINDEX(*rcptr));
      } else
      {
        succptr = stree->branchtab + GETBRANCHINDEX(*rcptr);
        printf("%s %lu\n",ISLARGE(*succptr) ? "Large" : "Small",
                         (Showuint) GETBRANCHINDEX(*rcptr));
      }
    }
  }
  if(final)
  {
    printf("rootchildren[~]=Leaf %lu\n",(Showuint) stree->textlen);
  }
  for(i=0; i<stree->nextfreeleafnum; i++)
  {
    printf("leaftab[%lu]=",(Showuint) i);
    SHOWINDEX(stree->leaftab[i]);
    printf("\n");
    (void) fflush(stdout);
  }
  printf(" Root:[");
  for(rcptr = stree->rootchildren; 
      rcptr <= stree->rootchildren + LARGESTCHARINDEX;
      rcptr++)
  {
    if(*rcptr != UNDEFINEDREFERENCE)
    {
      (void) putchar('(');
      if(ISLEAF(*rcptr))
      {
        leftpointer = stree->text + GETLEAFINDEX(*rcptr);
        showthesymbolstring(stdout,stree->sentinel,leftpointer,stree->sentinel);
        printf(",Leaf %lu)",(Showuint) GETLEAFINDEX(*rcptr));
      } else
      {
        succptr = stree->branchtab + GETBRANCHINDEX(*rcptr);
        GETBOTH(succdepth,headposition,succptr);
        leftpointer = stree->text + headposition;
        showthesymbolstring(stdout,stree->sentinel,leftpointer,leftpointer + succdepth - 1);
        printf(",%s %lu)",ISLARGE(*succptr) ? "Large" : "Small",
                          (Showuint) GETBRANCHINDEX(*rcptr));
      }
      (void) fflush(stdout);
    }
  }
  if(final)
  {
    printf(",(~,Leaf %lu)]\n",(Showuint) stree->textlen);
  } else
  {
    printf("]\n");
  }
  btptr = stree->branchtab + LARGEINTS; // skip the root
  printf("nodecount=%lu\n",(Showuint) stree->nodecount);
  for(i=UintConst(1); i < stree->nodecount; i++)
  {
    nodeaddress = BRADDR2NUM(stree,btptr);
    child = GETCHILD(btptr);
    brother = GETBROTHER(btptr);
    GETBOTH(depth,headposition,btptr);
    if(ISLARGE(*btptr))
    {
      printf(" L-Node %lu\"",(Showuint) nodeaddress);
      suffixlink = getlargelinkstree(stree,btptr,depth);
      btptr += LARGEINTS;
    } else
    {
      printf(" S-Node %lu\"",(Showuint) nodeaddress);
      suffixlink = nodeaddress + SMALLINTS;
      btptr += SMALLINTS;
    }
    showthesymbolstring(stdout,stree->sentinel,stree->text + headposition,
                                   stree->text + headposition + depth - 1);
    printf("\"(D=%lu,SN=%lu,SL=%lu,C=",(Showuint) depth,
                                       (Showuint) headposition,
                                       (Showuint) suffixlink);
    SHOWINDEX(child);
    printf(",B=");
    SHOWINDEX(brother);
    printf(")[");
    (void) fflush(stdout);
    succ = child;
    do 
    {
      (void) putchar('(');
      if(ISLEAF(succ))
      {
        leafindex = GETLEAFINDEX(succ);
        leftpointer = stree->text + depth + leafindex;
        showthesymbolstring(stdout,stree->sentinel,leftpointer,stree->sentinel);
        printf(",Leaf %lu)",(Showuint) leafindex);
        succ = LEAFBROTHERVAL(stree->leaftab[leafindex]);
      } else
      {
        succptr = stree->branchtab + GETBRANCHINDEX(succ);
        GETBOTH(succdepth,headposition,succptr);
        leftpointer = stree->text + depth + headposition;
        edgelen = succdepth - depth;
        showthesymbolstring(stdout,stree->sentinel,leftpointer,leftpointer + edgelen - 1);
        printf(",%s %lu)",ISLARGE(*succptr) ? "Large" : "Small",
                          (Showuint) GETBRANCHINDEX(succ));
        succ = GETBROTHER(succptr);
      }
    } while(!NILPTR(succ));
    printf("]\n");
    (void) fflush(stdout);
  }
}

/* 
   Check the following:
   (1) for each branching node there exist between 2 and 257 successors
   (2) for each branching node the list of successors is strictly ordered
       according to the first character of the edge label
   (3) there are no empty edge labels
   (4) there are \(n+1\) leaves and for each leaf there is exactly one 
       incoming edge
   (5) for each branching node (except for the root) there is exactly one
       incomming edge.
   (6) each suffix link point to a node whose depth is one smaller and
       whose headposition is either smaller or one larger than the previous node
   (7) the last branching node is not small
*/

void checkstree(Suffixtree *stree)
{
  Uint *largeptr, *btptr, *succptr, *slinkptr, lastsmall = 0, succdepth, 
       succheadposition, distance, succ, depth, leafindex, headposition, 
       edgelen, identitycount = 0, edgecount = 0,
       succcount, j, linkdepth, linkheadposition, *leafused, *branchused;
  Sint prevfirstchar, currentfirstchar;
#ifdef STREELARGE
  Uint largedepth = 0; 
#endif
  INITBITTAB(leafused,stree->textlen+1);
  INITBITTAB(branchused,stree->textlen+1);
  btptr = stree->branchtab; 
  while(btptr < stree->nextfreebranch)
  {
    succcount = 0;
    prevfirstchar = -1;
    GETBOTH(depth,headposition,btptr);
    succ = GETCHILD(btptr);
    do 
    {
      edgecount++;
      if(ISLEAF(succ))
      {
        leafindex = GETLEAFINDEX(succ);
        if(headposition == leafindex)
        {
          identitycount++;
        }
        if(ISLEAFUSED(leafindex))
        {
          fprintf(stderr,"Node %lu: more than one leaf edge to %lu\n",
                          (Showuint) BRADDR2NUM(stree,btptr),
                          (Showuint) leafindex);
          exit(EXIT_FAILURE);
        }
        SETLEAFUSED(leafindex);
        if(depth + leafindex < stree->textlen)
        {
          currentfirstchar = (Sint) stree->text[depth + leafindex];
        } else
        {
          currentfirstchar = LARGESTCHARINDEX + 1;
        }
        edgelen = stree->textlen + 1 - leafindex - depth;
        succ = LEAFBROTHERVAL(stree->leaftab[leafindex]);
      } else
      {
        succptr = stree->branchtab + GETBRANCHINDEX(succ);
        GETBOTH(succdepth,succheadposition,succptr);
        currentfirstchar = (Sint) stree->text[depth + succheadposition];
        edgelen = succdepth - depth;
        if(ISBRANCHUSED(succheadposition))
        {
          fprintf(stderr,"Node %lu: more than one edge to branch node %lu\n",
                          (Showuint) GETBRANCHINDEX(succ),
                          (Showuint) succheadposition);
          exit(EXIT_FAILURE);
        }
        SETBRANCHUSED(succheadposition);
        succ = GETBROTHER(succptr);
      }
      if(edgelen == 0)
      {
        fprintf(stderr,"Node %lu: outgoing '%s'-edge of length %lu\n",
                        (Showuint) BRADDR2NUM(stree,btptr),
                        showsymbol((SYMBOL) currentfirstchar),
                        (Showuint)edgelen);
        exit(EXIT_FAILURE);
      }
      if(prevfirstchar >= currentfirstchar)
      {
        fprintf(stderr,"Node %lu: '%s'-edge +",
                       (Showuint) BRADDR2NUM(stree,btptr),
                       showsymbol((SYMBOL) prevfirstchar));
        fprintf(stderr," '%s'-edge not in correct order\n",
                       showsymbol((SYMBOL) currentfirstchar));
        exit(EXIT_FAILURE);
      }
      prevfirstchar = currentfirstchar;
      succcount++;
    } while(!NILPTR(succ));
    if(succcount < (Uint) 2 || succcount > (Uint) (LARGESTCHARINDEX + 2))
    {
      fprintf(stderr,"Node %lu: %lu successors\n",
                (Showuint) BRADDR2NUM(stree,btptr),
                (Showuint) succcount);
      exit(EXIT_FAILURE);
    }
    NEXTNODE(btptr);
  }
  for(j=0; j<=stree->textlen; j++)
  {
    if(!ISLEAFUSED(j))
    {
      fprintf(stderr,"no leaf edge to %lu\n",(Showuint) j);
      exit(EXIT_FAILURE);
    }
  }
  btptr = stree->branchtab + LARGEINTS; 
  while(btptr < stree->nextfreebranch)
  {
    GETBOTH(depth,headposition,btptr);
    if(!ISBRANCHUSED(headposition))
    {
      fprintf(stderr,"no edge to branch node with headposition %lu\n",
                 (Showuint) headposition);
      exit(EXIT_FAILURE);
    }
    if(ISLARGE(*btptr))
    {
#ifdef STREELARGE
      if(!ISSMALLDEPTH(depth))
      {
        largedepth++;
      }
#endif
      lastsmall = 0;
      slinkptr = stree->branchtab + getlargelinkstree(stree,btptr,depth);
    } else
    {
      lastsmall = BRADDR2NUM(stree,btptr);
      slinkptr = btptr + SMALLINTS;
    }
    GETBOTH(linkdepth,linkheadposition,slinkptr);
    if(linkdepth + 1 != depth)
    {
      fprintf(stderr,"Node %lu(depth %lu) is linked to node (%lu,depth %lu)\n",
                      (Showuint) BRADDR2NUM(stree,btptr),
                      (Showuint) depth,
                      (Showuint) BRADDR2NUM(stree,slinkptr),
                      (Showuint) linkdepth);
      exit(EXIT_FAILURE);
    }
    if(headposition + 1 < linkheadposition)
    {
      fprintf(stderr,"Node %lu(headposition %lu) is linked to node"
                     " %lu(headposition %lu)\n",
                      (Showuint) BRADDR2NUM(stree,btptr),
                      (Showuint) headposition,
                      (Showuint) BRADDR2NUM(stree,slinkptr),
                      (Showuint) linkheadposition);
      exit(EXIT_FAILURE);
    }
    NEXTNODE(btptr);
  }
  if(lastsmall != 0)
  {
    fprintf(stderr,"Node %lu is the last node but it is small\n",
                   (Showuint) lastsmall);
    exit(EXIT_FAILURE);
  }
/*
  if(edgecount - identitycount -1  > (stree->textlen * 3)/2)
  {
    fprintf(stderr,"htabsize %lu >= %lu-%lu too large\n",
                  (stree->textlen * 3)/2,
                  (Showuint) edgecount,(Showuint) identitycount);
    exit(EXIT_FAILURE);
  }
*/
  FREESPACE(leafused);
  FREESPACE(branchused);
  DEBUG2(2,"#edgecount %lu identitycount %lu\n",
              (Showuint) edgecount,
              (Showuint) identitycount);
#ifdef STREELARGE
  DEBUG1(2,"#largedepth %lu\n",(Showuint) largedepth);
#endif
}

static void showsubtree(Suffixtree *stree,Uint indent,Uint *btptr)
{
  Uint *largeptr, *succptr, leafindex, succdepth, edgelen, succ, distance, 
       depth, headposition; 
  SYMBOL *leftpointer;

  GETBOTH(depth,headposition,btptr);
  succ = GETCHILD(btptr);
  do 
  {
    printf("%*.*s",(Fieldwidthtype) indent,(Fieldwidthtype) indent,"");
#ifdef SHOWLEAD
    SHOWINDEX(succ);
#endif 
    if(ISLEAF(succ))
    {
      leafindex = GETLEAFINDEX(succ);
      leftpointer = stree->text + depth + leafindex;
      showthesymbolstring(stdout,stree->sentinel,leftpointer,stree->sentinel);
      (void) putchar('\n');
      succ = LEAFBROTHERVAL(stree->leaftab[leafindex]);
    } else
    {
      succptr = stree->branchtab + GETBRANCHINDEX(succ);
      GETBOTH(succdepth,headposition,succptr);
      leftpointer = stree->text + depth + headposition;
      edgelen = succdepth - depth;
      showthesymbolstring(stdout,stree->sentinel,leftpointer,leftpointer + edgelen - 1);
      (void) putchar('\n');
      showsubtree(stree,indent+6,succptr);
      succ = GETBROTHER(succptr);
    } 
  } while(!NILPTR(succ));
}

void showstree(Suffixtree *stree)
{
  Uint *btptr, *rcptr, *largeptr, distance, headposition, succdepth;
  SYMBOL *leftpointer;

  for(rcptr = stree->rootchildren; 
      rcptr <= stree->rootchildren + LARGESTCHARINDEX;
      rcptr++)
  {
    if(*rcptr != UNDEFINEDREFERENCE)
    {
#ifdef SHOWLEAD
      SHOWINDEX(*rcptr);
#endif
      if(ISLEAF(*rcptr))
      {
        leftpointer = stree->text + GETLEAFINDEX(*rcptr);
        showthesymbolstring(stdout,stree->sentinel,leftpointer,stree->sentinel);
        (void) putchar('\n');
      } else
      {
        btptr = stree->branchtab + GETBRANCHINDEX(*rcptr);
        GETBOTH(succdepth,headposition,btptr);
        leftpointer = stree->text + headposition;
        showthesymbolstring(stdout,stree->sentinel,leftpointer,leftpointer + succdepth - 1);
        (void) putchar('\n');
        showsubtree(stree,UintConst(6),btptr);
      }
    }
  }
  printf("~\n");
}

void showstate(Suffixtree *stree)
{
  if(stree->headend == NULL)
  {
    printf(" head=%lu of depth %lu\n",
                (Showuint) BRADDR2NUM(stree,stree->headnode),
                (Showuint) stree->headnodedepth);
  } else
  {
    printf(" head=[%lu,",(Showuint) BRADDR2NUM(stree,stree->headnode));
    showthesymbolstring(stdout,stree->sentinel,stree->headstart,stree->headend);
    printf("] where headnode is of depth %lu\n",
                         (Showuint) stree->headnodedepth);
  }
  printf(" insertnode=");
  SHOWINDEX(stree->insertnode);
  printf("\n insertprev=");
  SHOWINDEX(stree->insertprev);
  printf("\n smallnode=%lu\n",(Showuint) stree->smallnode);
  printf(" largenode=%lu\n",(Showuint) stree->largenode);
  printf(" nodecount=%lu\n",(Showuint) stree->nodecount);
  if(stree->chainstart == NULL)
  {
    printf(" chainstart=UNDEFINEDREFERENCE\n");
  } else
  {
    printf(" chainstart=%lu\n",(Showuint) BRADDR2NUM(stree,stree->chainstart));
  }
  printf(" smallnotcompleted=%lu\n",(Showuint) stree->smallnotcompleted);
  printf(" nextfreebranch=%lu\n",
            (Showuint) BRADDR2NUM(stree,stree->nextfreebranch));
  printf(" nextfreeleafnum=%lu\n",(Showuint) stree->nextfreeleafnum);
  printf(" tail of length (%lu)=",(Showuint) (stree->sentinel-stree->tailptr));
  (void) fflush(stdout);
  showthesymbolstring(stdout,stree->sentinel,stree->tailptr,stree->sentinel);
  (void) putchar('\n');
  (void) fflush(stdout);
}

static void loc2stringstree(Suffixtree *stree,Stringtype *s,Location *loc)
{
  Branchinfo branchinfo;

  if(loc->nextnode.toleaf)
  {
    s->start = LEAFADDR2NUM(stree,loc->nextnode.address);
    s->length = stree->textlen - s->start - loc->remain + 1;
  } else
  {
    getbranchinfostree(stree,ACCESSDEPTH | ACCESSHEADPOS,&branchinfo,
                             loc->nextnode.address);
    s->start = branchinfo.headposition;
    s->length = branchinfo.depth - loc->remain;
  }
  if(s->length != loc->locstring.length)
  {
    fprintf(stderr,"s->length=%lu != %lu != loc->locstring.length\n",
                    (Showuint) s->length,
                    (Showuint) loc->locstring.length);
    exit(EXIT_FAILURE);
  }
  if(s->start != loc->locstring.start)
  {
    if(memcmp(stree->text+s->start,
              stree->text+loc->locstring.start,
              (size_t) s->length) != 0)
    {
      fprintf(stderr,"compare of strings failed: \"");
      (void) fwrite(stree->text+s->start,sizeof(SYMBOL),
                    (size_t) s->length,stderr);
      fprintf(stderr,"\" != \"");
      (void) fwrite(stree->text + loc->locstring.start,sizeof(SYMBOL),
                    (size_t) loc->locstring.length,stderr);
      fprintf(stderr,"\"\n");
      exit(EXIT_FAILURE);
    }
  }
}

void showlocation(FILE *fp,Suffixtree *stree,Location *loc)
{
  Stringtype lstr;

  fprintf(fp,"\"");
  loc2stringstree(stree,&lstr,loc);

  showthesymbolstring(fp,stree->sentinel,stree->text+lstr.start,
                             stree->text+lstr.start+lstr.length-1);
  fprintf(fp,"\"=(%lu,%lu,",(Showuint) loc->locstring.start,
                            (Showuint) loc->locstring.length);
  if(loc->remain > 0)
  {
    fprintf(fp,"Branch %lu,",(Showuint) BRADDR2NUM(stree,loc->previousnode));
    showthesymbolstring(fp,stree->sentinel,loc->firstptr,
                               loc->firstptr+loc->edgelen-loc->remain-1);
    fprintf(fp,",");
    if(loc->nextnode.toleaf)
    {
      showthesymbolstring(fp,stree->sentinel,loc->firstptr+loc->edgelen-loc->remain,
                                 stree->sentinel);
    } else
    {
      showthesymbolstring(fp,stree->sentinel,loc->firstptr+loc->edgelen-loc->remain,
                                 loc->firstptr+loc->edgelen-1);
    }
    fprintf(fp,",");
  } 
  if(loc->nextnode.toleaf)
  {
    fprintf(fp,"Leaf %lu",
            (Showuint) LEAFADDR2NUM(stree,loc->nextnode.address));
  } else
  {
    fprintf(fp,"Branch %lu",(Showuint) BRADDR2NUM(stree,loc->nextnode.address));
  }
  fprintf(fp,")");
}

#define SHOWREF(R)\
        ((Showuint) ((R).toleaf ? LEAFADDR2NUM(stree,(R).address)\
                                : BRADDR2NUM(stree,(R).address)))

static Sint comparelocs(Suffixtree *stree,Location *loc1,Location *loc2)
{
  if(loc1->remain != loc2->remain)
  {
    fprintf(stderr,"loc1->remain = %lu != %lu = loc2->remain\n",
              (Showuint) loc1->remain,
              (Showuint) loc2->remain);
    return -1;
  }
  if(loc1->nextnode.toleaf != loc2->nextnode.toleaf)
  {
    fprintf(stderr,"loc1->nextnode.toleaf = %s != %s = loc2->nextnode.toleaf\n",
                    SHOWBOOL(loc1->nextnode.toleaf),
                    SHOWBOOL(loc2->nextnode.toleaf));
    return -1;
  }
  if(loc1->nextnode.address != loc2->nextnode.address)
  {
    fprintf(stderr,"loc1->nextnode.address = %lu !="
                   " %lu = loc2->nextnode.address\n",
                    SHOWREF(loc1->nextnode),
                    SHOWREF(loc2->nextnode));
    return -1;
  }
  if(loc1->remain > 0)
  {
    if(loc1->firstptr != loc2->firstptr)
    {
      fprintf(stderr,"loc1->firstptr = %lu != %lu = loc2->firstptr\n",
                      (Showuint) (loc1->firstptr - stree->text),
                      (Showuint) (loc2->firstptr - stree->text));
      return -1;
    }
    if(loc1->edgelen != loc2->edgelen)
    {
      fprintf(stderr,"loc1->edgelen = %lu != %lu = loc2->edgelen\n",
                      (Showuint) loc1->edgelen,
                      (Showuint) loc2->edgelen);
      return -1;
    }
    if(loc1->previousnode != loc2->previousnode)
    {
      fprintf(stderr,"loc1->previousnode = %lu != %lu = loc2->previousnode\n",
                      (Showuint) BRADDR2NUM(stree,loc1->previousnode),
                      (Showuint) BRADDR2NUM(stree,loc2->previousnode));
      return -1;
    }
  }
  return 0;
}

void checklocation(Suffixtree *stree,Location *loc)
{
  Uint rescanlength;
  Stringtype lstr, llstr;
  SYMBOL *rest;
  Location scanprefixloc, rescanloc, linklocloc;

  loc2stringstree(stree,&lstr,loc);

  for(rescanlength = 0; rescanlength <= lstr.length; rescanlength++)
  {
    rest = scanprefixfromnodestree(stree,&scanprefixloc,stree->branchtab,
                                   stree->text+lstr.start,
                                   stree->text+lstr.start+lstr.length-1,
                                   rescanlength);
    if(rest != NULL)
    {
      fprintf(stderr,"string of location ");
      showlocation(stderr,stree,loc);
      fprintf(stderr,"\nscanned = ");
      showlocation(stderr,stree,&scanprefixloc);
      fprintf(stderr,"\nnot found\n");
      exit(EXIT_FAILURE);
    } 
    if(comparelocs(stree,loc,&scanprefixloc) == -1)
    {
      fprintf(stderr,"compare loc and scanprefixloc for string \"");
      (void) fwrite(stree->text+lstr.start,sizeof(SYMBOL),
                    (size_t) lstr.length,stderr);
      fprintf(stderr,"\"\n");
      showlocation(stderr,stree,loc);
      fprintf(stderr,"!=");
      showlocation(stderr,stree,&scanprefixloc);
      fprintf(stderr,"\n");
      exit(EXIT_FAILURE);
    }
  }
  rescanstree(stree,&rescanloc,stree->branchtab,stree->text+lstr.start,
              stree->text+lstr.start+lstr.length-1);
  if(comparelocs(stree,loc,&rescanloc) == -1)
  {
    fprintf(stderr,"compare loc and rescanloc for string \"");
    (void) fwrite(stree->text+lstr.start,sizeof(SYMBOL),
                  (size_t) lstr.length,stderr);
    fprintf(stderr,"\"\n");
    showlocation(stderr,stree,loc);
    fprintf(stderr,"!=");
    showlocation(stderr,stree,&rescanloc);
    fprintf(stderr,"\n");
    exit(EXIT_FAILURE);
  }
  if(!ROOTLOCATION(loc))
  {
    linklocstree(stree,&linklocloc,loc);
    loc2stringstree(stree,&llstr,&linklocloc);
    if(llstr.length + 1 != lstr.length || 
       memcmp(stree->text+lstr.start+1,stree->text+llstr.start,
              (size_t) llstr.length) != 0)
    {
      fprintf(stderr,"linkloc(");
      (void) fwrite(stree->text+lstr.start,sizeof(SYMBOL),
                    (size_t) lstr.length,stderr);
      fprintf(stderr,")=");
      (void) fwrite(stree->text+llstr.start,sizeof(SYMBOL),
                    (size_t) llstr.length,stderr);
      showlocation(stderr,stree,&linklocloc);
      fprintf(stderr," is wrong\n");
      exit(EXIT_FAILURE);
    }
  }
}

static void enumlocationssubtree(Suffixtree *stree,Uint *btptr,
                                  void(*processloc)(Suffixtree *stree,
                                                    Location *))
{
  Location loc;
  Uint leafindex, succ, *largeptr, distance, depth, headposition; 
  Branchinfo branchinfo;
  Sint r;

  GETBOTH(depth,headposition,btptr);
  succ = GETCHILD(btptr);
  do 
  {
#ifdef SHOWLEAD
    SHOWINDEX(succ);
#endif 
    if(ISLEAF(succ))
    {
      leafindex = GETLEAFINDEX(succ);
      loc.nextnode.toleaf = True;
      loc.nextnode.address = stree->leaftab + leafindex;
      loc.firstptr = stree->text + depth + leafindex;
      loc.edgelen = stree->textlen - (depth + leafindex) + 1;
      loc.previousnode = btptr;
      loc.locstring.start = leafindex;
      PROCESSALL(1,depth);
      succ = LEAFBROTHERVAL(stree->leaftab[leafindex]);
    } else
    {
      loc.nextnode.toleaf = False;
      loc.nextnode.address = stree->branchtab + GETBRANCHINDEX(succ);
      getbranchinfostree(stree,ACCESSDEPTH | ACCESSHEADPOS,&branchinfo,
                               loc.nextnode.address);
      loc.firstptr = stree->text + depth + branchinfo.headposition;
      loc.edgelen = branchinfo.depth - depth;
      loc.previousnode = btptr;
      loc.locstring.start = branchinfo.headposition;
      PROCESSALL(0,depth);
      enumlocationssubtree(stree,loc.nextnode.address,processloc);
      succ = GETBROTHER(loc.nextnode.address);
    } 
  } while(!NILPTR(succ));
}

void enumlocations(Suffixtree *stree,
                   void(*processloc)(Suffixtree *stree,Location *))
{
  Location loc;
  Uint leafindex, *rcptr;
  Branchinfo branchinfo;
  Sint r;

  loc.remain = 0;
  loc.nextnode.toleaf = False;
  loc.nextnode.address = stree->branchtab;
  loc.locstring.start = loc.locstring.length = 0;
  DEBUGCODE(3,showlocation(stdout,stree,&loc));
  DEBUG0(3,"\n");
  processloc(stree,&loc);
  for(rcptr = stree->rootchildren; 
      rcptr <= stree->rootchildren + LARGESTCHARINDEX;
      rcptr++)
  {
    if(*rcptr != UNDEFINEDREFERENCE)
    {
#ifdef SHOWLEAD
      SHOWINDEX(*rcptr);
#endif
      if(ISLEAF(*rcptr))
      {
        leafindex = GETLEAFINDEX(*rcptr);
        loc.nextnode.toleaf = True;
        loc.nextnode.address = stree->leaftab + leafindex;
        loc.firstptr = stree->text + leafindex;
        loc.edgelen = stree->textlen - leafindex + 1;
        loc.previousnode = stree->branchtab;
        loc.locstring.start = leafindex;
        PROCESSALL(1,0);
      } else
      {
        loc.nextnode.toleaf = False;
        loc.nextnode.address = stree->branchtab + GETBRANCHINDEX(*rcptr);
        getbranchinfostree(stree,ACCESSDEPTH | ACCESSHEADPOS,&branchinfo,
                                 loc.nextnode.address);
        loc.firstptr = stree->text + branchinfo.headposition;
        loc.edgelen = branchinfo.depth;
        loc.previousnode = stree->branchtab;
        loc.locstring.start = branchinfo.headposition;
        PROCESSALL(0,0);
        enumlocationssubtree(stree,loc.nextnode.address,processloc);
      }
    }
  }
}
#endif
