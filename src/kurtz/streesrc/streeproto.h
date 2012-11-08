/* 
  This file is generated. Do not edit.

  A Library for the Efficient Construction and Application of Suffix Trees

  Copyright (c) 2003 by Stefan Kurtz and The Institute for
  Genomic Research.  This is OSI Certified Open Source Software.
  Please see the file LICENSE for licensing information and
  the file ACKNOWLEDGEMENTS for names of contributors to the
  code base.
*/

#ifndef STREEPROTO_H
#define STREEPROTO_H

#ifdef __cplusplus
extern "C" {
#endif

Sint constructstree(Suffixtree *stree,SYMBOL *text,Uint textlen);
Sint constructmarkmaxstree(Suffixtree *stree,SYMBOL *text,Uint textlen);
Sint constructheadstree(Suffixtree *stree,SYMBOL *text,Uint textlen,void(*processhead)(Suffixtree *,Uint,void *),void *processheadinfo);
Sint constructprogressstree(Suffixtree *stree,SYMBOL *text,Uint textlen,void (*progress)(Uint,void *),void (*finalprogress)(void *),void *info);

void freestree(Suffixtree *stree);
#ifdef __cplusplus
}
#endif

#ifdef __cplusplus
extern "C" {
#endif
Uint getlargelinkstree(Suffixtree *stree,Bref btptr,Uint depth);
Uint getlargelinkstree(Suffixtree *stree,Bref btptr,Uint depth);
Uint getlargelinkstree(/*@unused@*/ Suffixtree *stree,Bref btptr,Uint depth);
void getleafinfostree(Suffixtree *stree,Leafinfo *leafinfo,Lref lptr);
void getbranchinfostree(Suffixtree *stree,Uint whichinfo,
                        Branchinfo *branchinfo,Bref btptr);
void getheadstringstree(Suffixtree *stree,Stringtype *str);
Uint getmaxtextlenstree(void);
void showpathstree(Suffixtree *stree,Bref bnode,
                   void (*showchar)(SYMBOL,void *),void *info);
void showsimplelocstree(Suffixtree *stree,Simpleloc *loc);
void showsimplelocliststree(Suffixtree *stree,ArraySimpleloc *ll);
void rootsucclocationsstree(Suffixtree *stree,ArraySimpleloc *ll);
void succlocationsstree(Suffixtree *stree,BOOL nosentinel,Simpleloc *loc,
                        ArraySimpleloc *ll);
/*@null@*/ SYMBOL *scanprefixfromnodestree(Suffixtree *stree,Location *loc,
                                           Bref btptr,SYMBOL *left,
                                           SYMBOL *right,Uint rescanlength);
/*@null@*/ SYMBOL *scanprefixstree(Suffixtree *stree,Location *outloc,
                                   Location *inloc,SYMBOL *left,
                                   SYMBOL *right,Uint rescanlength);
/*@null@*/SYMBOL *findprefixpathfromnodestree(Suffixtree *stree,
                                              ArrayPathinfo *path,
                                              Location *loc,
                                              Bref btptr,
                                              SYMBOL *left,
                                              SYMBOL *right,
                                              Uint rescanlength);
/*@null@*/ SYMBOL *findprefixpathstree(Suffixtree *stree,
                                       ArrayPathinfo *path,
                                       Location *outloc,
                                       Location *inloc,
                                       SYMBOL *left,
                                       SYMBOL *right,
                                       Uint rescanlength);
void rescanstree(Suffixtree *stree,Location *loc,
                 Bref btptr,SYMBOL *left,SYMBOL *right);
void linklocstree(Suffixtree *stree,Location *outloc,Location *inloc);
void showdepthtab(ArrayUint *dt);
void makedepthtabstree(ArrayUint *depthtab,Suffixtree *stree);
BOOL exactlytwoleavesstree(Suffixtree *stree,PairUint *twoleaves,Bref start);
Sint depthfirststree(Suffixtree *stree,Reference *startnode,
                     Sint (*processleaf)(Uint,Bref,void *),
                     BOOL (*processbranch1)(Bref,void *),
                     Sint (*processbranch2)(Bref,void *),
                     BOOL (*stoptraversal)(void *),void *stopinfo,void *info);
Sint makeleaflist(Suffixtree *stree,ArrayUint *leaflist,Reference *start);
void overallstree(Suffixtree *stree,BOOL skiproot,
                  void(*processnode)(Suffixtree *,Bref,Uint,Uint,void *),
                  void *info);
void overmaximalstree(Suffixtree *stree,
                      void(*processnode)(Suffixtree *,Bref,Uint,Uint,void *),
                      void *info);
void oversuccsstree(Suffixtree *stree,Bref bnode,
                    void(*processleaf)(Suffixtree *,Uint,void *),
                    void(*processbranch)(Suffixtree *,Bref,void *),
                    void *info);
Uint getleafcountstree(Suffixtree *stree,Bref nodeptr);
Sint addleafcountsstree(Suffixtree *stree);
/*@null@*/ Bref firstbranchingnode(Suffixtree *stree);
/*@null@*/ Bref nextbranchingnode(Suffixtree *stree,Bref bptr);
Lref firstleaf(Suffixtree *stree);
/*@null@*/ Lref nextleaf(Suffixtree *stree,Lref lptr);
/*@null@*/ Reference *firstnode(Suffixtree *stree,Reference *refspace);
/*@null@*/ Reference *nextnode(Suffixtree *stree,Reference *nref,
                               Reference *refspace);
/*@null@*/ Reference *firstsucc(Suffixtree *stree,Bref bptr,
                                Reference *refspace);
/*@null@*/ Reference *rightbrother(Suffixtree *stree,Reference *node);
Reference *firstnodedfs(Suffixtree *stree,DFSstate *dfsstate,
                        Reference *current);
/*@null@*/ Reference *nextnodedfs(Suffixtree *stree,Reference *current,
                                  DFSstate *dfsstate);
void showtable(Suffixtree *stree,BOOL final);
void checkstree(Suffixtree *stree);
void showstree(Suffixtree *stree);
void showstate(Suffixtree *stree);
void showlocation(FILE *fp,Suffixtree *stree,Location *loc);
void checklocation(Suffixtree *stree,Location *loc);
void enumlocations(Suffixtree *stree,
                   void(*processloc)(Suffixtree *stree,Location *));
#ifdef __cplusplus
}
#endif
#endif
