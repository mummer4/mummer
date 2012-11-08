/*
  Copyright (c) 2003 by Stefan Kurtz and The Institute for
  Genomic Research.  This is OSI Certified Open Source Software.
  Please see the file LICENSE for licensing information and
  the file ACKNOWLEDGEMENTS for names of contributors to the
  code base.
*/

Sint constructstree(Suffixtree *stree,SYMBOL *text,Uint textlen)
{

}

/*EE
  In some applications it is convenient to perform some extra computations
  during the suffix tree construction.  For this purpose, there are 
  two variations of \texttt{constructstree}:

  The following function 
  additionally marks all branching nodes which are maximal.
  A branching node \(\overline{v}\) is \emph{maximal} if 
  and only if there is 
  no suffix link to \(\overline{v}\). The maximal nodes can be 
  enumerated using the function \texttt{overmaximalstree}, see Section 
  \ref{Traverse}.
*/

Sint constructmarkmaxstree(Suffixtree *stree,SYMBOL *text,
                           Uint textlen)
{

}

/*EE
  The following variation
  constructs the suffix tree and calls for each \(j\in[0,n]\) the 
  function \texttt{processhead} with the second argument being \(j\),
  as soon as \(\HD{j}\) is computed.  
  The first argument of \texttt{processhead} is the current suffix tree,
  and the third argument is the pointer \texttt{processheadinfo}. 
*/

Sint constructheadstree(Suffixtree *stree,SYMBOL *text,
                        Uint textlen,
                        void(*processhead)(Suffixtree *,Uint,
                                           void *),
                        void *processheadinfo)
{
}

/*EE
  The following variation
  constructs the suffix tree for \texttt{text}. Additionally, whenever 
  \(textlen>\Numofcalls\), the function
  \texttt{progress} is called after \(i\cdot textlen/\Numofcalls\) 
  construction
  steps, for \(i\in[1,\Numofcalls]\). Moreover, the function 
  \texttt{finalprogress} is called at the end of the suffix tree 
  construction. When they are called, both functions are supplied
  with the argument \texttt{info}.
*/

Sint constructprogressstree(Suffixtree *stree,SYMBOL *text,
                            Uint textlen,
                            void (*progress)(Uint,void *),
                            void (*finalprogress)(void *),
                            void *info)
{

}

/*EE
  For example, one could define \texttt{progress} and \texttt{finalprogress}
  as follows, in order to write a line of 79 dots on the given file pointer,
  to show the progress of the suffix tree construction:
*/

void progresswithdot(Uint nextstep,void *info)
{
  fputc('.',stderr);
  fflush(stderr);
}

void finalprogress(void *info)
{
  fputc('\n',stderr);
}

/*EE
  The following function stores for each branching node 
  \(\overline{v}\) the \emph{lefcount},
  i.e.\ the number of leafs in the subtree below \(\overline{v}\). 
  The computation of the leafcounts is done in linear time. It can be 
  retrieved using the function \texttt{getleafcountstree}, described 
  below.
*/

void addleafcountsstree(Suffixtree *stree)
{

}

/*EE
  The following function frees the space allocated for \texttt{stree}.
*/

void freestree(Suffixtree *stree)
{

}

/*EE
  The following function
  returns the maximal length of an input string allowed for the
  suffix tree construction.
*/

Uint getmaxtextlenstree(void)
{

}

/*EE
 \subsection{Accessing the Stored Information}

  The following function stores the information for leaf 
  \texttt{lref} in the structure \texttt{leafinfo}.
*/

void getleafinfostree(Suffixtree *stree,Leafinfo *leafinfo,
                      Lref lptr)
{

}

/*EE
  The following function
  stores the information for branching node \texttt{bnode} in
  \texttt{branchinfo}. For efficiency reason, the function does not always 
  compute all 5 components stored for a branching node. Instead, the 
  components are delivered as specified by \texttt{whichinfo}. This 
  has the following possible bits: 
  \begin{center}
  \begin{tabular}{l}
  \texttt{ACCESSDEPTH}\\
  \texttt{ACCESSHEADPOS}\\
  \texttt{ACCESSSUFFIXLINK}\\
  \texttt{ACCESSFIRSTCHILD}\\
  \texttt{ACCESSBRANCHBROTHER}
  \end{tabular}
  \end{center}
*/

void getbranchinfostree(Suffixtree *stree,Uint whichinfo,
                        Branchinfo *branchinfo,Bref btptr)
{

}

/*EE
  The following function
  shows the path for the node \texttt{bnode}. A symbol \(a\) is shown 
  by applying the function \texttt{showchar} to it and the 
  additional argument \texttt{info}.
*/

void showpathstree(Suffixtree *stree,Bref bnode,
                   void (*showchar)(SYMBOL,void *),void *info)
{

}

/*EE
  The following function stores a start position and the length of 
  \(\HD{j}\), for the current \(j\) in \texttt{str}. The start 
  position is smaller than \(j\). 
*/

void getheadstringstree(Suffixtree *stree,Stringtype *str)
{

}

/*EE
  The following function returns \texttt{True}, if and only if the 
  branching node \texttt{bnode} has exactly two leaves, say with leaf 
  numbers \(i\) and \(j\), \(i<j\). In this case, \(i\) is stored in 
  \texttt{twoleaves->uint0} and \(j\) is stored in 
  \texttt{twoleaves->uint1}.
*/

BOOL exactlytwoleavesstree(Suffixtree *stree,
                           PairUint *twoleaves,Bref start)
{

}

/*EE
  The following function
  returns the leafcount for a node referenced by \texttt{nodeptr}.
  The function \texttt{addleafcountstree} must have been called
  exactly once before, in order to  compute the leafcounts.
*/

Uint getleafcountstree(Suffixtree *stree,Bref nodeptr)
{

}

/*EE
  The following function
  computes a table \texttt{depthtab}, such that 
  \texttt{depthtab->spaceUint[i]}
  holds the number of branching nodes in the given 
  suffix tree whose depth is
  exactly \(i\). The maximal depth of any branching node is
  \texttt{depthtab->nextfreeUint-1}.
*/

void makedepthtabstree(ArrayUint *depthtab,Suffixtree *stree)
{

}

/*EE
  \subsection{Traversing the Suffix Tree}\label{Traverse}

  The following function
  implements the function \emph{scanprefix}, restricted to the case that
  the starting location is a node location. \texttt{stree} is the suffix
  tree, \texttt{outloc} is the resulting location, \texttt{startnode} is
  the branching node, the scanning starts from, and \texttt{left}
  and \texttt{right} delimit the string, say \(s\), to be scanned. 
  The function returns \texttt{NULL}, if \(s\) is scanned completely.
  Otherwise, it points to a suffix of \(s\). 
  \texttt{rescanlength} is the length of a prefix of 
  \(s\) that already occurs in the suffix tree. It is always safe
  to set \texttt{rescanlength} to 0.
*/

SYMBOL *scanprefixfromnodestree(Suffixtree *stree,Location *loc,
                                Bref btptr,SYMBOL *left,
                                SYMBOL *right,Uint rescanlength)
{

}

/*EE
  The following function 
  implements the function \emph{scanprefix}.  \texttt{stree} is the suffix
  tree, \texttt{outloc} is the resulting location, \texttt{inloc} is the 
  location, the scanning starts from, and \texttt{left}
  and \texttt{right} delimit the string, say \(s\), to be scanned.
  The function returns \texttt{NULL}, if \(s\) is scanned completely. 
  Otherwise, it points to a suffix of \(s\).
  \texttt{rescanlength} is the length of a prefix of 
  \(s\) that already occurs in the suffix tree. It is always safe
  to set \texttt{rescanlength} to 0.
*/

SYMBOL *scanprefixstree(Suffixtree *stree,Location *outloc,
                        Location *inloc,SYMBOL *left,
                        SYMBOL *right,Uint rescanlength)
{

}

/*EE
  The following function
  is similar to the function \emph{scanprefixfromnodestree}. It additionally
  stores in the array \texttt{path} the sequence of branching nodes
  that has been visited during the scan to location \texttt{loc}. 
  The array \texttt{path} may not be empty, in which case the list of 
  visited branching nodes is appended to it.
*/

SYMBOL *findprefixpathfromnodestree(Suffixtree *stree,
                                    ArrayPathinfo *path,
                                    Location *loc,Bref btptr,
                                    SYMBOL *left,SYMBOL *right,
                                    Uint rescanlength)
{

}

/*EE
  The following function
  is similar to the function \emph{scanprefixstree}. It additionally
  stores in the array \texttt{path} the sequence of branching nodes
  that has been visited during the scan to location \texttt{inloc}. 
  The array \texttt{path} may not be empty, in which case the list 
  of visited branching nodes is appended to it.
*/

SYMBOL *findprefixpathstree(Suffixtree *stree,
                            ArrayPathinfo *path,
                            Location *outloc,
                            Location *inloc,
                            SYMBOL *left,SYMBOL *right,
                            Uint rescanlength);

/*EE
  The following function
  implements the function \emph{rescan}. \texttt{stree} is the suffix
  tree, \texttt{outloc} is the resulting location, \texttt{startnode} 
  is the branching node, the scanning starts from, and \texttt{left}
  and \texttt{right} delimit the string, say \(s\), to be rescanned.
*/

void rescanstree(Suffixtree *stree,Location *loc,
                 Bref btptr,SYMBOL *left,SYMBOL *right)
{

}

/*EE
  The following function 
  implements the function \emph{linkloc}. \texttt{stree} is the 
  suffix tree, \texttt{outloc} is the resulting location, and 
  \texttt{inloc} is the input location.
*/

void linklocstree(Suffixtree *stree,Location *outloc,
                  Location *inloc)
{

}

/*EE
  The following function
  applies the function \texttt{processnode} to all branching nodes 
  of the suffix tree \texttt{stree}. The \emph{root} is skipped, if 
  and only if \texttt{skiproot} is \texttt{True}. The arguments of 
  \texttt{processnode} are as follows: The suffix tree \texttt{stree}, 
  the branching node, its depth, its head position, and the 
  pointer \texttt{info}.
*/

void overallstree(Suffixtree *stree,BOOL skiproot,
                  void(*processnode)(Suffixtree *,Bref,
                                     Uint,Uint,void *),
                  void *info)
{

}

/*EE
  The following function
  applies the function \texttt{processnode} to all maximal branching 
  nodes of the suffix tree \texttt{stree}. The arguments of 
  \texttt{processnode} are as follows: The suffix tree \texttt{stree}, 
  the reference to the branching node, its depth, its head position, 
  and the pointer \texttt{info}.
*/

void overmaximalstree(Suffixtree *stree,
                      void(*processnode)(Suffixtree *,Bref,
                                         Uint,Uint,void *),
                      void *info)
{

}

/*EE
  The following function
  enumerates all immediate successors of the branching node \texttt{bnode}. 
  Suppose the leaf with leaf number \texttt{j} is a 
  successor of \texttt{bnode}. Then the function \texttt{processleaf} with 
  arguments \texttt{stree}, \texttt{j}, and \texttt{info} is called.
  Suppose the branching node \texttt{bsucc} is a successor of 
  \texttt{bnode}. Then the function \texttt{processbranch} is called with
  arguments \texttt{stree}, \texttt{bsucc}, and \texttt{info}.
*/

void oversuccsstree(Suffixtree *stree,Bref bnode,
                    void(*processleaf)(Suffixtree *,Uint,void *),
                    void(*processbranch)(Suffixtree *,Bref,void *),
                    void *info)
{
}

/*EE
  The following function
  performs a (possibly) limited depth first traversal of the suffix 
  tree rooted by \texttt{startnode}. 
  \texttt{startnode} can be a reference to 
  a leaf or a reference to a branching node. The nodes of the subtree are 
  enumerated in depth first left to right order. The argument 
  \texttt{processleaf}
  is not allowed to be \texttt{NULL}. \texttt{processbranch1} and 
  \texttt{processbranch2} can either be both \texttt{NULL} or 
  both different from \texttt{NULL}. 
  \begin{itemize}
  \item
  Each time a leaf, say \(\overline{v}\), with leaf number
  \texttt{j} is encountered, the function \texttt{processleaf} is called
  with arguments \texttt{j}, \texttt{lca}, and 
  \texttt{info}. \texttt{lca} is the longest common ancestor of
  \(\overline{v}\) and the previous leaf encountered during the traversal. 
  \texttt{lca} is \texttt{NULL}, if \(\overline{v}\) is the first leaf 
  encountered in the traversal. If \texttt{processleaf} returns a value 
  smaller than 0, then \texttt{depthfirststree} terminates with a 
  return value \texttt{-1}.
  \item
  Each time a branching node \texttt{bsucc} is visited for the 
  first time, the function \texttt{process\-branch1} is called with 
  arguments \texttt{bsucc} and \texttt{info}. If 
  \texttt{processbranch1} returns \texttt{False}, then 
  the entire subtree below \texttt{bsucc} is discarded. Otherwise the
  depth first traversal continues with the subtree rooted by \texttt{bsucc}.
  \item
  Each time a branching node \texttt{bsucc} 
  is visited for the second time (i.e.\
  the entire subtree below \texttt{bsucc} has been processed), the function 
  \texttt{processbranch2} is called with arguments \texttt{bsucc} and
  \texttt{info}. If \texttt{processbranch2} returns a value smaller than 
  0, then \texttt{depthfirststree} returns with value \texttt{-1}.
  \end{itemize}
  Either \texttt{processleaf} or \texttt{processbranch1} and 
  \texttt{processbranch2} can be \texttt{NULL}, in which case there is no 
  function call. If \texttt{stoptraversal} is not \texttt{NULL}, then
  after each call the function \texttt{stoptraversal} is applied
  to \texttt{stopinfo}. If the return value of this call is \texttt{True},
  then the depth first traversal stops.

  In case everything goes right, \texttt{depthfirststree} returns 0.
*/

Sint depthfirststree(Suffixtree *stree,Reference *startnode,
                     Sint (*processleaf)(Uint,Bref,void *),
                     BOOL (*processbranch1)(Bref,void *),
                     Sint (*processbranch2)(Bref,void *),
                     BOOL (*stoptraversal)(void *),
                     void *stopinfo,void *info)
{

}

