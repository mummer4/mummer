/*
  Copyright (c) 2003 by Stefan Kurtz and The Institute for
  Genomic Research.  This is OSI Certified Open Source Software.
  Please see the file LICENSE for licensing information and
  the file ACKNOWLEDGEMENTS for names of contributors to the
  code base.
*/

#ifndef STREEMAC_H
#define STREEMAC_H

#include "types.h"
#include "symboldef.h"
#include "errordef.h"

/*
  For each branching node we store five integers, see \cite{KUR:1999}.
  An or-combination of the following bits allow to access any subset of
  these five integers via the function \texttt{getbranchinfostree}.
*/

#define ACCESSDEPTH          UintConst(1)
#define ACCESSHEADPOS        (UintConst(1) << 1)
#define ACCESSSUFFIXLINK     (UintConst(1) << 2)
#define ACCESSFIRSTCHILD     (UintConst(1) << 3)
#define ACCESSBRANCHBROTHER  (UintConst(1) << 4)

/*
  The following macro simplifies calling the function for constructing
  suffix trees.
*/

#define CONSTRUCTSTREE(ST,TEXT,TEXTLEN,ACTION)\
        if(constructstree(ST,TEXT,TEXTLEN) != 0)\
        {\
          fprintf(stderr,"%s",messagespace());\
          ACTION;\
        }

/*
  The root of the suffix tree is stored in the first elements of the 
  table \texttt{branchtab}.
*/

#define ROOT(ST)            ((ST)->branchtab)

/*
  The following macro returns \texttt{True}, if and only if the given 
  location refers to the root.
*/

#define ROOTLOCATION(LOC)\
        (((LOC)->locstring.length == 0) ? True : False)

/*
  The following macros compute the index of a branch node and a leaf,
  relative to the first address and the first leaf, respectively.
*/

#define BRADDR2NUM(ST,A)      ((Uint) ((A) - ROOT(ST)))
#define LEAFADDR2NUM(ST,A)    ((Uint) ((A) - (ST)->leaftab))

//\Ignore{

#ifdef DEBUG
#define CHECKADDR(ST,A)\
        if((A).toleaf)\
        {\
          if(LEAFADDR2NUM(ST,(A).address) > (ST)->textlen)\
          {\
            printf("%s,%lu:",__FILE__,(Showuint) __LINE__);\
            printf("leafaddr = %lu invalid\n",\
                    (Showuint) LEAFADDR2NUM(ST,(A).address));\
            exit(EXIT_FAILURE);\
          }\
        } else\
        {\
          if(BRADDR2NUM(ST,(A).address) >= (ST)->nextfreebranchnum)\
          {\
            printf("%s,%lu:",__FILE__,(Showuint) __LINE__);\
            printf("branchaddr = %lu invalid\n",\
                    (Showuint) BRADDR2NUM(ST,(A).address));\
            exit(EXIT_FAILURE);\
          }\
        } 
#else
#define CHECKADDR(ST,A) /* Nothing */
#endif

//}

#endif
