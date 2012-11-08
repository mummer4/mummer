/*
  Copyright (c) 2003 by Stefan Kurtz and The Institute for
  Genomic Research.  This is OSI Certified Open Source Software.
  Please see the file LICENSE for licensing information and
  the file ACKNOWLEDGEMENTS for names of contributors to the
  code base.
*/

//\Ignore{

#ifndef SPACEDEF_H
#define SPACEDEF_H
#include "types.h"
#include "errordef.h"

#ifdef __cplusplus
  extern "C" {
#endif

/*@notnull@*/ void *allocandusespaceviaptr(char *file,Uint line,
                                           /*@null@*/ void *ptr,
                                           Uint size,Uint number);
/*@notnull@*/ char *dynamicstrdup(char *file,Uint line,char *source);
void freespaceviaptr(char *file,Uint line,void *ptr);
/*@null@*/ void *creatememorymapforfiledesc(char *file,Uint line,Sint fd,
                                            BOOL writemap,Uint numofbytes);
/*@null@*/ void *creatememorymap(char *file,Uint line,char *filename,
                                 BOOL writemap,Uint *numofbytes);
Sint deletememorymap(char *file,Uint line,void *mappedfile);

#ifdef __cplusplus
}
#endif

//}

/*
  This file defines macros to simplify the calls to the
  functions 
  \begin{itemize}
  \item
  \texttt{allocandusespaceviaptr}, 
  \item
  \texttt{freespaceviaptr},
  \item
  \texttt{dynamicstrdup}, 
  \item
  \texttt{creatememorymapforfiledesc}, 
  \item
  \texttt{creatememorymap}, 
  \item
  \texttt{delete\-memorymap}.
  \end{itemize}

  \begin{enumerate}
  \item
    The first parameter to \texttt{ALLOCSPACE} is \texttt{NULL}
    or a pointer to a previously
    allocated space block. 
  \item
  The second argument of the macro is the type of the space block
  to be allocated.
  \item
  The third argument is the number of elements of that type to be 
  allocated space for.
  \end{enumerate}
*/

#define ALLOCSPACE(S,T,N)\
        (T *) allocandusespaceviaptr(__FILE__,(Uint) __LINE__,\
                                     S,(Uint) sizeof(T),N)

/*
  The macro \texttt{FREESPACE} frees the space pointed to by \texttt{P},
  if this is not \texttt{NULL}. It also sets the 
  pointer to \texttt{NULL}.
*/

#define FREESPACE(P)\
        if((P) != NULL)\
        {\
          freespaceviaptr(__FILE__,(Uint) __LINE__,P);\
          P = NULL;\
        }

/*
  The remaining macros call the corresponding function with
  the filename and the line number where the function call 
  appears.
*/

#define DYNAMICSTRDUP(S)\
        dynamicstrdup(__FILE__,(Uint) __LINE__,S)

#define CREATEMEMORYMAP(F,WM,NB)\
        creatememorymap(__FILE__,(Uint) __LINE__,F,WM,NB)

#define CREATEMEMORYMAPFORFILEDESC(FD,WM,NB)\
        creatememorymapforfiledesc(__FILE__,(Uint) __LINE__,FD,WM,NB)

#define DELETEMEMORYMAP(MF)\
        deletememorymap(__FILE__,(Uint) __LINE__,(void *) MF)

//\Ignore{

#endif

//}
