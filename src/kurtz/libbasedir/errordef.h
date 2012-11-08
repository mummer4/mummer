/*
  Copyright (c) 2003 by Stefan Kurtz and The Institute for
  Genomic Research.  This is OSI Certified Open Source Software.
  Please see the file LICENSE for licensing information and
  the file ACKNOWLEDGEMENTS for names of contributors to the
  code base.
*/

//\Ignore{

#ifndef ERRORDEF_H
#define ERRORDEF_H
#include <stdio.h>
#include <stdlib.h>
#include "types.h"

#ifdef __cplusplus
  extern "C" {
#endif

char *messagespace(void);
Sint maxerrormsg(void);

#ifdef __cplusplus
}
#endif

//}

/*
  This file contains some macros to write error messages into a
  buffer returned by the function \texttt{messagespace}. 
*/

/*
  There is a generic macro \texttt{GENERROR} (definition
  not given) that checks if the 
  result of the computation \texttt{C} exceed the value returned 
  by the function \texttt{maxerrormessage}. If so, then a 
  corresponding error message is written to stderr.
*/

//\IgnoreLatex{

#ifdef DEBUG
#define THROWERRORLINE\
        DEBUG2(1,"# throw error message in %s line %lu\n",__FILE__,\
                                                (Showuint) __LINE__)
#else
#define THROWERRORLINE /* Nothing */
#endif

#define GENERROR(C);\
        THROWERRORLINE;\
        if((C) >= maxerrormsg())\
        {\
          fprintf(stderr,"file %s, line %lu: "\
                         "space for errormessage too small\n",\
                  __FILE__,(Showuint) __LINE__);\
          exit(EXIT_FAILURE);\
        }

/*
  The different error macros call \texttt{sprintf} with the corresponding
  number of arguments.
*/


#define ERROR0(F)\
        GENERROR(sprintf(messagespace(),F))

#define ERROR1(F,A1)\
        GENERROR(sprintf(messagespace(),F,A1))

#define ERROR2(F,A1,A2)\
        GENERROR(sprintf(messagespace(),F,A1,A2))

#define ERROR3(F,A1,A2,A3)\
        GENERROR(sprintf(messagespace(),F,A1,A2,A3))

#define ERROR4(F,A1,A2,A3,A4)\
        GENERROR(sprintf(messagespace(),F,A1,A2,A3,A4))

#define ERROR5(F,A1,A2,A3,A4,A5)\
        GENERROR(sprintf(messagespace(),F,A1,A2,A3,A4,A5))

//}


/*
  The following is a macro to show the usage line for all programs
  which have options and an indexname as the last argument.
*/

#define USAGEOUT\
        ERROR2("Usage: %s options indexname\n"\
               "%s -help shows possible options",\
                argv[0],argv[0]);

/*
  The following is the standard message in the main function. It shows the
  program name and the error message as returned by the function
  \texttt{messagespace}.
*/

#define STANDARDMESSAGE\
        fprintf(stderr,"%s: %s\n",argv[0],messagespace());\
        return EXIT_FAILURE

#define SIMPLESTANDARDMESSAGE\
        fprintf(stderr,"%s\n",messagespace());\
        return EXIT_FAILURE

/*
  The following is a generell error message which leads to a termination
  of the program.
*/

#define NOTSUPPOSED\
        fprintf(stderr,"%s: line %lu: This case is not supposed to occur\n",\
                       __FILE__,(Showuint) __LINE__);\
        exit(EXIT_FAILURE)

/*
  The following macro checks a ptr. If it is \texttt{NULL}, then the 
  program terminates with an error.
*/

#ifdef DEBUG
#define NOTSUPPOSEDTOBENULL(PTR)\
        if((PTR) == NULL)\
        {\
          NOTSUPPOSED;\
        }
#else
#define NOTSUPPOSEDTOBENULL(PTR) /* Nothing */
#endif

//\Ignore{

#endif

//}
