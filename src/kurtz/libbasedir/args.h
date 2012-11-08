/*
  Copyright (c) 2003 by Stefan Kurtz and The Institute for
  Genomic Research.  This is OSI Certified Open Source Software.
  Please see the file LICENSE for licensing information and
  the file ACKNOWLEDGEMENTS for names of contributors to the
  code base.
*/

//\Ignore{

#ifndef ARGS_H
#define ARGS_H
#include <stdlib.h>
#include <stdio.h>

//}

/*
  This header file implements three macros for checking argument numbers
  and to parse integers and floats from strings.
*/

/*
  The following macro checks if the number of argument is exactly
  \texttt{N}. Otherwise, an error message is thrown.
*/

#define CHECKARGNUM(N,S)\
        if (argc != N)\
        {\
          fprintf(stderr,"Usage: %s %s\n",argv[0],S);\
          exit(EXIT_FAILURE);\
        }

/*
  The following scans an integer \texttt{readint} from a string.
*/

#define PARSEINTARG(S)\
        if(sscanf(S,"%ld",&readint) != 1 || readint < 0)\
        {\
          fprintf(stderr,"invalid argument \"%s\": " \
                         "non-negative number expected\n",S);\
          exit(EXIT_FAILURE);\
        }

/*
  The following scans a floating point value \texttt{readfloat} from a 
  string.
*/

#define PARSEFLOATARG(S)\
        if(sscanf(S,"%f",&readfloat) != 1)\
        {\
          fprintf(stderr,"invalid argument \"%s\":"\
                         " floating point number expected\n",S);\
          exit(EXIT_FAILURE);\
        }

//\IgnoreLatex{

#endif

//}
