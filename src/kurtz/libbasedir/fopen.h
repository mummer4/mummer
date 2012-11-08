/*
  Copyright (c) 2003 by Stefan Kurtz and The Institute for
  Genomic Research.  This is OSI Certified Open Source Software.
  Please see the file LICENSE for licensing information and
  the file ACKNOWLEDGEMENTS for names of contributors to the
  code base.
*/

//\Ignore{

#ifndef FOPEN_H
#define FOPEN_H
#include <stdio.h>
#include <stdlib.h>

//}

/*
  This file defines macros for opening and writing files via
  file pointers.
*/

#define FILEOPEN(FP,FILENAME,MODE)\
        if ((FP = fopen(FILENAME,MODE)) == NULL)\
        {\
          fprintf(stderr,"(%s,%lu): Cannot open file \"%s\"\n",\
                  __FILE__,(Showuint) __LINE__,FILENAME);\
          exit(EXIT_FAILURE);\
        }

#define FPBINWRITE(FP,BUF,SIZE)\
        if(fwrite(BUF,SIZE,1,FP) != 1)\
        {\
          fprintf(stderr,"(%s,%lu): fwrite failed\n",__FILE__,\
                  (Showuint) __LINE__);\
          exit(EXIT_FAILURE);\
        }

//\Ignore{

#endif

//}
