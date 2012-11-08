/*
  Copyright (c) 2003 by Stefan Kurtz and The Institute for
  Genomic Research.  This is OSI Certified Open Source Software.
  Please see the file LICENSE for licensing information and
  the file ACKNOWLEDGEMENTS for names of contributors to the
  code base.
*/

//\IgnoreLatex{


#ifndef VISIBLE_H
#define VISIBLE_H
#include "types.h"

//}

/*
  This header file defines some constants and macros
  to test for visibility of a character (ASCII code is between
  33 and 126) and to show this character.
*/

/*
  The smallest visible character is the blank with code 33.
*/

#define LOWESTVISIBLE        33

/*
  The largest visible character is the tilde with code 126.
*/

#define HIGHESTVISIBLE      126

/*
  Check if character is invisible according to the definition from above.
*/

#define INVISIBLE(C)        ((C) < (Uchar) LOWESTVISIBLE ||\
                             (C) > (Uchar) HIGHESTVISIBLE)

/*
  Rescale characters denoted by numbers starting at 0 to
  the visible ASCII characters.
*/

#define VISIBLECHAR(I)      ((char)((I)+LOWESTVISIBLE))

/*
  Reverse the previous operation.
*/

#define INVISIBLECHAR(C)    ((Sint)((C)-LOWESTVISIBLE))

/*
  The following macro prints a character to a file pointer.
  If the character is not visible, then it is shown as the 
  corresponding ASCII-number with a prepended backslash.
*/

#define SHOWCHARFP(FP,C)\
        if(INVISIBLE(C))\
        {\
          fprintf(FP,"\\%lu",(Showuint) (C));\
        } else\
        {\
          (void) putc((Fputcfirstargtype) (C),FP);\
        }

/*
  The following macro is a variation of the previous macro,
  always showing the output to standard out.
*/

#define SHOWCHAR(C) SHOWCHARFP(stdout,C)

//\IgnoreLatex{

#endif

//}
