/*
  Copyright (c) 2003 by Stefan Kurtz and The Institute for
  Genomic Research.  This is OSI Certified Open Source Software.
  Please see the file LICENSE for licensing information and
  the file ACKNOWLEDGEMENTS for names of contributors to the
  code base.
*/

//\Ignore{

#ifndef INTBITS_H
#define INTBITS_H
#include <limits.h>
#include "types.h"
#include "errordef.h"
#include "spacedef.h"

//}

/*
  This file contains some definitions manipulating bitvectors represented
  by a \texttt{Uint}. In the comment lines we use $w$ for the word size
  and \texttt{\symbol{94}} for exponentiation of the previous character.
*/

#define INTWORDSIZE\
        (UintConst(1) << LOGWORDSIZE)     // # of bits in Uint = w
#define FIRSTBIT\
        (UintConst(1) << (INTWORDSIZE-1)) // \(10^{w-1}\)
#define ISBITSET(S,I)\
        (((S) << (I)) & FIRSTBIT)         // is \(i\)th bit set?
#define ITHBIT(I)\
        (FIRSTBIT >> (I))                 // \(0^{i}10^{w-i-1}\) 
#define SECONDBIT\
        (FIRSTBIT >> 1)                   // \(010^{w-2}\)
#define THIRDBIT\
        (FIRSTBIT >> 2)                   // \(0010^{w-3}\)
#define FIRSTTWOBITS\
        (UintConst(3) << (INTWORDSIZE-2)) // \(11^{w-2}\)
#define EXCEPTFIRSTBIT\
        (~FIRSTBIT)                       // \(01^{w-1}\)
#define EXCEPTFIRSTTWOBITS\
        (EXCEPTFIRSTBIT >> 1)             // \(001^{w-2}\)
#define EXCEPTFIRSTTHREEBITS\
        (EXCEPTFIRSTBIT >> 2)             // \(0001^{w-3}\)
#define DIVWORDSIZE(I)\
        ((I) >> LOGWORDSIZE)              // \((I) div w\)
#define MODWORDSIZE(I)\
        ((I) & (INTWORDSIZE-1))           // \((I) mod w\)
#define MULWORDSIZE(I)\
        ((I) << LOGWORDSIZE)              // \((I) * w\)

/*
  The following macro allocates a bitarray of \texttt{N} bits. All bits
  are off.
*/

#define INITBITTAB(TAB,N)\
        {\
          Uint *tabptr, tabsize = 1 + DIVWORDSIZE(N);\
          TAB = ALLOCSPACE(NULL,Uint,tabsize);\
          for(tabptr = TAB; tabptr < (TAB) + tabsize; tabptr++)\
          {\
            *tabptr = 0;\
          }\
        }

/*
  The following macro inititalizes a bitarray such tha all bits
  are off.
*/

#define CLEARBITTAB(TAB,N)\
        {\
          Uint *tabptr, tabsize = 1 + DIVWORDSIZE(N);\
          for(tabptr = TAB; tabptr < TAB + tabsize; tabptr++)\
          {\
            *tabptr = 0;\
          }\
        }


/*
  \texttt{SETIBIT(TAB,I)} sets the \texttt{I}-th bit in bitarray 
  \texttt{TAB} to 1.
*/

#define SETIBIT(TAB,I)    (TAB)[DIVWORDSIZE(I)] |= ITHBIT(MODWORDSIZE(I))

/*
  \texttt{UNSSETIBIT(TAB,I)} sets the \texttt{I}-th bit in bitarray 
  \texttt{TAB} to 0.
*/

#define UNSETIBIT(TAB,I)  (TAB)[DIVWORDSIZE(I)] &= ~(ITHBIT(MODWORDSIZE(I)))

/*
  \texttt{ISIBITSET(TAB,I)} checks if the \texttt{I}-th bit in bitarray 
  \texttt{TAB} is 1.
*/

#define ISIBITSET(TAB,I)  ((TAB)[DIVWORDSIZE(I)] & ITHBIT(MODWORDSIZE(I)))

//\Ignore{

#ifdef __cplusplus
  extern "C" {
#endif
  
char *intbits2string(Uint bs);

#ifdef __cplusplus
}
#endif

#endif

//}
