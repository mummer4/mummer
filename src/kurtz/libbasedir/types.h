/*
  Copyright (c) 2003 by Stefan Kurtz and The Institute for
  Genomic Research.  This is OSI Certified Open Source Software.
  Please see the file LICENSE for licensing information and
  the file ACKNOWLEDGEMENTS for names of contributors to the
  code base.
*/

//\IgnoreLatex{

#ifndef TYPES_H
#define TYPES_H
#include <sys/types.h>
#include <limits.h>
#include <stdlib.h>
#include <stdio.h>

#ifdef MSWINDOWS
typedef void * caddr_t;
#endif

/*
  Some rules about types:
  - do not use Ulong, these are not portable.
  - do not use the constants, UINT_MAX, INT_MAX and INT_MIN
  The following are the assumptions about the types:
  - size(Uint) >= 4
  - size(Sint) >= 4
  - size(Ushort) = 2
  - size(Sshort) = 2
  No other assumptions are to be made.
*/

//}

/*
  This file contains some basic type definition.
*/

typedef unsigned char  Uchar;         // \Typedef{Uchar}
typedef unsigned short Ushort;        // \Typedef{Ushort}

/*
  The following is the central case distinction to accomodate
  code for 32 bit integers and 64 bit integers.
*/

#ifdef SIXTYFOURBITS

typedef unsigned long  Uint;          // \Typedef{Uint}
typedef signed   long  Sint;          // \Typedef{Sint}
#define LOGWORDSIZE    6              // base 2 logarithm of wordsize
#define UintConst(N)   (N##UL)        // unsigned integer constant

#else

typedef unsigned int  Uint;          // \Typedef{Uint}
typedef signed   int  Sint;          // \Typedef{Sint}
#define LOGWORDSIZE    5             // base 2 logarithm of wordsize
#define UintConst(N)   (N##U)        // unsigned integer constant

#endif

/*
  Type of unsigned integer in \texttt{printf}.
*/

typedef unsigned long Showuint;     // \Typedef{Showuint}

/*
  Type of signed integer in \texttt{printf}.
*/

typedef signed long Showsint;       // \Typedef{Showsint}

/*
  Type of integer in \texttt{scanf}.
*/

typedef signed long Scaninteger;    // \Typedef{Scaninteger}

/*
  Argument of a function from \texttt{ctype.h}.
*/

typedef int Ctypeargumenttype;      // \Typedef{Ctypeargumenttype}

/*
  Return type of \texttt{fgetc} and \texttt{getc}.
*/

typedef int Fgetcreturntype;        // \Typedef{Fgetcreturntype}

/*
  Type of first argument of \texttt{fputc}.
*/

typedef int Fputcfirstargtype;      // \Typedef{Fputsfirstargtype}  

/*
  Return type of \texttt{strcmp}.
*/

typedef int Strcmpreturntype;       // \Typedef{Strcmpreturntype}

/*
  Type of a file descriptor.
*/

typedef int Filedesctype;           // \Typedef{Filedesctype}

/*
  Return type of \texttt{qsort} function.
*/

typedef int Qsortcomparereturntype; // \Typedef{Qsortcomparefunction}

/*
  Return type of \texttt{sprintf} function.
*/

typedef int Sprintfreturntype;     // \Typedef{Sprintfreturntype} 

/*
  Type of fieldwidth in \texttt{printf} format string.
*/

typedef int Fieldwidthtype;         // \Typedef{Fieldwidthtype}

/*
  Type of \texttt{argc}-parameter in main.
*/

typedef int Argctype;               // \Typedef{Argctype}

/*
  Return type of \texttt{getrlimit}
*/

typedef int Getrlimitreturntype;    // \Typedef{Getrlimitreturntype}

#ifdef WITHSYSCONF
typedef int Sysconfargtype;         // \Typedef{Sysconfargtype}
#endif

/*
  The following macros define some basic division, multiplication,
  and modulo operations on unsigned integers.
*/

#define DIV2(N)      ((N) >> 1)
#define DIV4(N)      ((N) >> 2)
#define DIV8(N)      ((N) >> 3)
#define MULT2(N)     ((N) << 1)
#define MULT4(N)     ((N) << 2)
#define MULT8(N)     ((N) << 3)
#define MOD2(N)      ((N) & 1)
#define MOD4(N)      ((N) & 3)
#define MOD8(N)      ((N) & 7)

//\IgnoreLatex{

#define CHECKTYPESIZE(T,OP,S)\
        if(sizeof(T) OP (S))\
        {\
          DEBUG4(1,"# sizeof(%s) %s (%ld bytes,%ld bits) as epected\n",\
                  #T,#OP,(Showsint) sizeof(T),\
                         (Showsint) (CHAR_BIT * sizeof(T)));\
        } else\
        {\
          fprintf(stderr,"typesize constraint\n");\
          fprintf(stderr,"  sizeof(%s) = (%ld bytes,%ld bits) %s %lu bytes\n",\
                  #T,\
                  (Showsint) sizeof(T),\
                  (Showsint) (CHAR_BIT * sizeof(T)),\
                  #OP,\
                  (Showuint) (S));\
          fprintf(stderr,"does not hold\n");\
          exit(EXIT_FAILURE);\
        }

/*
  The following function checks some type constraints 
*/

#define CHECKALLTYPESIZES\
        CHECKTYPESIZE(char,==,(size_t) 1)\
        CHECKTYPESIZE(short,==,(size_t) 2)\
        CHECKTYPESIZE(int,==,(size_t) 4)\
        CHECKTYPESIZE(long,>=,(size_t) 4)\
        CHECKTYPESIZE(void *,>=,(size_t) 4)

//}

/*
  Here is a prototype for the main function.
*/

#define MAINFUNCTION int main(Argctype argc,char *argv[])

//\IgnoreLatex{

#ifndef __cplusplus
int mkstemp(char *);
#endif

//}

/*
  A type for boolean values defined as a constant to allow 
  checking if it has been defined previously.
*/

#ifndef BOOL
#define BOOL unsigned char
#endif

#ifndef False
#define False ((BOOL) 0)
#endif

#ifndef True
#define True ((BOOL) 1)
#endif

/*
  Show a boolean value as a string or as a character 0 or 1.
*/

#define SHOWBOOL(B) ((B) ? "True" : "False")
#define SHOWBIT(B)  ((B) ? '1' : '0')

/*
  Pairs, triples, and quadruples of unsigned integers.
*/

typedef struct 
{
  Uint uint0, uint1;
} PairUint;                // \Typedef{PairUint}

typedef struct 
{
  Uint uint0, uint1, uint2;
} ThreeUint;               // \Typedef{ThreeUint}

typedef struct 
{
  Uint uint0, uint1, uint2, uint3;
} FourUint;                // \Typedef{FourUint}

//\IgnoreLatex{

/*
  A list is stored with its start position in some space block 
  and its length.
*/

typedef struct 
{
  Uint start, length;
} Listtype;                // \Typedef{Listtype}

/*
  A string is just a list.
*/

typedef Listtype Stringtype;    // \Typedef{Stringtype}

/*
  The default type for length-values is unsigned int.
*/

#ifndef LENGTHTYPE
#define LENGTHTYPE Uint
#endif

/*
  The default number of bytes in a bitvector used for dynamic programming
  is 4.
*/

#ifndef DPBYTESINWORD
#define DPBYTESINWORD 4
#endif

/*
  The number of bytes in a dynamic programming bitvector determines the type
  of integers, the dp-bits are stored in.
*/

#if DPBYTESINWORD == 1
typedef unsigned char DPbitvector;          // \Typedef{DPbitvector}
#else
#if DPBYTESINWORD == 2
typedef unsigned short DPbitvector;
#else
#if DPBYTESINWORD == 4
typedef unsigned int DPbitvector;
#else
#if DPBYTESINWORD == 8
typedef unsigned long long DPbitvector;
#endif
#endif
#endif
#endif

typedef unsigned int DPbitvector4;          // \Typedef{DPbitvector4}

#if (LOGWORDSIZE==6)
typedef unsigned long DPbitvector8;         // \Typedef{DPbitvector8}
#endif

//}

/*
  The following type stores filenames and the length of the corresponding
  files.
*/

typedef struct
{
  char *filenamebuf;    // pointer to a copy of a filename
  Uint filelength;      // the length of the corresponding file
} Fileinfo;             // \Typedef{Fileinfo}

/*
  The following is the type of the comparison function
  to be provided to the function \texttt{qsort}.
*/

typedef int (*Qsortcomparefunction)(const void *,const void *);

//\IgnoreLatex{

#endif

//}
