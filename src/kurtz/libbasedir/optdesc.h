/*
  Copyright (c) 2003 by Stefan Kurtz and The Institute for
  Genomic Research.  This is OSI Certified Open Source Software.
  Please see the file LICENSE for licensing information and
  the file ACKNOWLEDGEMENTS for names of contributors to the
  code base.
*/

//\Ignore{

#ifndef OPTDESC_H
#define OPTDESC_H
#include "types.h"
#include "errordef.h"

//}

/*
  This file defines some macros to manipulate and a type to store
  information about options.
*/

/*
  The macro \texttt{OPTION} fills the component of a structure of 
  type \texttt{OptionDesc}. This macro should not be used. It
  is mainly for compatibility of the code.
*/

#define OPTION(S,N,O,D)\
        S[N].optname = O;\
        S[N].description = D;\
        S[N].optval = N;\
        S[N].isalreadyset = False;\
        S[N].declared = True

/*
  The macro \texttt{ADDOPTION} has a similar effect, but adding the 
  option is done via the function \texttt{addoption}. This is the
  recommended way to add an option.
*/

#define ADDOPTION(N,O,D)\
        if(addoption(&options[0],(Uint) NUMOFOPTIONS,(Uint) N,O,D) != 0)\
        {\
          return -1;\
        }

/*
  The following macro checks if an option is already set.
*/

#define ISSET(N)  (options[N].isalreadyset)

/*
  The following macro checks if an option is used. If not, then
  an error is reported.
*/

#define OPTIONMANDATORY(A)\
        if(!ISSET(A))\
        {\
          ERROR1("option %s is mandatory",options[A].optname);\
          return -1;\
        }

/*
  The following four macros check if all options \texttt{B} 
  (\texttt{C}, \texttt{D}, \texttt{E}) are used, whenever option 
  \texttt{A} is used.
*/

#define OPTIONIMPLY(A,B)\
        if(ISSET(A) && !ISSET(B))\
        {\
          ERROR2("option %s requires option %s",\
                  options[A].optname,options[B].optname);\
          return -1;\
        }

#define OPTIONIMPLYEITHER2(A,B,C)\
        if(ISSET(A))\
        {\
          if(!ISSET(B) && !ISSET(C))\
          {\
            ERROR3("option %s requires either option %s or %s",\
                    options[A].optname,options[B].optname,\
                    options[C].optname);\
            return -1;\
          }\
        }

#define OPTIONIMPLYEITHER3(A,B,C,D)\
        if(ISSET(A))\
        {\
          if(!ISSET(B) && !ISSET(C) && !ISSET(D))\
          {\
            ERROR4("option %s requires one of the options %s, %s, %s",\
                    options[A].optname,\
                    options[B].optname,\
                    options[C].optname,\
                    options[D].optname);\
            return -1;\
          }\
        }

#define OPTIONIMPLYEITHER4(A,B,C,D,E)\
        if(ISSET(A))\
        {\
          if(!ISSET(B) && !ISSET(C) && !ISSET(D) && !ISSET(E))\
          {\
            ERROR5("option %s requires one of the options %s, %s, %s, %s",\
                    options[A].optname,\
                    options[B].optname,\
                    options[C].optname,\
                    options[D].optname,\
                    options[E].optname);\
            return -1;\
          }\
        }

/*
  The following macro checks if two options are used simultaneously.
  If so, then an error is reported.
*/

#define OPTIONEXCLUDE(A,B)\
        if(ISSET(A) && ISSET(B))\
        {\
          ERROR2("option %s and option %s exclude each other",\
                  options[A].optname,options[B].optname);\
          return -1;\
        }

/*
  An option is described by the following type.
*/

typedef struct 
{
  char *optname,             // the option string, begins with -
       *description;         // help text describing purpose of option
  Uint optval;               // the unique number of an option
  BOOL isalreadyset,         // has the option already been set?
       declared;             // is the option declared by
                             // a line of the form \texttt{ADDOPTION}
} OptionDescription;         // \Typedef{OptionDescription}

//\Ignore{

#endif 

//}
