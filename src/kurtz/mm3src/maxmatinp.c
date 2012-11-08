/*
  Copyright (c) 2003 by Stefan Kurtz and The Institute for
  Genomic Research.  This is OSI Certified Open Source Software.
  Please see the file LICENSE for licensing information and
  the file ACKNOWLEDGEMENTS for names of contributors to the
  code base.
*/

//\IgnoreLatex{

#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include "chardef.h"
#include "spacedef.h"
#include "protodef.h"
#include "debugdef.h"
#include "maxmatdef.h"

//}

/*EE
  This module contains all functions parsing the subject and 
  query sequences.
*/

/*EE
  For each sequence in the \texttt{Multiseq}-record, the dynamic
  array \texttt{startdesc} stores the positions in the dynamic array
  \texttt{descspace} where the sequence description starts. The 
  following macro checks if enough memory has been allocated
  for \texttt{startdesc}. If not, then this is done by incrementing
  the size of the array by 128 entries. Finally, the appropriate
  entry in \texttt{startdesc} is assigned the correct value.
*/

#define STORESTARTDESC\
        if(multiseq->numofsequences >= allocatedstartdesc)\
        {\
          allocatedstartdesc += 128;\
          multiseq->startdesc\
            = ALLOCSPACE(multiseq->startdesc,Uint,allocatedstartdesc);\
        }\
        multiseq->startdesc[multiseq->numofsequences]\
          = multiseq->descspace.nextfreeUchar

/*EE
  The following function scans a string containing the content of
  a multiple fasta formatted file. The parameter are as follows:
  \begin{enumerate}
  \item
  \texttt{multiseq} is the \texttt{Multiseq}-record to store the
  scanned information in.
  \item
  \texttt{filename} is the information from which the file
  contents was read.
  \item
  \texttt{replacewildcardchar} is the character used to
  replace a wildcard (then it should be different from the
  characters occuring in DNA sequences) or 0 if wildcards are not 
  replaced. 
  \item
  \texttt{input} points to the inputstring to be scanned,
  \item
  \texttt{inputlen} is the length of the input.
  \end{enumerate}
  Each sequence description begins with the symbol 
  \texttt{>}.  If it does, then this symbol is skipped. The rest of the 
  line up to the first white space character is stored in
  \texttt{descspace}. Otherwise, the rest of the line is discarded. 
  The remaining lines (until the next symbol \texttt{>} or
  the end of the input string) are scanned for alphanumeric
  characters which make up the sequence. White spaces are ignored.
  Upper case characters are transformed to lower case.
  The input string must contain at least one sequence.
  In case of a error, an negative error code is returned.
  In case of success, the return code is 0.
*/

Sint scanmultiplefastafile (Multiseq *multiseq,
                            char *filename,
                            Uchar replacewildcardchar,
                            Uchar *input,
                            Uint inputlen)
{
  Uchar *inputptr,      // points to a suffix of the input
        *newptr,        // points to the transform,ed
        tmpchar;        // temporary character
  Uint allocatedstartdesc = 0;  // num of characters allocated for startdesc
  BOOL indesc = False,          // inside description part of sequence
       copydescription = False; // currently copying the description 

  fprintf(stderr,"# reading input file \"%s\" ",filename);
  initmultiseq (multiseq);
  multiseq->originalsequence = NULL;

  newptr = multiseq->sequence = input;
  for (inputptr = input; inputptr < input + inputlen; inputptr++)
  {
    if (indesc)
    {
      if(copydescription)
      {
        if(isspace((Ctypeargumenttype) *inputptr))
        {
          copydescription = False;
          STOREINARRAY (&multiseq->descspace, 
                        Uchar,
                        4096,
                        (Uchar) '\n');
        } else
        {
          STOREINARRAY (&multiseq->descspace, 
                        Uchar,
                        4096,
                        *inputptr);
        }
      } 
      if (*inputptr == '\n')
      {
        indesc = False;
      }
    }
    else
    {
      if (*inputptr == FASTASEPARATOR)
      {
        STORESTARTDESC;
        if (multiseq->numofsequences > 0)
        {
          STOREINARRAY (&multiseq->markpos, 
                        Uint,
                        128,
                        (Uint) (newptr - multiseq->sequence));

          *newptr++ = SEPARATOR;
        }
        multiseq->numofsequences++;
        indesc = True;
        copydescription = True;
      }
      else
      {
        tmpchar = *inputptr;
        if (!isspace ((Ctypeargumenttype) tmpchar))
        {
          tmpchar = (Uchar) tolower ((Ctypeargumenttype) tmpchar);
          if (replacewildcardchar != 0)  // replace wildcards
          {
            switch (tmpchar)
            {
              case 'a':
              case 'c':
              case 'g':
              case 't':
                break;
              default:
                /*
                   fprintf(stderr,"filename %s: replace '%c' by '%c'\n",
                   filename,tmpchar,replacewildcardchar);
                 */
                tmpchar = replacewildcardchar;
            }
          }
#ifdef WARNINGIFNONUCLEOTIDES
          else
          {
            switch (tmpchar)
            {
              case 'a':
              case 'c':
              case 'g':
              case 't':
              case 's':
              case 'w':
              case 'r':
              case 'y':
              case 'm':
              case 'k':
              case 'b':
              case 'd':
              case 'h':
              case 'v':
              case 'n':
                break;
              default:
                fprintf (stderr,
                         "Unexpected character '%c\' in string %s\n",
                         tmpchar, filename);
                tmpchar = 'n';
            }
          }
#endif
          *newptr++ = tmpchar;
        }
      }
    }
  }
  STORESTARTDESC;
  if (multiseq->numofsequences == 0)
  {
    ERROR0 ("no sequences in multiple fasta file");
    return -2;
  }
  multiseq->totallength = (Uint) (newptr - multiseq->sequence);
  if(multiseq->totallength == 0)
  {
    ERROR0 ("empty sequence in multiple fasta file");
    return -3;
  }
  fprintf(stderr,"of length %lu\n",(Showuint) multiseq->totallength);
  return 0;
}

/*EE
  The following function reads the subject and queryfile and
  delivers the parsed multiple sequences in the corresponding
  \texttt{Multiseq}-records. The files are read via memory mapping.
  The subject file must contain exactly one sequence.
  Both files cannot be empty. The parameter \texttt{matchnucleotidesonly}
  is true iff if the programm was called with option \texttt{-n},
  which means that only \texttt{ACGTs} are matched. If an error occurs, 
  then the function delivers a negative error code. Otherwise the error
  code is 0.
*/

Sint getmaxmatinput (Multiseq *subjectmultiseq,
                     BOOL matchnucleotidesonly,
                     char *subjectfile)
{
  Uint filelen;
  Uchar *filecontent;

  filecontent = CREATEMEMORYMAP (subjectfile, True, &filelen);
  if (filecontent == NULL || filelen == 0)
  {
    ERROR2("cannot open file \"%s\" or file \"%s\" is empty",subjectfile,
                                                             subjectfile);
    return -1;
  }
  if (scanmultiplefastafile (subjectmultiseq, subjectfile,
                             matchnucleotidesonly ? MMREPLACEMENTCHARSUBJECT
                                                  : 0,
                             filecontent, 
                             filelen) != 0)
  {
    return -2;
  }
  DEBUG1 (2, "subject of length %lu=",
          (Showuint) subjectmultiseq->totallength);
  DEBUGCODE (2, (void) fwrite (subjectmultiseq->sequence, sizeof (Uchar),
                               (size_t) subjectmultiseq->totallength, stdout));
  DEBUG0 (2, "\n");
  return 0;
}
