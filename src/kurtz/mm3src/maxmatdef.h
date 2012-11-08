/*
  Copyright (c) 2003 by Stefan Kurtz and The Institute for
  Genomic Research.  This is OSI Certified Open Source Software.
  Please see the file LICENSE for licensing information and
  the file ACKNOWLEDGEMENTS for names of contributors to the
  code base.
*/

//\IgnoreLatex{

#ifndef MAXMATDEF_H
#define MAXMATDEF_H
#include <limits.h>
#include "chardef.h"
#include "multidef.h"

//}

/*
  This file defines some constants and types for computing maximal 
  matches using suffix trees.
*/

/*
  The following two characters replace wildcard
  characters in a DNA sequences whenever the option \texttt{-n}
  (for match only nucleotides) is used. One is for the
  subject string and the other for the query string
*/

#define MMREPLACEMENTCHARSUBJECT (WILDCARD-2)
#define MMREPLACEMENTCHARQUERY   (WILDCARD-3)

/*
  The maximal number of query files.
*/

#define MAXNUMOFQUERYFILES       32

/*
  The following type contains all information
  derived from parsing the arguments of the program
*/

typedef struct
{
  BOOL showstring,              // show the matching string
       reversecomplement,       // compute matches on reverse strand
       forward,                 // compute matches on forward strand
       fourcolumn,              // always use 4 column format
       showreversepositions,    // give reverse pos. rel. to orig. string
       showsequencelengths,     // show length of sequences on header line
       matchnucleotidesonly,    // match ONLY acgt's
       cmaxmatch,               // compute all maximal matches
       cmumcand,                // compute reference-unique maximal matches
       cmum;                    // compute real matches unique in both sequences
  Uint minmatchlength,          // minimal length of a match to be reported
       numofqueryfiles;         // number of query files
  char program[PATH_MAX+1],     // the path of the program
       subjectfile[PATH_MAX+1], // filename of the subject-sequence
       queryfilelist[MAXNUMOFQUERYFILES][PATH_MAX+1]; 
                                // filenames of the query-sequences
} MMcallinfo;                   // \Typedef{MMcallinfo}

/*
  Functions processing a maximal match are of the following type.
*/

typedef Sint (*Processmatchfunction)
             (void *,Uint,Uint,Uint,Uint); // \Typedef{Processmatchfunction}

//\IgnoreLatex{

#endif

//}
