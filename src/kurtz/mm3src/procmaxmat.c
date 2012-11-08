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
#include <string.h>
#include <ctype.h>
#include "types.h"
#include "debugdef.h"
#include "errordef.h"
#include "protodef.h"
#include "spacedef.h"
#include "streedef.h"
#include "maxmatdef.h"

//}

/*EE
  This file contains functions to appropriately
  call the function \texttt{findmumcandidates} and \texttt{findmaxmatches}
  and to process their result according to the options given by the user.
*/

/*EE
  The following structure contains all information
  required while computing and processing the matches.
*/

typedef struct
{
  Suffixtree stree;            // the suffix tree of the subject-sequence
  Multiseq *subjectmultiseq,   // reference to multiseq of subject
           querymultiseq;      // the Multiseq record of the queries
  ArrayMUMcandidate mumcandtab;// a table containing MUM-candidates
                               // when option \texttt{-mum} is on
  Uint minmatchlength,         // minimum length of a match
       maxdesclength,          // maximum length of a description
       currentquerylen;        // length of the current query sequence
  BOOL showstring,             // is option \texttt{-s} on?
       showsequencelengths,    // is option \texttt{-L} on?
       showreversepositions,   // is option \texttt{-c} on?
       forward,                // compute forward matches
       fourcolumn,             // is option \texttt{-F} on?
       reversecomplement,      // compute reverse complement matches
       cmumcand,               // compute MUM candidates
       cmum,                   // compute MUMs
       currentisrcmatch;       // true iff currently rc-matches are computed
} Matchprocessinfo;            // \Typedef{Matchprocessinfo}

//\IgnoreLatex{

/*
  The following is the type for the functions to finding maximal matches or
  MUM candidates.
*/

typedef Sint (*Findmatchfunction)(Suffixtree *,
                                  Uint,
                                  Processmatchfunction,
                                  void *,
                                  Uchar *,
                                  Uint,
                                  Uint);

/*
  The following function is imported from \texttt{findmumcand.c}.
*/

Sint findmumcandidates(Suffixtree *stree,
                       Uint minmatchlength,
                       Processmatchfunction processmatch,
                       void *processinfo,
                       Uchar *query,
                       Uint querylen,
                       Uint seqnum);

/*
  The following function is imported from \texttt{findmaxmat.c}
*/

Sint findmaxmatches(Suffixtree *stree,
                    Uint minmatchlength,
                    Processmatchfunction processmatch,
                    void *processinfo,
                    Uchar *query,
                    Uint querylen,
                    Uint seqnum);

/*
  The following function is imported from \texttt{maxmatinp.c}
*/

Sint scanmultiplefastafile (Multiseq *multiseq,
                            char *filename,
                            Uchar replacewildcardchar,
                            Uchar *input,
                            Uint inputlen);

//}

/*EE
  The following code fragement is used when computing MUMs.
  It calls the function \texttt{mumuniqueinquery}.
  The MUM-candidates are stored in the dynamic 
  array \texttt{mumcandtab}.  After the real MUMs are output, the table 
  of MUM-candidates are declared to be empty.
*/

#define PROCESSREALMUMS\
        if(matchprocessinfo->cmum)\
        {\
          if(mumuniqueinquery(info,\
                              matchprocessinfo->showstring ?\
                                 showseqandmaximalmatch :\
                                 showmaximalmatch,\
                              &matchprocessinfo->mumcandtab) != 0)\
          {\
            return -2;\
          }\
          matchprocessinfo->mumcandtab.nextfreeMUMcandidate = 0;\
        }

/*
  The following function assigns for a given base the corresponding
  complement character. It also handles wild card characters 
  appropriately.
*/

#define ASSIGNMAXMATCOMPLEMENT(VL,VR)\
        if((VR) >= MMREPLACEMENTCHARQUERY)\
        {\
          VL = VR;\
        } else\
        {\
          switch (VR)\
          {\
            case 'a':\
              VL = (Uchar) 't';\
              break;\
            case 'c':\
              VL = (Uchar) 'g';\
              break;\
            case 'g':\
              VL = (Uchar) 'c';\
              break;\
            case 't':\
              VL = (Uchar) 'a';\
              break;\
            case 'r':                       /* a or g */\
              VL = (Uchar) 'y';\
              break;\
            case 'y':                       /* c or t */\
              VL = (Uchar) 'r';\
              break;\
            case 's':                       /* c or g */\
              VL = (Uchar) 's';\
              break;\
            case 'w':                       /* a or t */\
              VL = (Uchar) 'w';\
              break;\
            case 'm':                       /* a or c */\
              VL = (Uchar) 'k';\
              break;\
            case 'k':                       /* g or t */\
              VL = (Uchar) 'm';\
              break;\
            case 'b':                       /* c, g or t */\
              VL = (Uchar) 'v';\
              break;\
            case 'd':                       /* a, g or t */\
              VL = (Uchar) 'h';\
              break;\
            case 'h':                       /* a, c or t */\
              VL = (Uchar) 'd';\
              break;\
            case 'v':                       /* a, c or g */\
              VL = (Uchar) 'b';\
              break;\
            default:                        /* anything */\
              VL = (Uchar) 'n';\
              break;\
          }\
        }

/*
  The following function computes the Watson-Crick complement
  of the sequence pointed to by \texttt{seq}. The sequnce is
  of length \texttt{seqlen}. The result is computed in-place.
*/

static void wccSequence (Uchar *seq,
                         Uint seqlen)
{
  Uchar *front, *back, tmp = 0;

  for (front = seq, back = seq + seqlen - 1; 
       front <= back; front++, back--)
  {
    ASSIGNMAXMATCOMPLEMENT (tmp, *front);
    ASSIGNMAXMATCOMPLEMENT (*front, *back);
    *back = tmp;
  }
}

/*
  The following function shows the sequence description of sequence
  number \texttt{seqnum} in \texttt{multiseq}.
*/

static void showsequencedescription(Multiseq *multiseq, Uint maxdesclength,
                                    Uint seqnum)
{
  Uint i, desclength = DESCRIPTIONLENGTH(multiseq,seqnum);
  Uchar *desc = DESCRIPTIONPTR(multiseq,seqnum);

  for(i=0; i<desclength; i++)
  {
    if(isspace((Ctypeargumenttype) desc[i]))
    {
      break;
    }
    (void) putchar((Fputcfirstargtype) desc[i]);
  }
  if(desclength < maxdesclength)
  {
    for(i=0; i < maxdesclength - desclength; i++)
    {
      (void) putchar(' ');
    }
  }
}

/*
  The following function shows the sequence header for the
  sequence with number \texttt{seqnum} in the \texttt{Multiseq}-record. 
  \texttt{showsequencelengths} is true if the option \texttt{-L} was
  set. \texttt{currentisrcmatch} is True iff if the current matches
  to be reported (if any) are reverse complemented matches.
  \texttt{seqlen} is the length of the sequence for which the header
  is shown.
*/

static void showsequenceheader(Multiseq *multiseq,
                               BOOL showsequencelengths,
                               BOOL currentisrcmatch,
                               Uint seqnum,
                               Uint seqlen)
{

  printf("%c ",FASTASEPARATOR);
  showsequencedescription(multiseq,0,seqnum);
  if(currentisrcmatch)
  {
    printf(" Reverse");
  }
  if(showsequencelengths)
  {
    printf("  Len = %lu",(Showuint) seqlen);
  }
  printf("\n");
}

/*
  The following is one of three functions to further process 
  a MUM-candidate specfied by its length, its start position
  in the subject sequence, the number of the query in which
  it matches as well as the start of the match in the query.
  The function \texttt{showmaximalmatch} simply shows the 
  relevant information as a triple of three integers.
*/

static Sint showmaximalmatch (void *info,
                              Uint matchlength,
                              Uint subjectstart,
                              /*@unused@*/ Uint seqnum,
                              Uint querystart)
{
  Matchprocessinfo *matchprocessinfo = (Matchprocessinfo *) info;

  if(matchprocessinfo->subjectmultiseq->numofsequences == UintConst(1)
     &&
     !matchprocessinfo->fourcolumn)
  {
    printf ("%8lu  ", (Showuint) (subjectstart+1));
  } else
  {
    PairUint pp;

    if(pos2pospair(matchprocessinfo->subjectmultiseq,&pp,subjectstart) != 0)
    {
      return -1;
    }
    printf("  ");
    showsequencedescription(matchprocessinfo->subjectmultiseq,
                            matchprocessinfo->maxdesclength,
                            pp.uint0);
    printf ("  %8lu  ",(Showuint) (pp.uint1+1));
  }
  if(matchprocessinfo->currentisrcmatch && 
     matchprocessinfo->showreversepositions)
  {
    printf ("%8lu  ", (Showuint) (matchprocessinfo->currentquerylen - 
                                  querystart));
  } else
  {
    printf ("%8lu  ", (Showuint) (querystart+1));
  }
  printf ("%8lu\n", (Showuint) matchlength);
  return 0;
}

/*
  The following function shows the same information as the previous
  function and additionally the sequence information.
*/

static Sint showseqandmaximalmatch (void *info,
                                    Uint matchlength,
                                    Uint subjectstart,
                                    Uint seqnum,
                                    Uint querystart)
{
  Matchprocessinfo *matchprocessinfo = (Matchprocessinfo *) info;

  (void) showmaximalmatch (info,
                           matchlength,
                           subjectstart,
                           seqnum,
                           querystart);
  if (fwrite (matchprocessinfo->subjectmultiseq->sequence + subjectstart, 
              sizeof (Uchar), 
              (size_t) matchlength,
              stdout) != (size_t) matchlength)
  {
    ERROR1 ("cannot output string of length %lu", (Showuint) matchlength);
    return -1;
  }
  printf ("\n");
  return 0;
}

/*
  The following function stores the information about a MUM-candidate
  in the next free position of the dynamic array \texttt{mumcandtab}.
*/

static Sint storeMUMcandidate (void *info,
                               Uint matchlength,
                               Uint subjectstart,
                               Uint seqnum,
                               Uint querystart)
{
  Matchprocessinfo *matchprocessinfo = (Matchprocessinfo *) info;
  MUMcandidate *mumcandptr;

  DEBUG4(3,"storeMUMcandiate %lu %lu %lu %lu\n",
            (Showuint) matchlength,
            (Showuint) subjectstart,
            (Showuint) seqnum,
            (Showuint) querystart);
  GETNEXTFREEINARRAY(mumcandptr,
                     &matchprocessinfo->mumcandtab,
                     MUMcandidate,1024);
  mumcandptr->mumlength = matchlength;
  mumcandptr->dbstart = subjectstart;
  mumcandptr->queryseq = seqnum;
  mumcandptr->querystart = querystart;
  return 0;
}

/*
  The following function searches for forward and reverse complemented
  MUM-candidates (if necessary) in the current query of length
  \texttt{querylen}. The number of the query sequence is 
  \texttt{querylen}.
*/

static Sint findmaxmatchesonbothstrands(void *info,Uint seqnum,
                                        Uchar *query,Uint querylen)
{
  Matchprocessinfo *matchprocessinfo = (Matchprocessinfo *) info;
  Processmatchfunction processmatch;
  Findmatchfunction findmatchfunction;

  if(matchprocessinfo->cmum)
  {
    processmatch = storeMUMcandidate;
    findmatchfunction = findmumcandidates;
  } else
  {
    if(matchprocessinfo->showstring)
    {
      processmatch = showseqandmaximalmatch;
    } else
    {
      processmatch = showmaximalmatch;
    }
    if(matchprocessinfo->cmumcand)
    {
      findmatchfunction = findmumcandidates;
    } else
    {
      findmatchfunction = findmaxmatches;
    }
  }
  matchprocessinfo->currentquerylen = querylen;
  if(matchprocessinfo->forward)
  {
    showsequenceheader(&matchprocessinfo->querymultiseq,
                       matchprocessinfo->showsequencelengths,
                       False,
                       seqnum,
                       querylen);
    matchprocessinfo->currentisrcmatch = False;
    if(findmatchfunction(&matchprocessinfo->stree,
                         matchprocessinfo->minmatchlength,
                         processmatch,
                         info,
                         query,
                         querylen,
                         seqnum) != 0)
    {
      return -1;
    }
    PROCESSREALMUMS;
  }
  if(matchprocessinfo->reversecomplement)
  {
    if(matchprocessinfo->cmum)
    {
      matchprocessinfo->mumcandtab.nextfreeMUMcandidate = 0;
    }
    showsequenceheader(&matchprocessinfo->querymultiseq,
                       matchprocessinfo->showsequencelengths,
                       True,
                       seqnum,
                       querylen);
    wccSequence(query,querylen);
    matchprocessinfo->currentisrcmatch = True;
    if(findmatchfunction(&matchprocessinfo->stree,
                         matchprocessinfo->minmatchlength,
                         processmatch,
                         info,
                         query,
                         querylen,
                         seqnum) != 0)
    {
      return -2;
    }
    PROCESSREALMUMS;
  }
  return 0;
}

static Sint getmaxdesclen(Multiseq *multiseq)
{
  Uint desclen, maxdesclen, seqnum;
  if(multiseq->numofsequences == 0)
  {
    ERROR0("multiple sequence contains 0 sequences");
    return -1;
  }
  maxdesclen = DESCRIPTIONLENGTH(multiseq,0);
  for(seqnum = UintConst(1); seqnum < multiseq->numofsequences; seqnum++)
  {
    desclen = DESCRIPTIONLENGTH(multiseq,seqnum);
    if(desclen > maxdesclen)
    {
      maxdesclen = desclen;
    }
  }
  return (Sint) maxdesclen;
}

/*EE
  The following function constructs the suffix tree,
  initializes the \texttt{Matchprocessinfo}-record appropriately,
  initializes the dynamic array \texttt{mumcandtab} (if necessary),
  and then iterates the function \texttt{findmaxmatchesonbothstrands}
  over all sequences in \texttt{querymultiseq}. Finally, the space
  allocated for the suffix tree and the space for \texttt{mumcandtab}
  is freed.
*/

Sint procmaxmatches(MMcallinfo *mmcallinfo,Multiseq *subjectmultiseq)
{
  Matchprocessinfo matchprocessinfo;
  Uint filenum, filelen;
  Sint retcode;
  Uchar *filecontent;

  fprintf(stderr,"# construct suffix tree for sequence of length %lu\n",
           (Showuint) subjectmultiseq->totallength);
  fprintf(stderr,"# (maximum reference length is %lu)\n",
           (Showuint) getmaxtextlenstree());
  fprintf(stderr,"# (maximum query length is %lu)\n",
          (Showuint) ~((Uint)0));
  if(constructprogressstree (&matchprocessinfo.stree,
                             subjectmultiseq->sequence,
                             subjectmultiseq->totallength,
                             NULL,
                             NULL,
                             NULL) != 0)
  {
    return -1;
  }
  fprintf(stderr,"# CONSTRUCTIONTIME %s %s %.2f\n",
         &mmcallinfo->program[0],&mmcallinfo->subjectfile[0],
         getruntime());
  matchprocessinfo.subjectmultiseq = subjectmultiseq;
  matchprocessinfo.minmatchlength = mmcallinfo->minmatchlength;
  matchprocessinfo.showstring = mmcallinfo->showstring;
  matchprocessinfo.showsequencelengths = mmcallinfo->showsequencelengths;
  matchprocessinfo.showreversepositions = mmcallinfo->showreversepositions;
  matchprocessinfo.forward = mmcallinfo->forward;
  matchprocessinfo.fourcolumn = mmcallinfo->fourcolumn;
  matchprocessinfo.cmum = mmcallinfo->cmum;
  matchprocessinfo.cmumcand = mmcallinfo->cmumcand;
  matchprocessinfo.reversecomplement = mmcallinfo->reversecomplement;
  if(mmcallinfo->cmum)
  {
    INITARRAY(&matchprocessinfo.mumcandtab,MUMcandidate);
  }
  retcode = getmaxdesclen(subjectmultiseq);
  if(retcode < 0)
  {
    return -2;
  }
  matchprocessinfo.maxdesclength = (Uint) retcode;
  for(filenum=0; filenum < mmcallinfo->numofqueryfiles; filenum++)
  {
    filecontent = CREATEMEMORYMAP (mmcallinfo->queryfilelist[filenum],
                                   True, 
                                   &filelen);
    if (filecontent == NULL || filelen == 0)
    {
      ERROR2("cannot open file \"%s\" or file \"%s\" is empty",
              mmcallinfo->queryfilelist[filenum],
              mmcallinfo->queryfilelist[filenum]);
      return -3;
    }
    if (scanmultiplefastafile (&matchprocessinfo.querymultiseq, 
                               mmcallinfo->queryfilelist[filenum],
                               mmcallinfo->matchnucleotidesonly 
                                ? MMREPLACEMENTCHARQUERY 
                                : 0,
                               filecontent, filelen) != 0)
    {
      return -4;
    }
    fprintf(stderr,
	    "# matching query-file \"%s\"\n# against subject-file \"%s\"\n",
            mmcallinfo->queryfilelist[filenum],
            mmcallinfo->subjectfile);
    if (overallsequences (False,
                          &matchprocessinfo.querymultiseq,
                          (void *) &matchprocessinfo,
                          findmaxmatchesonbothstrands) != 0)
    {
      return -5;
    }
    freemultiseq(&matchprocessinfo.querymultiseq);
  }
  if(mmcallinfo->cmum)
  {
    FREEARRAY(&matchprocessinfo.mumcandtab,MUMcandidate);
  }
  freestree (&matchprocessinfo.stree);
  return 0;
}
