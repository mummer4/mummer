/*
  Copyright (c) 2003 by Stefan Kurtz and The Institute for
  Genomic Research.  This is OSI Certified Open Source Software.
  Please see the file LICENSE for licensing information and
  the file ACKNOWLEDGEMENTS for names of contributors to the
  code base.
*/

//\Ignore{

#include "types.h"
#include "mumcand.h"

//}

/*EE
  This module contains functions to extract from a table of
  MUM-candidates those which are also unique in the query
  sequence.
*/

/*
  The following function compares two MUM-candidates. The MUM-candidate
  with smaller \texttt{dbstart}-value comes first.
  If both MUMs have the same \texttt{dbstart}-value, then the MUM-candidate
  with the larger length comes first.
*/

static Sint compareMUMcandidates(MUMcandidate *p,MUMcandidate *q)
{
  if(p->dbstart == q->dbstart)
  {
    return (p->mumlength < q->mumlength) ? 1 : -1;
  }
  return (p->dbstart > q->dbstart) ? 1 : -1;
}

/*
  Sort all MUM-candidates according by increasing \texttt{dbstart}-value
  and decreasing length.
*/

static void sortMUMcandidates(ArrayMUMcandidate *mumcand)
{
  qsort(mumcand->spaceMUMcandidate,
        (size_t) mumcand->nextfreeMUMcandidate,
        sizeof(MUMcandidate),
        (Qsortcomparefunction) compareMUMcandidates);
}

/*EE
  Output all MUM candidates that are unique in the query sequence.
  These are the MUMs. The MUM-candidates are stored in table
  \texttt{mumcand}. The MUM is processed further by the function
  \texttt{processmum} which takes the \texttt{processinfo}
  as an argument.
*/

Sint mumuniqueinquery(void *processinfo,
                      Sint (*processmum)(void *,Uint,Uint,Uint,Uint),
                      ArrayMUMcandidate *mumcand)
{
  if(mumcand->nextfreeMUMcandidate > 0)
  {
    Uint currentright, dbright = 0;
    MUMcandidate *mumcandptr;
    BOOL ignorecurrent, ignoreprevious = False;

    sortMUMcandidates(mumcand);
    for(mumcandptr = mumcand->spaceMUMcandidate;
        mumcandptr < mumcand->spaceMUMcandidate + 
                     mumcand->nextfreeMUMcandidate;
        mumcandptr++)
    {
      ignorecurrent = False;
      currentright = mumcandptr->dbstart + mumcandptr->mumlength - 1;
      if(dbright > currentright)
      {
        ignorecurrent = True;
      } else
      {
        if(dbright == currentright)
        {
          ignorecurrent = True;
          if(!ignoreprevious && (mumcandptr-1)->dbstart == mumcandptr->dbstart)
          {
            ignoreprevious = True;
          }
        } else
        {
          dbright = currentright;
        }
      }
      if(mumcandptr > mumcand->spaceMUMcandidate && !ignoreprevious)
      {
        if(processmum(processinfo, 
                      (mumcandptr-1)->mumlength, 
                      (mumcandptr-1)->dbstart, 
                      (mumcandptr-1)->queryseq, 
                      (mumcandptr-1)->querystart) != 0)
        {
          return -1;
        }
      }
      ignoreprevious = ignorecurrent;
    }
    if(!ignoreprevious)
    {
      mumcandptr = mumcand->spaceMUMcandidate + 
                   mumcand->nextfreeMUMcandidate - 1;
      if(processmum(processinfo, 
                    mumcandptr->mumlength, 
                    mumcandptr->dbstart, 
                    mumcandptr->queryseq,
                    mumcandptr->querystart) != 0)
      {
        return -1;
      }
    }
  }
  return 0;
}
