/*
  Copyright (c) 2003 by Stefan Kurtz and The Institute for
  Genomic Research.  This is OSI Certified Open Source Software.
  Please see the file LICENSE for licensing information and
  the file ACKNOWLEDGEMENTS for names of contributors to the
  code base.
*/

//\Ignore{

#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <string.h>
#include <limits.h>
#include <unistd.h>
#include "types.h"
#include "fopen.h"
#include "debugdef.h"
#include "errordef.h"
#include "megabytes.h"

//}

/*EE
  This module defines functions for handling debug levels and
  other information related to producing debugging messages.
  The debug mechanism is only available if the \texttt{DEBUG} compiler 
  flag is used.
*/

#ifdef DEBUG

static Sint debuglevel = 0;        // the value of \texttt{DEBUGLEVEL}
static BOOL debugwhere = False;    // the value of \texttt{DEBUGWHERE}
/*@null@*/ static FILE 
           *debugfileptr = NULL;  // the file pointer to show the debug info

/*EE
  The following function returns the \texttt{DEBUGLEVEL}.
*/

Sint getdebuglevel(void)
{
  return debuglevel;
}

/*EE
  The following function returns the value of \texttt{DEBUGWHERE}.
*/

BOOL getdebugwhere(void)
{
  return debugwhere;
}

/*EE
  The following function sets the debug level by looking up the 
  environment variable \texttt{DEBUGLEVEL}. Moreover, the environment 
  variable \texttt{DEBUGWHERE} is read and \texttt{debugwhere} is set
  accordingly.
*/

void setdebuglevel(void)
{
  char *envstring;

  debugfileptr = stdout;
  if((envstring = getenv("DEBUGLEVEL")) != NULL)
  {
    if(!(strlen(envstring) == (size_t) 1 && 
       isdigit((Ctypeargumenttype) *envstring)))
    {
      fprintf(stderr,"environment variable DEBUGLEVEL=%s, ",envstring);
      fprintf(stderr,"it must be a digit between 0 and %lu\n",
              (Showuint) MAXDEBUGLEVEL);
      exit(EXIT_FAILURE);
    }
    if ((debuglevel = atoi(envstring)) > MAXDEBUGLEVEL)
    {
      fprintf(stderr,"environment variable DEBUGLEVEL=%s, ",envstring);
      fprintf(stderr,"it must be a digit between 0 and %lu\n",
              (Showuint) MAXDEBUGLEVEL);
      exit(EXIT_FAILURE);
    }
  }
  if((envstring = getenv("DEBUGWHERE")) != (char *) NULL)
  {
    if(strcmp(envstring,"on") == 0)
    {
      debugwhere = True;
    } else
    {
      if(strcmp(envstring,"off") == 0)
      {
        debugwhere = False;
      } else
      {
        fprintf(stderr,"environment variable DEBUGWHERE=%s, ",envstring);
        fprintf(stderr,"it must be set to \"on\" or \"off\"\n");
        exit(EXIT_FAILURE);
      }
    }
  }
  CHECKALLTYPESIZES
#ifdef WITHSYSCONF
  showmemsize();
#endif
}

/*EE
  The following function opens the given filename for writing the debug 
  messages  to. It also sets the debug level. 
  This function is called only very rarely. If only \texttt{setdebuglevel}
  is called, then the output goes to standard output.
*/

void setdebuglevelfilename(char *filename)
{
  FILEOPEN(debugfileptr,filename,"w");
  setdebuglevel();
}

/*EE
  The following function looks up the output pointer. 
*/

FILE *getdbgfp(void)
{
  if(debugfileptr == NULL)
  {
    fprintf(stderr,"DEBUGLEVELSET not called\n");
    exit(EXIT_FAILURE);
  }
  return debugfileptr;
}

/*EE
  The following function closes the debug output pointer, if it is not 
  standard out.
*/

void debugclosefile(void)
{
  if(debugfileptr == NULL)
  {
    fprintf(stderr,"cannot close debugfileptr\n");
    exit(EXIT_FAILURE);
  }
  if(fclose(debugfileptr) != 0)
  {
    NOTSUPPOSED;
  }
}

#endif  /* DEBUG */
