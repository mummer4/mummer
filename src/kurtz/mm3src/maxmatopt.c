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
#include "types.h"
#include "optdesc.h"
#include "debugdef.h"
#include "errordef.h"
#include "protodef.h"
#include "maxmatdef.h"

//}

/*EE
  This file contains functions to parse the possible 
  options of \texttt{maxmat3} and to appropriately initialize 
  the \texttt{mmcallinfo}-record according to the given options.
*/

//\IgnoreLatex{

/*
  The default value for the minimal unique match length.
*/

#define DEFAULTMINUNIQUEMATCHLEN 20

//}

/*EE
  The following type declares symbolic constants for the options.
*/

typedef enum
{
  OPTMUM = 0,
  OPTMUMCAND,
  OPTMUMREF,
  OPTMAXMATCH,
  OPTMATCHNUCLEOTIDESONLY,
  OPTLEASTLENGTH,
  OPTCOMPUTEBOTHDIRECTIONS,
  OPTONLYREVERSECOMPLEMENT,
  OPTSHOWSTRING,
  OPTSHOWREVERSEPOSITIONS,
  OPTFOURCOLUMN,
  OPTSHOWSEQUENCELENGTHS,
  OPTH,
  OPTHELP,
  NUMOFOPTIONS
} Optionnumber;

/*
  The following function stores the help-text for the option \texttt{-l}.
  This is necessary, since the text depends on the value of the
  symbolic constant \texttt{DEFAULTMINUNIQUEMATCHLEN};
*/

static void makeleastlengthtext(char *spacefortext)
{
  sprintf(spacefortext,"set the minimum length of a match\n"
                       "if not set, the default value is %lu",
                       (Showuint) DEFAULTMINUNIQUEMATCHLEN);
}

/*
  The following function shows a usage line including the
  possible options.
*/

static void showusage(char *program,OptionDescription *options,
                      Uint numofoptions)
{
  printf("Usage: %s [options] <reference-file> <query-files>\n\n"
         "Find and output (to stdout) the positions and length of all\n"
         "sufficiently long maximal matches of a substring in\n"
         "<query-file> and <reference-file>\n\n",program);
  printf("Options:\n");
  showoptions(stdout,program,options,numofoptions);
}

/*EE
  The following function declares the possible options
  in a record \texttt{options}. It then ananlyzes the \texttt{argv}-vector
  step by step. If everything is okay, 0 is returned and the 
  \texttt{mmcallinfo} is correctly initialized.
  Otherwise, a negative value is returned.
*/

Sint parsemaxmatoptions(MMcallinfo *mmcallinfo,Argctype argc, char **argv)
{
  OptionDescription options[NUMOFOPTIONS];   // store the options
  Sint optval;         // neg. return val. if error, otherwise option number
  Uint argnum;         // pointer to argv
  Scaninteger readint; // temporary integer to read value from string
  char leastlengthtext[128+1];

  DEBUGLEVELSET;
  initoptions(&options[0],(Uint) NUMOFOPTIONS);
  ADDOPTION(OPTMUM,"-mum",
            "compute maximal matches that are unique in both sequences");
  ADDOPTION(OPTMUMREF,"-mumreference",
	    "compute maximal matches that are unique in the reference-\n"
            "sequence but not necessarily in the query-sequence (default)");
  ADDOPTION(OPTMUMCAND,"-mumcand",
            "same as -mumreference");
  ADDOPTION(OPTMAXMATCH,"-maxmatch",
	    "compute all maximal matches regardless of their uniqueness");
  ADDOPTION(OPTMATCHNUCLEOTIDESONLY,"-n",
            "match only the characters a, c, g, or t\n"
            "they can be in upper or in lower case");
  makeleastlengthtext(&leastlengthtext[0]);
  ADDOPTION(OPTLEASTLENGTH,"-l",&leastlengthtext[0]);
  ADDOPTION(OPTCOMPUTEBOTHDIRECTIONS,"-b",
            "compute forward and reverse complement matches");
  ADDOPTION(OPTONLYREVERSECOMPLEMENT,"-r",
            "only compute reverse complement matches");
  ADDOPTION(OPTSHOWSTRING,"-s",
            "show the matching substrings");
  ADDOPTION(OPTSHOWREVERSEPOSITIONS,"-c",
            "report the query-position of a reverse complement match\n"
            "relative to the original query sequence");
  ADDOPTION(OPTFOURCOLUMN,"-F",
	    "force 4 column output format regardless of the number of\n"
	    "reference sequence inputs");
  ADDOPTION(OPTSHOWSEQUENCELENGTHS,"-L",
            "show the length of the query sequences on the header line");
  ADDOPTION(OPTH,"-h",
	    "show possible options");
  ADDOPTION(OPTHELP,"-help",
            "show possible options");
  mmcallinfo->showstring = False;
  mmcallinfo->reversecomplement = False;
  mmcallinfo->forward = True;
  mmcallinfo->showreversepositions = False;
  mmcallinfo->fourcolumn = False;
  mmcallinfo->showsequencelengths = False;
  mmcallinfo->matchnucleotidesonly = False;
  mmcallinfo->cmum = False;
  mmcallinfo->cmumcand = False;
  mmcallinfo->cmaxmatch = False;
  mmcallinfo->minmatchlength = (Uint) DEFAULTMINUNIQUEMATCHLEN;

  if(argc == 1)
  {
    showusage(argv[0],&options[0],(Uint) NUMOFOPTIONS);
    return 1;
  }

  for(argnum = UintConst(1); argnum < (Uint) argc && argv[argnum][0] == '-'; 
      argnum++)
  {
    optval = procoption(options,(Uint) NUMOFOPTIONS,argv[argnum]);
    if(optval < 0)
    {
      return -1;
    }
    switch(optval)
    {
      case OPTSHOWSTRING:
        mmcallinfo->showstring = True; 
        break;
      case OPTCOMPUTEBOTHDIRECTIONS:
        mmcallinfo->reversecomplement = True; 
        break;
      case OPTSHOWREVERSEPOSITIONS:
        mmcallinfo->showreversepositions = True; 
        break;
      case OPTLEASTLENGTH:  // additionally check the length parameter
        argnum++;
        if(argnum > (Uint) (argc-2))
        {
          ERROR1("missing argument for option %s",
                  options[OPTLEASTLENGTH].optname);
          return -2;
        }
        if(sscanf(argv[argnum],"%ld",&readint) != 1 || readint <= 0)
        {
          ERROR2("argument %s for option %s is not a positive integer",
                  argv[argnum],options[OPTLEASTLENGTH].optname);
          return -3;
        }
        mmcallinfo->minmatchlength = (Uint) readint;
        break;
      case OPTFOURCOLUMN:
	mmcallinfo->fourcolumn = True;
	break;
      case OPTSHOWSEQUENCELENGTHS:
        mmcallinfo->showsequencelengths = True; 
        break;
      case OPTMATCHNUCLEOTIDESONLY:
        mmcallinfo->matchnucleotidesonly = True; 
        break;
      case OPTONLYREVERSECOMPLEMENT:
        mmcallinfo->forward = False; 
        mmcallinfo->reversecomplement = True; 
        break;
      case OPTMAXMATCH:
	mmcallinfo->cmaxmatch = True;
	break;
      case OPTMUMREF:
      case OPTMUMCAND:
        mmcallinfo->cmumcand = True;
        break;
      case OPTMUM:
        mmcallinfo->cmum = True;
        break;
      case OPTH:
      case OPTHELP:
        showusage(argv[0],&options[0],(Uint) NUMOFOPTIONS);
        return 1;
    }
  }
  if(argnum > (Uint) (argc-2))
  {
    ERROR0("missing file arguments");
    return -4;
  }
  if(safestringcopy(&mmcallinfo->program[0],argv[0],PATH_MAX) != 0)
  {
    return -5;
  }
  if(safestringcopy(&mmcallinfo->subjectfile[0],argv[argnum],PATH_MAX) != 0)
  {
    return -6;
  }
  for(argnum++, mmcallinfo->numofqueryfiles = 0; 
      argnum < (Uint) argc; mmcallinfo->numofqueryfiles++, argnum++)
  {
    if(mmcallinfo->numofqueryfiles >= (Uint) MAXNUMOFQUERYFILES)
    {
      ERROR1("too many query files, maximal number is %lu",
              (Showuint) MAXNUMOFQUERYFILES);
      return -7;
    }
    if(safestringcopy(&mmcallinfo->queryfilelist
                       [mmcallinfo->numofqueryfiles][0],
                      argv[argnum],PATH_MAX) != 0)
    {
      return -8;
    }
  }
  /*
    verify that mum options are not interchanged
  */
  OPTIONEXCLUDE(OPTMUM,OPTMUMCAND);
  OPTIONEXCLUDE(OPTMUM,OPTMUMREF);
  OPTIONEXCLUDE(OPTMUM,OPTMAXMATCH);
  OPTIONEXCLUDE(OPTMUMCAND,OPTMAXMATCH);
  OPTIONEXCLUDE(OPTMUMREF,OPTMAXMATCH);
  if ( mmcallinfo->cmaxmatch )
    {
      mmcallinfo->cmum = False;
      mmcallinfo->cmumcand = False;
    }
  else if ( mmcallinfo->cmum )
    {

    }
  else /* default to cmumcand */
    {
      mmcallinfo->cmumcand = True;
    }
  /*
    verify that the options -b and -r are not used at the same time
  */
  OPTIONEXCLUDE(OPTCOMPUTEBOTHDIRECTIONS,OPTONLYREVERSECOMPLEMENT);
  /*
    verify that -c is only used in combination with either -b or -r
  */
  OPTIONIMPLYEITHER2(OPTSHOWREVERSEPOSITIONS,
                     OPTCOMPUTEBOTHDIRECTIONS,OPTONLYREVERSECOMPLEMENT);
  return 0;
}
