/*
  Copyright (c) 2003 by Stefan Kurtz and The Institute for
  Genomic Research.  This is OSI Certified Open Source Software.
  Please see the file LICENSE for licensing information and
  the file ACKNOWLEDGEMENTS for names of contributors to the
  code base.
*/

#include <stdio.h>
#include <stdlib.h>
#include <sys/types.h>
#include "debugdef.h"
#include "args.h"
#include "types.h"
#include "streedef.h"
#include "protodef.h"
#include "streeacc.h"
#include "errordef.h"
#include "spacedef.h"

MAINFUNCTION
{
  Uchar *text;
  Uint textlen;
  Suffixtree stree;

  DEBUGLEVELSET;

  CHECKARGNUM(2,"filename");
  text = (Uchar *) CREATEMEMORYMAP(argv[1],False,&textlen);
  if(text == NULL)
  {
    STANDARDMESSAGE;
  }
  CONSTRUCTSTREE(&stree,text,textlen,return EXIT_FAILURE);
#ifdef DEBUG
  enumlocations(&stree,checklocation);
#endif
  freestree(&stree);
  if(DELETEMEMORYMAP(text) != 0)
  {
    STANDARDMESSAGE;
  }
  return EXIT_SUCCESS;
}
