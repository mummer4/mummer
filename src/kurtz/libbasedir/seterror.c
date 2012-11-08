/*
  Copyright (c) 2003 by Stefan Kurtz and The Institute for
  Genomic Research.  This is OSI Certified Open Source Software.
  Please see the file LICENSE for licensing information and
  the file ACKNOWLEDGEMENTS for names of contributors to the
  code base.
*/

//\Ignore{

#include <stdio.h>
#include "types.h"

//}

/*EE
  This module implements a simple mechanism to write error messages
  into a global buffer, and to output this buffer when required. 
  We also maintain a global error code.
*/

/*
  The buffer to write the error message is of the following size.
*/

#define MAXERRORMSG 1024

static Sint errorcode = 0;
static char errormessage[MAXERRORMSG+1];

/*EE
  The following function returns the size of the buffer, and thus
  the maximal length of an error message.
*/

Sint maxerrormsg(void)
{
  return (Sint) MAXERRORMSG;
}

/*EE
  The following function returns a reference to the buffer for storing
  or retrieving the error message.
*/

char *messagespace(void)
{
  return errormessage; // write message with sprintf into this array
}

/*EE
  The following function sets the error code. It should be
  called in the innermost function in which the error occurs.
*/

void seterror(Sint code)
{
  errorcode = code;
}

/*EE
  The following function returns the error code.
*/

Sint geterror(void)
{
  return errorcode;
}

/*EE
  The following function resets the error code to 0.
*/

void reseterror(void)
{
  errorcode = 0;
}
