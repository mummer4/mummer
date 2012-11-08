/*
  Copyright (c) 2003 by Stefan Kurtz and The Institute for
  Genomic Research.  This is OSI Certified Open Source Software.
  Please see the file LICENSE for licensing information and
  the file ACKNOWLEDGEMENTS for names of contributors to the
  code base.
*/

//\Ignore{

#include <stdio.h>
#include <time.h>
#include "types.h"

//}

/*
  This module implements function to measures the running time
  of programs or program parts.
*/

/*
  The following values store the the clockticks at start time 
  and stop time of the clock.
*/

static clock_t startclock, 
               stopclock;

/*EE
  The following function initializes the clock.
*/

void initclock(void)
{ 
  startclock = clock(); 
}

/*EE
  The following function delivers the time since the 
  clock was initialized. The time is reported in seconds
  as a floating point value.
*/

double getruntime(void)
{
   stopclock = clock();
   return (stopclock-startclock) / (double) CLOCKS_PER_SEC;
}

/*EE
  The following function delivers the clock ticks betwenn 
  \texttt{startclock} to \texttt{stopclock}.
*/

Uint getclockticks(void)
{
   stopclock = clock();
   return (Uint) (stopclock-startclock);
}
