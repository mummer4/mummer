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
#include <sys/time.h>
#include <sys/resource.h>
#include <unistd.h>
#include <string.h>
#include "types.h"
#include "errordef.h"
#include "megabytes.h"
#include "debugdef.h"

//}

/*EE
  This file contains functions to store pointers to dynamically allocated
  spaceblocks, and to maintain the number of cells and their size in
  each block. The arguments \texttt{file} and \texttt{line} (if they occur)
  are always the filename and the linenumber the function is called from.
  To supply these arguments, we recommend to call the corresponding 
  functions via some useful macros, as defined in the file 
  \texttt{spacedef.h}. 
  \begin{enumerate}
  \item
  The function \texttt{allocspaceviaptr} should be called
  via the macro \texttt{ALLOCSPACE}.
  \item
  The function \texttt{freespaceviaptr} should be called
  via the macro \texttt{FREESPACE}.
  \item
  The function \texttt{dynamicstrdup} should be called
  via the macro \texttt{DYNAMICSTRDUP}.
  \end{enumerate}
*/

/*EE
  The following is a general macro to print fatal error messages.
*/

#define ALLOCVIAFATAL(M) \
        fprintf(stderr,"file \"%s\", line %lu: "\
                       " allocandusespaceviaptr(%lu,%lu) failed:%s\n",\
                       file,(Showuint) line,(Showuint) size,\
                       (Showuint) number,M);\
        exit(EXIT_FAILURE)

typedef struct
{
  void *spaceptr;      // ptr to the spaceblock
  Uint sizeofcells,    // size of cells of the block
       numberofcells;  // number of cells in the block
  char *fileallocated; // the filenames where the block was allocated
  Uint lineallocated;  // the linenumber where the
} Blockdescription;

/*@null@*/ static Blockdescription *blocks = NULL;

static Uint numberofblocks = 0, // numberofblocks
            nextfreeblock = 0,  // index of next free block
            currentspace = 0,   // currently allocated num of bytes
            spacepeak = 0;      // maximally allocated num of bytes

/*
  The following two tables store important information to
  generate meaningfull error messages.
*/

/*
  The following function sets the soft limit on the data size to the hard
  limit, if it is not already identical to this.
  This is necessary for the DEC alpha platform, where the soft limit
  is too small. \texttt{setmaxspace} is called for the first time, 
  space is allocated.
*/

static void setmaxspace(void)
{
  Getrlimitreturntype rc;
  struct rlimit rls;

  /*@ignore@*/
  if((rc = getrlimit(RLIMIT_DATA,&rls)) != 0)
  /*@end@*/
  {
    fprintf(stderr,"cannot find rlimit[RLIMIT_DATA]\n");
    exit(EXIT_FAILURE);
  }

  if(rls.rlim_cur < rls.rlim_max)
  {
    rls.rlim_cur = rls.rlim_max;
    /*@ignore@*/
    if((rc = setrlimit(RLIMIT_DATA, &rls)) != 0)
    /*@end@*/
    {
      fprintf(stderr,"cannot set rlimit[RLIMIT_DATA]\n");
      exit(EXIT_FAILURE);
    }
  }
}

/*
  The following two functions \texttt{addspace} and \texttt{subtractspace} 
  maintain the variables \texttt{currentspace} and \texttt{spacepeak}.
*/

static void addspace(Uint space)
{
  if(currentspace == 0)
  {
    setmaxspace();
    currentspace = space;
  } else
  {
    currentspace += space;
  }
  if(currentspace > spacepeak)
  {
    spacepeak = currentspace;
    DEBUG1(2,"# spacepeak = %.2f reached\n",MEGABYTES(spacepeak));
  }
}

static void subtractspace(Uint space)
{
  currentspace -= space;
}

/*EE
  The following function allocates \texttt{number} cells of \texttt{size}
  for a given pointer \texttt{ptr}. If this is \texttt{NULL}, then the next 
  free block is used. Otherwise, we look for the block number corresponding 
  to \texttt{ptr}. If there is none, then the program exits with exit code 1. 
*/

/*@notnull@*/ void *allocandusespaceviaptr(char *file,Uint line, 
                                           /*@null@*/ void *ptr,
                                           Uint size,Uint number)
{
  Uint i, blocknum;

  DEBUG2(2,"\n# allocandusespaceviaptr(file=%s,line=%lu)\n",file,
                      (Showuint) line);
  if(nextfreeblock > 0)
  {
    NOTSUPPOSEDTOBENULL(blocks);
    for(blocknum=0; blocknum < nextfreeblock; blocknum++)
    {
      if(blocks[blocknum].spaceptr == ptr)
      {
        break;
      }
    }
  } else
  {
    blocknum = 0;
  }
  if(blocknum == nextfreeblock)
  {
    if(ptr == NULL)
    {
      nextfreeblock += 64;
      blocks = (Blockdescription *) realloc(blocks,
                                            (size_t) (sizeof(Blockdescription)* 
                                                      nextfreeblock));
      if(blocks == NULL)
      {
        ALLOCVIAFATAL("not enough space for the block descriptions available");
      }
      for(i=blocknum; i < nextfreeblock; i++)
      {
        blocks[i].spaceptr = NULL;
        blocks[i].sizeofcells = 0;
        blocks[i].numberofcells = 0;
      }
    } else
    {
      ALLOCVIAFATAL("cannot find space block");
    }
  }
  NOTSUPPOSEDTOBENULL(blocks);
  subtractspace(blocks[blocknum].numberofcells * blocks[blocknum].sizeofcells);
  addspace(size*number);
  DEBUG3(2,"# allocandusespaceviaptr:block %lu: %lu cells of size %lu\n",
            (Showuint) blocknum,(Showuint) number,(Showuint) size);
  DEBUG2(2,"# previously allocated for this block (%lu,%lu)\n",
            (Showuint) blocks[blocknum].numberofcells,
            (Showuint) blocks[blocknum].sizeofcells);
  blocks[blocknum].numberofcells = number;
  blocks[blocknum].sizeofcells = size;
  blocks[blocknum].fileallocated = file;
  blocks[blocknum].lineallocated = line;
  if(blocks[blocknum].spaceptr == NULL)
  {
    numberofblocks++;
  }
  if((blocks[blocknum].spaceptr 
      = realloc(blocks[blocknum].spaceptr,(size_t) (size*number))) == NULL)
  {
    ALLOCVIAFATAL("not enough memory");
  }
  DEBUG0(2,"# allocandusespaceviaptr Okay\n");
  NOTSUPPOSEDTOBENULL(blocks[blocknum].spaceptr);
  return blocks[blocknum].spaceptr;
}

/*EE
  The following function makes a copy of a 0-terminated string pointed to by 
  \texttt{source}. 
*/

/*@notnull@*/ char *dynamicstrdup(char *file,Uint line,char *source)
{
  Uint len;
  char *dest;

  NOTSUPPOSEDTOBENULL(source);
  len = (Uint) strlen(source);
  dest = (char *) allocandusespaceviaptr(file,line,NULL,(Uint) sizeof(char),
                                         len+1);
  strcpy(dest,source);
  return dest;
}

/*EE
  The following function frees the space for the given pointer 
  \texttt{ptr}. This cannot be \texttt{NULL}. 
*/

void freespaceviaptr(char *file,Uint line,void *ptr)
{
  Uint blocknum;

  DEBUG2(2,"\n# freespaceviaptr(file=%s,line=%lu):\n",file,(Showuint) line);
  if(ptr == NULL)
  {
    fprintf(stderr,"freespaceviaptr(file=%s,line=%lu): Cannot free NULL-ptr\n",
                    file,(Showuint) line);
    exit(EXIT_SUCCESS);
  }
  NOTSUPPOSEDTOBENULL(blocks);
  for(blocknum=0; blocknum < nextfreeblock; blocknum++)
  {
    if(blocks[blocknum].spaceptr == ptr)
    {
      break;
    }
  }
  if(blocknum == nextfreeblock)
  {
    fprintf(stderr,"freespaceviaptr(file=%s,line=%lu): "
                   " cannot find space block\n",
            file,(Showuint) line);
    exit(EXIT_FAILURE);
  }
  free(blocks[blocknum].spaceptr);
  subtractspace(blocks[blocknum].numberofcells * blocks[blocknum].sizeofcells);
  DEBUG3(2,"# freespaceviaptr:block %lu: %lu cells of size %lu\n",
            (Showuint) blocknum,
            (Showuint) blocks[blocknum].numberofcells,
            (Showuint) blocks[blocknum].sizeofcells);
  DEBUG2(2,"# this block was allocated in file \"%s\", line %lu\n",
            blocks[blocknum].fileallocated,
            (Showuint) blocks[blocknum].lineallocated);
  blocks[blocknum].numberofcells = 0;
  blocks[blocknum].sizeofcells = 0;
  blocks[blocknum].fileallocated = NULL;
  blocks[blocknum].lineallocated = 0;
  blocks[blocknum].spaceptr = NULL;
  if(numberofblocks == 0)
  {
    NOTSUPPOSED;
  }
  numberofblocks--;
}

//\IgnoreLatex{

/*EE
  The following function frees the space for all main memory blocks
  which have not already been freed.
*/

void wrapspace(void)
{
  Uint blocknum;

  DEBUG0(2,"# wrapspace\n");
  NOTSUPPOSEDTOBENULL(blocks);
  for(blocknum=0; blocknum < nextfreeblock; blocknum++)
  {
    if(blocks[blocknum].spaceptr != NULL)
    {
      DEBUG3(2,"# free block %lu: %lu cells of size %lu\n",
                (Showuint) blocknum,
                (Showuint) blocks[blocknum].numberofcells,
                (Showuint) blocks[blocknum].sizeofcells);
      free(blocks[blocknum].spaceptr);
      blocks[blocknum].spaceptr = NULL;
    }
    subtractspace(blocks[blocknum].sizeofcells * 
                  blocks[blocknum].numberofcells);
    blocks[blocknum].sizeofcells = 0;
    blocks[blocknum].numberofcells = 0;
    blocks[blocknum].fileallocated = NULL;
    blocks[blocknum].lineallocated = 0;
  }
}

//}

/*EE
  The following function prints a list of block numbers 
  which have not been freed. For each block number the filename
  and line number in which the call appears allocating which 
  allocated this block.
*/

void activeblocks(void)
{
  Uint blocknum;

  NOTSUPPOSEDTOBENULL(blocks);
  for(blocknum=0; blocknum < nextfreeblock; blocknum++)
  {
    if(blocks[blocknum].spaceptr != NULL)
    {
      fprintf(stderr,"# active block %lu: ",(Showuint) blocknum);
      fprintf(stderr,"allocated in file \"%s\", line %lu\n",
              blocks[blocknum].fileallocated,
              (Showuint) blocks[blocknum].lineallocated);
    }
  }
}

/*EE
  The following function checks if all blocks previously allocated, have 
  explicitely been freed. If there is a block that was not freed, then 
  an error is reported accordingly. We recommend to call this function
  before the program terminates. This easily allows to discover 
  space leaks.
*/

void checkspaceleak(void)
{
  Uint blocknum;

  NOTSUPPOSEDTOBENULL(blocks);
  for(blocknum=0; blocknum < nextfreeblock; blocknum++)
  {
    if(blocks[blocknum].spaceptr != NULL)
    {
      fprintf(stderr,"space leak: main memory for block %lu not freed\n",
              (Showuint) blocknum);
      fprintf(stderr,"%lu cells of size %lu\n",
              (Showuint) blocks[blocknum].numberofcells,
              (Showuint) blocks[blocknum].sizeofcells);
      fprintf(stderr,"allocated: ");
      if(blocks[blocknum].fileallocated == NULL)
      {
        fprintf(stderr,"cannot identify\n");
      } else
      {
        fprintf(stderr,"file \"%s\", line %lu\n",
               blocks[blocknum].fileallocated,
               (Showuint) blocks[blocknum].lineallocated);
      }
      exit(EXIT_FAILURE);
    }
  }
  if(numberofblocks > 0)
  {
    fprintf(stderr,"space leak: number of blocks = %u\n",numberofblocks);
    exit(EXIT_FAILURE);
  } 
  free(blocks);
  blocks = NULL;
  numberofblocks = 0;
  nextfreeblock = 0;
  currentspace = 0;
  spacepeak = 0;
}

/*EE
  The following function shows the space peak in megabytes on \texttt{stderr}.
*/

void showspace(void)
{
  fprintf(stderr,"# space peak in megabytes: %.2f\n",MEGABYTES(spacepeak));
}

/*EE
  The following function returns the space peak in bytes.
*/

Uint getspacepeak(void)
{
  return spacepeak;
}

/*EE
  The following function delivers the space limit of the machine 
  in megabytes. This only works if the variable
  \texttt{WITHSYSCONF} is defined. This is currently the case for
  Linux and Solaris.
*/

#ifdef WITHSYSCONF
#ifdef DEBUG
void showmemsize(void)
{
  Sint pagesize = (Sint) sysconf((Sysconfargtype) _SC_PAGESIZE);
  Sint physpages = (Sint) sysconf((Sysconfargtype) _SC_PHYS_PAGES);

  DEBUG1(1,"# pagesize = %ld\n",(Showsint) pagesize);
  DEBUG1(1,"# number of physical pages = %ld\n",(Showsint) physpages);
  DEBUG1(1,"# memory size = %.0f MB\n",MEGABYTES(pagesize * physpages));
}
#endif
#endif
