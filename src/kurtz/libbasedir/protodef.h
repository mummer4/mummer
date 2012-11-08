/*
  Copyright (c) 2003 by Stefan Kurtz and The Institute for
  Genomic Research.  This is OSI Certified Open Source Software.
  Please see the file LICENSE for licensing information and
  the file ACKNOWLEDGEMENTS for names of contributors to the
  code base.
*/

#ifndef PROTODEF_H
#define PROTODEF_H

#include <sys/types.h>
#include <stdio.h>
#include <stdlib.h>
#include "types.h"
#include "optdesc.h"
#include "multidef.h"
#include "mumcand.h"
#ifdef __cplusplus
extern "C" {
#endif
Sint mumuniqueinquery(void *processinfo,
                      Sint (*processmum)(void *,Uint,Uint,Uint,Uint),
                      ArrayMUMcandidate *mumcand);
void initclock(void);
double getruntime(void);
Uint getclockticks(void);
Sint getdebuglevel(void);
BOOL getdebugwhere(void);
void setdebuglevel(void);
void setdebuglevelfilename(char *filename);
FILE *getdbgfp(void);
void debugclosefile(void);
Sint simplefileOpen(char *filename,Uint *numofbytes);
/*@null@*/ void *creatememorymapforfiledesc(char *file,Uint line,Sint fd,
                                            BOOL writemap,Uint numofbytes);
/*@null@*/ void *creatememorymap(char *file,Uint line,char *filename,
                                 BOOL writemap,Uint *numofbytes);
Sint deletememorymap(char *file,Uint line,void *mappedfile);
void mmcheckspaceleak(void);
Sint mmwrapspace(void);
void mmshowspace(void);
Uint mmgetspacepeak(void);
void initmultiseq(Multiseq *multiseq);
void freemultiseq(Multiseq *multiseq);
Sint overallsequences(BOOL rcmode,Multiseq *multiseq,void *applyinfo,
                      Sint(*apply)(void *,Uint,Uchar *,Uint));
Sint getrecordnum(Uint *recordseps,Uint numofrecords,Uint totalwidth,
                  Uint position);
Sint getseqnum(Multiseq *multiseq,Uint position);
Sint pos2pospair(Multiseq *multiseq,PairUint *pos,Uint position);
void initoptions(OptionDescription *options,Uint numofoptions);
Sint addoption(OptionDescription *options,Uint numofoptions,
               Uint optnum,char *optname,char *optdesc);
Sint procoption(OptionDescription *opt,Uint numofopt,char *optstring);
void showoptions(FILE *outfp,char *program,OptionDescription *opt,
                 Uint numofopt);
void showoptionswithoutexclude(FILE *outfp,char *program,
                               OptionDescription *opt,
                               Sint *excludetab,Uint numofopt);
Sint checkdoubleexclude(Uint numofopts,OptionDescription *opt,
                        Sint *excludetab,Uint len);
Sint checkexclude(OptionDescription *opt,Sint *excludetab,Uint len);
void showexclude(OptionDescription *opt,Sint *excludetab,Uint len);
Sint safestringcopy(char *dest,char *source,Sint maxlen);
Sint maxerrormsg(void);
char *messagespace(void);
void seterror(Sint code);
Sint geterror(void);
void reseterror(void);
/*@notnull@*/ void *allocandusespaceviaptr(char *file,Uint line,
                                           /*@null@*/ void *ptr,
                                           Uint size,Uint number);
/*@notnull@*/ char *dynamicstrdup(char *file,Uint line,char *source);
void freespaceviaptr(char *file,Uint line,void *ptr);
void wrapspace(void);
void activeblocks(void);
void checkspaceleak(void);
void showspace(void);
Uint getspacepeak(void);
void showmemsize(void);
#ifdef __cplusplus
}
#endif
#endif
