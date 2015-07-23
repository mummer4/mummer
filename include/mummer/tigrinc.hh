#ifndef  __TIGRINC_HH
#define  __TIGRINC_HH


#include  <cstdio>
#include  <cstdlib>
#include  <cmath>
#include  <cstring>
#include  <cctype>
#include  <climits>
#include  <cfloat>
#include  <ctime>
#include  <cassert>
#include  <cerrno>
#include  <unistd.h>


#ifndef  EXIT_FAILURE
  #define  EXIT_FAILURE  -1
#endif
#ifndef  EXIT_SUCCESS
  #define  EXIT_SUCCESS  0
#endif

#define  INCR_SIZE 10000
#define  SMALL_INIT_SIZE 100
#define  INIT_SIZE 10000
#define  MAX_LINE 1024


FILE *  File_Open  (const char *, const char *);
void *  Safe_calloc  (size_t, size_t);
void *  Safe_malloc  (size_t);
void *  Safe_realloc  (void *, size_t);
char  Complement  (char);
bool CompareIUPAC (char, char);
bool  Read_String  (FILE *, char * &, long int &, char [], bool);
void  Reverse_Complement (char S [], long int Lo, long int Hi);

#endif
