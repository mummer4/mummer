/* Programmer:  A. Delcher
*     Written:  7 September 1998
*        File:  ~delcher/Align/annotate.cc
*
*  This program reads the output from the  gaps  program and
*  adds alignment information to it.  Reads frag info from
*  the file named on the '>' lines of the gap file.  The word
*  "reverse" after that name will make the query sequence be reverse
*  complemented.
*/

#include <mummer/tigrinc.hh>
#include <vector>
using namespace std;

#define  FIELD_LEN  20
#define  MAX_ALIGN  10000
#define  MAX_LINE_LEN  200
#define  MAX_NAME_LEN  500
#define  PREFIX_LEN  10
#define  WIDTH  60

#ifndef  EXIT_FAILURE
  #define  EXIT_FAILURE  -1
#endif
#ifndef  EXIT_SUCCESS
  #define  EXIT_SUCCESS  0
#endif


void  Show_Alignment (char [], long int, char [], long int);


char  Line [MAX_LINE_LEN];
FILE  * Gaps_With_Errors_File;



int main  (int argc, char * argv [])
 
  {
   FILE  * Gap_File, * Data_File;
   char  * Data1, * Data2;
   char  Name [MAX_NAME_LEN], Tag [MAX_NAME_LEN];
   char  Olap_String [FIELD_LEN], Gap_String1 [FIELD_LEN], Gap_String2 [FIELD_LEN];
   int  Tag_Len;
   long int  Input_Size, Len2, Line_Len;
   //   long int Len1;
   long int  Start1, Start2, Len, Olap, Gap1, Gap2;
   long int  i;

   if  (argc < 3)
       {
        printf ("Usage: %s <gapfile> <datafile> \n", argv [0]);
        exit (EXIT_FAILURE);
       }


   Data_File = File_Open (argv [2], "r");
   Data2 = (char *) Safe_malloc (INIT_SIZE);
   Input_Size = INIT_SIZE;

   Gaps_With_Errors_File = File_Open ("witherrors.gaps", "w");

   Read_String (Data_File, Data2, Input_Size, Name, false);
   fclose (Data_File);
   Data_File = NULL;
   Len2 = strlen (Data2 + 1);

   Gap_File = File_Open (argv [1], "r");

   Data1 = (char *) Safe_malloc (INIT_SIZE);
   Input_Size = INIT_SIZE;

   while  (fgets (Line, MAX_LINE_LEN, Gap_File) != NULL)
     {
      printf ("%s", Line);
      for  (i = 0;  isspace (Line [i]) && Line [i] != '\n';  i ++)
        ;
      if  (Line [i] == '>')
          {
           Tag [0] = '\0';
           sscanf (Line + i + 1, "%s %s", Name, Tag);
	   if ( ! strcmp (Name, "Wrap") )
	     continue;

           Data_File = File_Open (Name, "r");
           Read_String (Data_File, Data1, Input_Size, Name, false);
           fclose (Data_File);

           //           Len1 = strlen (Data1 + 1);

           Tag_Len = strlen (Tag);
           for  (i = 0;  i < Tag_Len;  i ++)
             Tag [i] = tolower (Tag [i]);
           if  (strcmp (Tag, "reverse") == 0)
	     Reverse_Complement (Data2, 1, Len2);
           fprintf (Gaps_With_Errors_File, "%s", Line);
           continue;
          }
      else if  (! isdigit (Line [i]))
          {
           fprintf (Gaps_With_Errors_File, "%s", Line);
           continue;
          }
      sscanf (Line, "%ld %ld %ld %s %s %s", & Start1, & Start2,
          & Len, Olap_String, Gap_String1, Gap_String2);

      if  (Gap_String1 [0] == '-' || Gap_String2 [0] == '-')
          {
           fprintf (Gaps_With_Errors_File, "%s", Line);
           continue;
          }
      Gap1 = strtol (Gap_String1, NULL, 10);
      Gap2 = strtol (Gap_String2, NULL, 10);
      if  (isalpha (Olap_String [0]))
          Olap = 0;
        else
          Olap = strtol (Olap_String, NULL, 10);

      Start1 -= Gap1 + PREFIX_LEN - Olap;
      Start2 -= Gap2 + PREFIX_LEN - Olap;

      Line_Len = strlen (Line);
      Line [Line_Len - 1] = '\0';

      Show_Alignment (Data1 + Start1 - 1, Gap1 - Olap + 2 * PREFIX_LEN,
                      Data2 + Start2 - 1, Gap2 - Olap + 2 * PREFIX_LEN);

     }
   
   fclose (Gaps_With_Errors_File);

   return  0;
  }


template<typename T>
class Matrix_t
{
public:
  Matrix_t()
  { clear(); }

  void clear()
  {
    d_m.clear();
    nRows_m = nCols_m = 0;
  }

  void resize(long nRows, long nCols)
  {
    nRows_m = nRows;
    nCols_m = nCols;
    if ( (unsigned long)nRows * (unsigned long)nCols > d_m.size() )
      {
        try {
          d_m.resize(nRows*nCols);
        } catch (...) {
          d_m.clear();
          throw;
        }
      }
  }

  inline T & operator()(long row, long col)
  { return d_m[nCols_m*row+col]; }

private:

  vector<T> d_m;
  long nRows_m, nCols_m;
};

void  Show_Alignment (char A [], long int M, char B [], long int N)

//  Print the alignment between strings  A [1 .. M]  and  B [1 .. N] .

  {
   static Matrix_t<int> D;
   static Matrix_t<char> Op;
   static vector<char> Show_A;
   static vector<char> Show_B;
   int  Errors, Tmp;
   long int  i, j, Ct;

   if  (M >= MAX_ALIGN || N >= MAX_ALIGN)
       {
        printf ("\n   *** Too long ***\n\n");
        fprintf (Gaps_With_Errors_File, "%s %7s\n", Line, "-");
        return;
       }

   try {
     D.resize(M+1,N+1);
     Op.resize(M+1,N+1);
     Show_A.resize(M+N+2);
     Show_B.resize(M+N+2);
   } catch (...) {
     printf ("\n   *** Too long ***\n\n");
     fprintf (Gaps_With_Errors_File, "%s %7s\n", Line, "-");
     return;
   }

   D (0,0) = 0;
   Op (0,0) = 'a';
   for  (i = 1;  i <= M;  i ++)
     {
      D (i,0) = i + 1;
      Op (i,0) = 'i';
     }
   for  (j = 1;  j <= N;  j ++)
     {
      D (0,j) = j + 1;
      Op (0,j) = 'd';
     }

   for  (i = 1;  i <= M;  i ++) // for each row (A)
     for  (j = 1;  j <= N;  j ++) // for each col (B)
       {
        D (i,j) = D (i,j - 1) + 1 + (Op (i,j - 1) == 'd' ? 0 : 1);
        Op (i,j) = 'd';
        Tmp = D (i - 1,j) + 1 + (Op (i - 1,j) == 'i' ? 0 : 1);
        if  (Tmp < D (i,j))
            {
             D (i,j) = Tmp;
             Op (i,j) = 'i';
            }
        Tmp = D (i - 1,j - 1) + (A [i] == B [j] ? 0 : 1);
        if  (Tmp < D (i,j))
            {
             D (i,j) = Tmp;
             Op (i,j) = 'a';
            }
       }

   Ct = 0;
   i = M;
   j = N;
   while  (i > 0 || j > 0)
     switch  (Op (i,j))
       {
        case  'a' :                      // align
          Show_A [Ct] = A [i --];
          Show_B [Ct ++] = B [j --];
          break;
        case  'i' :                      // insert into A
          Show_A [Ct] = A [i --];
          Show_B [Ct ++] = '.';
          break;
        case  'd' :                      // delete from A
          Show_A [Ct] = '.';
          Show_B [Ct ++] = B [j --];
          break;
        default :
          printf (" *** i = %ld   j = %ld   Op (i,j) = %c\n",
                      i, j, Op (i,j));
          exit (EXIT_FAILURE);
       }
   Errors = 0;
   for  (i = 0;  i < Ct;  i ++)
     if  (Show_A [i] != Show_B [i])
         Errors ++;
   printf ("    Errors = %d\n", Errors);
   fprintf (Gaps_With_Errors_File, "%s %7d\n", Line, Errors);
   do
     {
      printf ("T:  ");
      for  (i = Ct - 1;  i >= 0 && i >= Ct - WIDTH;  i --)
        putchar (Show_A [i]);
      putchar ('\n');
      printf ("S:  ");
      for  (i = Ct - 1;  i >= 0 && i >= Ct - WIDTH;  i --)
        putchar (Show_B [i]);
      putchar ('\n');
      printf ("    ");
      for  (i = Ct - 1;  i >= 0 && i >= Ct - WIDTH;  i --)
        putchar (Show_A [i] == Show_B [i] ? ' ' : '^');
      putchar ('\n');
      Ct -= WIDTH;
     }  while  (Ct > 0);
   return;
  }
