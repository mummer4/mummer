/* Programmer:  A. Delcher
*        File:  ~delcher/TIGR/gaps.cc
*
*  This program reads a list of unique matches between 2 strings and
*  outputs first the longest consistent set of matches, and
*  then all the other matches.
*/

#include <mummer/tigrinc.hh>

const bool  ALLOW_WRAPAROUND = false;


struct  Match
  {
   long int  Start1, Start2, Len;
   long int  Simple_Score, Wrap_Score;
   long int  Simple_From, Wrap_From;
   long int  Simple_Adj, Wrap_Adj;
   int  Good : 1;
   int  Wrap_Here : 1;
  };


int  Cmp  (const void * A, const void * B)

//  Return how  A  and  B  compare if converted to  Match

  {
   if  (((Match *) A) -> Start2 < ((Match *) B) -> Start2)
       return  -1;
   else if  (((Match *) A) -> Start2 == ((Match *) B) -> Start2)
       return  0;
     else
       return  1;
  }


inline
long int  Max  (long int A, long int B)

//  Return the larger of  A  and  B .

  {
   if  (A < B)
       return  B;
     else
       return  A;
  }


int main  (int argc, char * argv [])

  {
   Match  * A = NULL;
   char  * File_Name;
   int  Used_Wrap, Use_Reverse_Complement;
   long int  N = 0;
   long int  Adj, Olap, Olap1, Olap2, Prev, S1, S2, Len;
   long int  i, j, Best, Next;


   if  (argc < 2)
       {
        fprintf (stderr, "ERROR:  Specify first genome filename\n");
        exit (-1);
       }
   File_Name = argv [1];

   Use_Reverse_Complement = (argc > 2 && strcmp (argv [2], "-r") == 0);

   while  (scanf ("%ld %ld %ld", & S1, & S2, & Len) != EOF)
     {
      A = (Match *) Safe_realloc (A, (N + 1) * sizeof (Match));
      A [N] . Start1 = S1;
      A [N] . Start2 = S2;
      A [N] . Len = Len;
      A [N] . Good = false;
      A [N] . Wrap_Here = false;
      N ++;
     }

   qsort (A, N, sizeof (Match), Cmp);

   for  (i = 0;  i < N;  i ++)
     {
      A [i] . Simple_Score = A [i] . Wrap_Score = A [i] . Len;
      A [i] . Simple_Adj = A [i] . Wrap_Adj = 0;
      A [i] . Simple_From = A [i] . Wrap_From = -1;
      for  (j = 0;  j < i;  j ++)
        {
         Olap2 = A [j] . Start2 + A [j] . Len - A [i] . Start2;
         Olap = Max (0, Olap2);
         if  (A [j] . Simple_Score + A [i] . Len - Olap > A [i] . Wrap_Score)
             {
              A [i] . Wrap_From = j;
              A [i] . Wrap_Score = A [j] . Simple_Score + A [i] . Len - Olap;
              A [i] . Wrap_Adj = Olap;
              A [i] . Wrap_Here = true;
             }
         Olap1 = A [j] . Start1 + A [j] . Len - A [i] . Start1;
         Olap = Max (Olap, Olap1);
         if  (A [j] . Simple_Score + A [i] . Len - Olap > A [i] . Simple_Score)
             {
              A [i] . Simple_From = j;
              A [i] . Simple_Score = A [j] . Simple_Score + A [i] . Len - Olap;
              A [i] . Simple_Adj = Olap;
             }
         if  (A [j] . Wrap_Score + A [i] . Len - Olap >= A [i] . Wrap_Score)
             {
              A [i] . Wrap_From = j;
              A [i] . Wrap_Score = A [j] . Wrap_Score + A [i] . Len - Olap;
              A [i] . Wrap_Adj = Olap;
              A [i] . Wrap_Here = false;
             }
        }
     }

   if  (ALLOW_WRAPAROUND)
       {
        Best = 0;
        for  (i = 1;  i < N;  i ++)
          if  (A [i] . Wrap_Score > A [Best] . Wrap_Score)
              Best = i;
        Used_Wrap = false;
        for  (i = Best;  i >= 0;  i = Next)
          {
           A [i] . Good = true;
           if  (Used_Wrap)
               {
                Next = A [i] . Simple_From;
                A [i] . Wrap_Here = false;
               }
             else
               Next = A [i] . Wrap_From;
           if  (A [i] . Wrap_Here)
               Used_Wrap = true;
          }
       }
     else
       {
        Best = 0;
        for  (i = 1;  i < N;  i ++)
          if  (A [i] . Simple_Score > A [Best] . Simple_Score)
              Best = i;
        for  (i = Best;  i >= 0;  i = A [i] . Simple_From)
          A [i] . Good = true;
       }

   printf ("> %s %s Consistent matches\n", File_Name,
                  Use_Reverse_Complement ? "reverse" : "");
   Prev = -1;
   for  (i = 0;  i < N;  i ++)
     if  (A [i] . Good)
         {
          if  (Prev == -1)
              printf ("%8ld %8ld %6ld %7s %6s %6s\n",
                  A [i] . Start1, A [i] . Start2, A [i] . Len,
                  "none", "-", "-");
          else if  (ALLOW_WRAPAROUND)
              {
               Adj = A [i] . Wrap_Adj;
               if  (A [i] . Wrap_Here)
                   printf ("> Wrap around\n");
               
               printf ("%8ld %8ld %6ld",
                   A [i] . Start1 + Adj, A [i] . Start2 + Adj,
                   A [i] . Len - Adj);
               if  (Adj == 0)
                   printf (" %7s", "none");
                 else
                   printf (" %7ld", - Adj);
               if  (A [i] . Wrap_Here)
                   printf (" %6s %6ld\n",
                       "-",
                       A [i] . Start2 + Adj - A [Prev] . Start2 - A [Prev] . Len);
                 else
                   printf (" %6ld %6ld\n",
                       A [i] . Start1 + Adj - A [Prev] . Start1 - A [Prev] . Len,
                       A [i] . Start2 + Adj - A [Prev] . Start2 - A [Prev] . Len);
              }
            else
              {
               Adj = A [i] . Simple_Adj;
               printf ("%8ld %8ld %6ld",
                   A [i] . Start1 + Adj, A [i] . Start2 + Adj,
                   A [i] . Len - Adj);
               if  (Adj == 0)
                   printf (" %7s", "none");
                 else
                   printf (" %7ld", - Adj);
               printf (" %6ld %6ld\n",
                   A [i] . Start1 + Adj - A [Prev] . Start1 - A [Prev] . Len,
                   A [i] . Start2 + Adj - A [Prev] . Start2 - A [Prev] . Len);
              }
          Prev = i;
         }

   printf ("> %s %s Other matches\n", File_Name,
                  Use_Reverse_Complement ? "reverse" : "");
   Prev = -1;
   for  (i = 0;  i < N;  i ++)
     if  (! A [i] . Good)
         {
          if  (Prev == -1)
              printf ("%8ld %8ld %6ld %7s %6s %6s\n",
                  A [i] . Start1, A [i] . Start2, A [i] . Len,
                  "none", "-", "-");
            else
              {
               if  (A [i] . Simple_From == Prev)
                   Adj = A [i] . Simple_Adj;
                 else
                   Adj = 0;
               printf ("%8ld %8ld %6ld",
                   A [i] . Start1 + Adj, A [i] . Start2 + Adj,
                   A [i] . Len - Adj);
               if  (Adj == 0)
                   printf (" %7s", "none");
                 else
                   printf (" %7ld", - Adj);
               if  (A [i] . Simple_From == Prev)
                   printf (" %6ld %6ld\n",
                       A [i] . Start1 + Adj - A [Prev] . Start1 - A [Prev] . Len,
                       A [i] . Start2 + Adj - A [Prev] . Start2 - A [Prev] . Len);
                 else
                   printf (" %6s %6s\n", "-", "-");
              }
          Prev = i;
         }

   return  0;
  }
