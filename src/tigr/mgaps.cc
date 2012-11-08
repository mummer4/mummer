/* Programmer:  A. Delcher
*        File:  mgaps.cc
*
*  This program reads lists of unique matches between a sequence of strings
*  and a reference string.  For each string in the sequence, it clusters
*  the matches together into groups that may represent longer, inexact
*  matches.
*/


#include  "tigrinc.hh"


const int  DEFAULT_FIXED_SEPARATION = 5;
const long int  DEFAULT_MAX_SEPARATION = 1000;
const long int  DEFAULT_MIN_OUTPUT_SCORE = 200;
const double  DEFAULT_SEPARATION_FACTOR = 0.05;


struct  Match_t
  {
   long int  Start1, Start2, Len;
   long int  Simple_Score;
   long int  Simple_From;
   long int  Simple_Adj;
   int  cluster_id : 30;
   unsigned int  Good : 1;
   unsigned int  Tentative : 1;
  };


inline
long int  Max  (long int A, long int B)

//  Return the larger of  A  and  B .

  {
   if  (A < B)
       return  B;
     else
       return  A;
  }


static int  Check_Labels = FALSE;
static int  Fixed_Separation = DEFAULT_FIXED_SEPARATION;
static long int  Max_Separation = DEFAULT_MAX_SEPARATION;
static long int  Min_Output_Score = DEFAULT_MIN_OUTPUT_SCORE;
static double  Separation_Factor = DEFAULT_SEPARATION_FACTOR;
static int  * UF = NULL;
static int  Use_Extents = FALSE;
  // If TRUE use end minus start as length of cluster instead of
  // sum of component lengths


static int  By_Start2
    (const void * A, const void * B);
static int  By_Cluster
    (const void * A, const void * B);
static void  Filter_Matches
    (Match_t * A, int & N);
static int  Find
    (int a);
static void  Parse_Command_Line
    (int argc, char * argv []);
static void  Process_Matches
    (Match_t * A, int N, char * label);
static int  Process_Cluster
    (Match_t * A, int N, char * label);
static void  Union
    (int a, int b);
static void  Usage
    (char * command);




int  main
    (int argc, char * argv [])

  {
   Match_t  * A = NULL;
   char  line [MAX_LINE];
   char  save [MAX_LINE];
   int  first = TRUE;
   int  header_line_ct = 0;
   long int  S1, S2, Len;
   long int  N = 0, Size = 0;

   Parse_Command_Line  (argc, argv);

   Size = 500;
   A = (Match_t *) Safe_malloc (Size * sizeof (Match_t));
   UF = (int *) Safe_malloc (Size * sizeof (int));

   while  (fgets (line, MAX_LINE, stdin) != NULL)
     {
      if  (line [0] == '>')
          {
           if  (first)
               first = FALSE;
             else
               Process_Matches (A, N, save);
           N = 0;
           strcpy (save, line);
           if  (Check_Labels && (++ header_line_ct % 2 == 0))
               assert (strstr (line, "Reverse") != NULL);
          }
      else if  (sscanf (line, "%ld %ld %ld", & S1, & S2, & Len) == 3)
          {
           if  (N >= Size - 1)
               {
                Size *= 2;
                A = (Match_t *) Safe_realloc (A, Size * sizeof (Match_t));
                UF = (int *) Safe_realloc (UF, Size * sizeof (int));
               }
           N ++;
           A [N] . Start1 = S1;
           A [N] . Start2 = S2;
           A [N] . Len = Len;
           A [N] . Good = FALSE;
           A [N] . Tentative = FALSE;
          }
     }

   Process_Matches (A, N, save);


#if  0
   printf ("> Other matches\n");
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
#endif

   return  0;
  }



static int  By_Start2
    (const void * A, const void * B)

//  Return how  A  and  B  compare if converted to  Match_t
//  based on  Start2  value.  If  Start2  values are equal use
//  Start1  values for comparison.

  {
   Match_t  * x, * y;

   x = (Match_t *) A;
   y = (Match_t *) B;

   if  (x -> Start2 < y -> Start2)
       return  -1;
   else if  (x -> Start2 > y -> Start2)
       return  1;
   else if  (x -> Start1 < y -> Start1)
       return  -1;
   else if  (x -> Start1 > y -> Start1)
       return  1;
     else
       return  0;
  }



static int  By_Cluster
    (const void * A, const void * B)

//  Return how  A  and  B  compare if converted to  Match_t
//  first based on  cluster_id  value, then by  Start2  value,
//  then by  Start1  value.

  {
   Match_t  * x, * y;

   x = (Match_t *) A;
   y = (Match_t *) B;

   if  (x -> cluster_id < y -> cluster_id)
       return -1;
   else if  (x -> cluster_id > y -> cluster_id)
       return  1;
   else if  (x -> Start2 < y -> Start2)
       return  -1;
   else if  (x -> Start2 > y -> Start2)
       return  1;
   else if  (x -> Start1 < y -> Start1)
       return  -1;
   else if  (x -> Start1 > y -> Start1)
       return  1;
     else
       return  0;
  }



static void  Filter_Matches
    (Match_t * A, int & N)

//  Remove from  A [0 .. (N - 1)]  any matches that are internal to a repeat,
//  e.g., if seq1 has 27 As and seq2 has 20 then the first and
//  last matches will be kept, but the 6 matches in the middle will
//  be eliminated.  Also combine overlapping matches on the same
//  diagonal.  Pack all remaining matches into the front of  A  and
//  reduce the value of  N  if any matches are removed.
//  Matches in  A  *MUST* be sorted by  Start2  value.

  {
   int  i, j;

   for  (i = 0;  i < N;  i ++)
     A [i] . Good = TRUE;

   for  (i = 0;  i < N - 1;  i ++)
     {
      int  i_diag, i_end;

      if  (! A [i] . Good)
          continue;

      i_diag = A [i] . Start2 - A [i] . Start1;
      i_end = A [i] . Start2 + A [i] . Len;

      for  (j = i + 1;  j < N && A [j] . Start2 <= i_end;  j ++)
        {
         int  olap;
         int  j_diag;

         assert (A [i] . Start2 <= A [j] . Start2);

         if  (! A [j] . Good)
             continue;

         j_diag = A [j] . Start2 - A [j] . Start1;
         if  (i_diag == j_diag)
             {
              int  j_extent;

              j_extent = A [j] . Len + A [j] . Start2 - A [i] . Start2;
              if  (j_extent > A [i] . Len)
                  {
                   A [i] . Len = j_extent;
                   i_end = A [i] . Start2 + j_extent;
                  }
              A [j] . Good = FALSE;
             }
         else if  (A [i] . Start1 == A [j] . Start1)
             {
              olap = A [i] . Start2 + A [i] . Len - A [j] . Start2;
              if  (A [i] . Len < A [j] . Len)
                  {
                   if  (olap >=  A [i] . Len / 2)
                       {
                        A [i] . Good = FALSE;
                        break;
                       }
                  }
              else if  (A [j] . Len < A [i] . Len)
                  {
                   if  (olap >= A [j] . Len / 2)
                       {
                        A [j] . Good = FALSE;
                       }
                  }
                else
                  {
                   if  (olap >= A [i] . Len / 2)
                       {
                        A [j] . Tentative = TRUE;
                        if  (A [i] . Tentative)
                            {
                             A [i] . Good = FALSE;
                             break;
                            }
                       }
                  }
             }
         else if  (A [i] . Start2 == A [j] . Start2)
             {
              olap = A [i] . Start1 + A [i] . Len - A [j] . Start1;
              if  (A [i] . Len < A [j] . Len)
                  {
                   if  (olap >=  A [i] . Len / 2)
                       {
                        A [i] . Good = FALSE;
                        break;
                       }
                  }
              else if  (A [j] . Len < A [i] . Len)
                  {
                   if  (olap >= A [j] . Len / 2)
                       {
                        A [j] . Good = FALSE;
                       }
                  }
                else
                  {
                   if  (olap >= A [i] . Len / 2)
                       {
                        A [j] . Tentative = TRUE;
                        if  (A [i] . Tentative)
                            {
                             A [i] . Good = FALSE;
                             break;
                            }
                       }
                  }
             }
        }
     }

   for  (i = j = 0;  i < N;  i ++)
     if  (A [i] . Good)
         {
          if  (i != j)
              A [j] = A [i];
          j ++;
         }
   N = j;

   for  (i = 0;  i < N;  i ++)
     A [i] . Good = FALSE;

   return;
  }



static int  Find
    (int a)

//  Return the id of the set containing  a  in  UF .

  {
   int  i, j, k;

   if  (UF [a] < 0)
       return  a;

   for  (i = a;  UF [i] > 0;  i = UF [i])
     ;
   for  (j = a;  UF [j] != i;  j = k)
     {
      k = UF [j];
      UF [j] = i;
     }

   return  i;
  }



static void  Parse_Command_Line
    (int argc, char * argv [])

//  Get options and parameters from command line with  argc
//  arguments in  argv [0 .. (argc - 1)] .

  {
   int  ch, errflg = FALSE;
   char  * p;

   optarg = NULL;

   while  (! errflg
             && ((ch = getopt (argc, argv, "Cd:ef:l:s:")) != EOF))
     switch  (ch)
       {
        case  'C' :
          Check_Labels = TRUE;
          break;

        case  'd' :
          Fixed_Separation = strtol (optarg, & p, 10);
          break;

        case  'e' :
          Use_Extents = TRUE;
          break;

        case  'f' :
          Separation_Factor = strtod (optarg, & p);
          break;

        case  'l' :
          Min_Output_Score = strtol (optarg, & p, 10);
          break;

        case  's' :
          Max_Separation = strtol (optarg, & p, 10);
          break;

        case  '?' :
          fprintf (stderr, "Unrecognized option -%c\n", optopt);

        default :
          errflg = TRUE;
       }

   if  (errflg || optind != argc)
       {
        Usage (argv [0]);
        exit (EXIT_FAILURE);
       }

   return;
  }



static int  Process_Cluster
    (Match_t * A, int N, char * label)

//  Process the cluster of matches in  A [0 .. (N - 1)]  and output them
//  after a line containing  label .  Return the number of clusters
//  printed.

  {
   long int  adj, total, hi, lo, extent, score;
   int  best, prev;
   int  print_ct = 0;
   int  i, j, k;
     
   do
     {
      for  (i = 0;  i < N;  i ++)
        {
         A [i] . Simple_Score = A [i] . Len;
         A [i] . Simple_Adj = 0;
         A [i] . Simple_From = -1;
         for  (j = 0;  j < i;  j ++)
           {
            long int  Pen;
            long int  Olap, Olap1, Olap2;

            Olap1 = A [j] . Start1 + A [j] . Len - A [i] . Start1;
            Olap = Max (0, Olap1);
            Olap2 = A [j] . Start2 + A [j] . Len - A [i] . Start2;
            Olap = Max (Olap, Olap2);

            // penalize off diagonal matches
            Pen = Olap + abs ( (A [i] . Start2 - A [i] . Start1) -
                               (A [j] . Start2 - A [j] . Start1) );

            if  (A [j] . Simple_Score + A [i] . Len - Pen > A [i] . Simple_Score)
                {
                 A [i] . Simple_From = j;
                 A [i] . Simple_Score = A [j] . Simple_Score + A [i] . Len - Pen;
                 A [i] . Simple_Adj = Olap;
                }
           }
        }

      best = 0;
      for  (i = 1;  i < N;  i ++)
        if  (A [i] . Simple_Score > A [best] . Simple_Score)
            best = i;
      total = 0;
      hi = LONG_MIN;
      lo = LONG_MAX;
      for  (i = best;  i >= 0;  i = A [i] . Simple_From)
        {
         A [i] . Good = TRUE;
         total += A [i] . Len;
         if  (A [i] . Start1 + A [i] . Len > hi)
             hi = A [i] . Start1 + A [i] . Len;
         if  (A [i] . Start1 < lo)
             lo = A [i] . Start1;
        }
      extent = hi - lo;

      if  (Use_Extents)
          score = extent;
        else
          score = total;
      if  (score >= Min_Output_Score)
          {
           print_ct ++;
           fputs (label, stdout);
           prev = -1;
           for  (i = 0;  i < N;  i ++)
             if  (A [i] . Good)
                 {
                  if  (prev == -1)
                      printf ("%8ld %8ld %6ld %7s %6s %6s\n",
                          A [i] . Start1, A [i] . Start2, A [i] . Len,
                          "none", "-", "-");
                    else
                      {
                       adj = A [i] . Simple_Adj;
                       printf ("%8ld %8ld %6ld",
                           A [i] . Start1 + adj, A [i] . Start2 + adj,
                           A [i] . Len - adj);
                       if  (adj == 0)
                           printf (" %7s", "none");
                         else
                           printf (" %7ld", - adj);
                       printf (" %6ld %6ld\n",
                           A [i] . Start1 + adj - A [prev] . Start1 - A [prev] . Len,
                           A [i] . Start2 + adj - A [prev] . Start2 - A [prev] . Len);
                      }
                  prev = i;
                 }
           label = "#\n";
          }

      for  (i = k = 0;  i < N;  i ++)
        if  (! A [i] . Good)
            {
             if  (i != k)
                 {
                  A [k] = A [i];
                 }
             k ++;
            }
      N = k;
     }  while  (N > 0);

   return  print_ct;
  }




static void  Process_Matches
    (Match_t * A, int N, char * label)

//  Process matches  A [1 .. N]  and output them after
//  a line containing  label .

  {
   long int  cluster_size, sep;
   int  print_ct = 0;
   int  a, b, i, j;

   if  (N <= 0)
       {
        fputs (label, stdout);
        return;
       }

   //  Use Union-Find to create connected-components based on
   //  separation and similar diagonals between matches

   for  (i = 1;  i <= N;  i ++)
     UF [i] = -1;

   qsort (A + 1, N, sizeof (Match_t), By_Start2);

   Filter_Matches (A + 1, N);

   for  (i = 1;  i < N;  i ++)
     {
      long int  i_end = A [i] . Start2 + A [i] . Len;
      long int  i_diag = A [i] . Start2 - A [i] . Start1;

      for  (j = i + 1;  j <= N;  j ++)
        {
         long int  diag_diff;

         sep = A [j] . Start2 - i_end;
         if  (sep > Max_Separation)
             break;

         diag_diff = abs ((A [j] . Start2 - A [j] . Start1) - i_diag);
         if  (diag_diff <= Max (Fixed_Separation, long (Separation_Factor * sep)))
             {
              a = Find (i);
              b = Find (j);
              if  (a != b)
                  Union (a, b);
             }
        }
     }

   //  Set the cluster id of each match

   for  (i = 1;  i <= N;  i ++)
     A [i] . cluster_id = Find (i);

   qsort (A + 1, N, sizeof (Match_t), By_Cluster);

   for  (i = 1;  i <= N;  i += cluster_size)
     {

      for  (j = i + 1;  j <= N && A [i] . cluster_id == A [j] . cluster_id;  j ++)
        ;
      cluster_size = j - i;
      print_ct += Process_Cluster (A + i, cluster_size, label);
      if  (print_ct > 0)
          label = "#\n";
     }

#if  0
   //  Find the largest cluster

   j = 1;
   for  (i = 2;  i <= N;  i ++)
     if  (UF [i] < UF [j])
         j = i;

   //  j is now the largest cluster

   cluster_size = - UF [j];

   for  (i = k = 1;  i <= N;  i ++)
     if  (Find (i) == j)
         {
          if  (i != k)
              {
               Match_t  tmp = A [i];

               A [i] = A [k];
               A [k] = tmp;
              }
          k ++;
         }

   assert (cluster_size == k - 1);
#endif

   if  (print_ct == 0)
       fputs (label, stdout);

   return;
  }



static void  Union
    (int a, int b)

//  Union the sets whose id's are  a  and  b  in  UF .

  {
   assert (UF [a] < 0 && UF [b] < 0);

   if  (UF [a] < UF [b])
       {
        UF [a] += UF [b];
        UF [b] = a;
       }
     else
       {
        UF [b] += UF [a];
        UF [a] = b;
       }

   return;
  }



static void  Usage
    (char * command)

//  Print to stderr description of options and command line for
//  this program.   command  is the command that was used to
//  invoke it.

  {
   fprintf (stderr,
           "USAGE:  %s [-d <DiagDiff>] [-f <DiagFactor>] [-l <MatchLen>]\n"
           "        [-s <MaxSeparation>]\n"
           "\n"
           "Clusters MUMs based on diagonals and separation.\n"
           "Input is from stdin in format produced by mummer.\n"
           "Ouput goes to stdout.\n"
           "\n"
           "Options:\n"
           "-C       Check that fasta header labels alternately have \"Reverse\"\n"
           "-d num   Fixed diagonal difference to join matches\n"
           "-e       Use extent of match (end - start) rather than sum of piece\n"
           "         lengths to determine length of cluster\n"
           "-f num   Fraction of separation for diagonal difference\n"
           "-l num   Minimum length of cluster match\n"
           "-s num   Maximum separation between matches in cluster\n",
           command);

   return;
  }



