/* Programmer:  A. Delcher
*        File:  mgaps.cc
*
*  This program reads lists of unique matches between a sequence of strings
*  and a reference string.  For each string in the sequence, it clusters
*  the matches together into groups that may represent longer, inexact
*  matches.
*/

#include <cassert>
#include <algorithm>
#include <vector>
#include <iostream>
#include <iomanip>
#include <string>

#include  "tigrinc.hh"

namespace mummer_mgaps {

const int  DEFAULT_FIXED_SEPARATION = 5;
const long int  DEFAULT_MAX_SEPARATION = 1000;
const long int  DEFAULT_MIN_OUTPUT_SCORE = 200;
const double  DEFAULT_SEPARATION_FACTOR = 0.05;


// Union find data structures. Indices MUST be in the range [1, s],
// where s is the size given to reset (no checks).
struct UnionFind {
  std::vector<int> m_UF;

  void reset(size_t s) {
    m_UF.resize(0);
    m_UF.resize(s + 1, -1);
    assert(s >= 1 && m_UF[1] == -1);
  }

  int find(int a) { //  Return the id of the set containing  a  in  UF .
    if(m_UF [a] < 0)
     return  a;

   int i;
   for(i = a;  m_UF [i] > 0;  i = m_UF [i])
     ;
   for(int k, j = a;  m_UF [j] != i;  j = k) {
      k = m_UF [j];
      m_UF [j] = i;
   }
   return  i;
  }


  void union_sets(int a, int b) { //  Union the sets whose id's are  a  and  b  in  UF .
    if(a == b) return;
    assert (m_UF [a] < 0 && m_UF [b] < 0);

   if(m_UF [a] < m_UF [b]) {
     m_UF [a] += m_UF [b];
     m_UF [b] = a;
   } else {
     m_UF [b] += m_UF [a];
     m_UF [a] = b;
   }
  }
};


struct  Match_t {
  long int Start1, Start2, Len;
  long int Simple_Score;
  long int Simple_From;
  long int Simple_Adj;
  unsigned int  cluster_id : 30;
  unsigned int  Good : 1;
  unsigned int  Tentative : 1;
  Match_t() = default;
  Match_t(long int S1, long int S2, long int L) : Start1(S1), Start2(S2), Len(L), Good(true), Tentative(false) { }
};


static bool     Check_Labels      = false;
static int      Fixed_Separation  = DEFAULT_FIXED_SEPARATION;
static long int Max_Separation    = DEFAULT_MAX_SEPARATION;
static long int Min_Output_Score  = DEFAULT_MIN_OUTPUT_SCORE;
static double   Separation_Factor = DEFAULT_SEPARATION_FACTOR;
static bool     Use_Extents       = false;
  // If true use end minus start as length of cluster instead of
  // sum of component lengths


static void Filter_Matches(Match_t * A, int & N);
static int  Process_Matches(Match_t * A, UnionFind& UF, int N, std::vector<std::vector<Match_t>>& clusters);
static int  Process_Cluster(Match_t * A, int N, std::vector<std::vector<Match_t>>& cluster);


bool By_Start2(const Match_t& A, const Match_t& B) {
  return (A.Start2 < B.Start2) || (A.Start2 == B.Start2 && A.Start1 < B.Start1);
}

bool By_Cluster(const Match_t& A, const Match_t& B) {
  return (A.cluster_id < B.cluster_id) ||
    (A.cluster_id == B.cluster_id && By_Start2(A, B));
}

static void  Filter_Matches(Match_t * A, int & N) {

//  Remove from  A [0 .. (N - 1)]  any matches that are internal to a repeat,
//  e.g., if seq1 has 27 As and seq2 has 20 then the first and
//  last matches will be kept, but the 6 matches in the middle will
//  be eliminated.  Also combine overlapping matches on the same
//  diagonal.  Pack all remaining matches into the front of  A  and
//  reduce the value of  N  if any matches are removed.
//  Matches in  A  *MUST* be sorted by  Start2  value.
  for  (int i = 0;  i < N - 1;  i ++) {
    if  (! A[i].Good) continue;

    const int i_diag = A[i].Start2 - A[i].Start1;
    int       i_end  = A[i].Start2 + A[i].Len;

    for  (int j = i + 1;  j < N && A[j].Start2 <= i_end;  j ++) {
      assert (A[i].Start2 <= A[j].Start2);
      if  (! A[j].Good) continue;
      int j_diag = A[j].Start2 - A[j].Start1;
      if  (i_diag == j_diag) {
        int  j_extent = A[j].Len + A[j].Start2 - A[i].Start2;
        if  (j_extent > A[i].Len) {
          A[i].Len = j_extent;
          i_end = A[i].Start2 + j_extent;
        }
        A[j].Good = false;
      } else if  (A[i].Start1 == A[j].Start1) {
        int olap = A[i].Start2 + A[i].Len - A[j].Start2;
        if  (A[i].Len < A[j].Len) {
          if  (olap >=  A[i].Len / 2) {
            A[i].Good = false;
            break;
          }
        } else if  (A[j].Len < A[i].Len) {
          if  (olap >= A[j].Len / 2)
            A[j].Good = false;
        } else {
          if  (olap >= A[i].Len / 2) {
            A[j].Tentative = true;
            if  (A[i].Tentative) {
              A[i].Good = false;
              break;
            }
          }
        }
      } else if  (A[i].Start2 == A[j].Start2) {
        int olap = A[i].Start1 + A[i].Len - A[j].Start1;
        if  (A[i].Len < A[j].Len) {
          if  (olap >=  A[i].Len / 2) {
            A[i].Good = false;
            break;
          }
        } else if  (A[j].Len < A[i].Len) {
          if  (olap >= A[j].Len / 2)
            A[j].Good = false;

        } else {
          if  (olap >= A[i].Len / 2) {
            A[j].Tentative = true;
            if  (A[i].Tentative) {
              A[i].Good = false;
              break;
            }
          }
        }
      }
    }
  }

  int j = 0;
  for  (int i = 0;  i < N;  i ++) {
    if  (A[i].Good) {
      if  (i != j)
        A[j] = A[i];
      j ++;
    }
  }
  N = j;
}

// static std::vector<Match_t> Process_Cluster(Match_t * A, int N) {
//   std::vector<Match_t> res;
//   Process_Cluster(A, N, res);
//   return res;
// }

static int  Process_Cluster(Match_t * A, int N, std::vector<std::vector<Match_t>>& out) {
//  Process the cluster of matches in  A [0 .. (N - 1)]  and output them
//  after a line containing  label .  Return the number of clusters
//  printed.
  int  count = 0;

  do {
    for  (int i = 0;  i < N;  i ++) {
      A [i] . Simple_Score = A [i] . Len;
      A [i] . Simple_Adj = 0;
      A [i] . Simple_From = -1;
      for  (int j = 0;  j < i;  j ++) {
        const long int Olap1 = A [j] . Start1 + A [j] . Len - A [i] . Start1;
        const long int Olap2 = A [j] . Start2 + A [j] . Len - A [i] . Start2;
        const long int Olap = std::max(std::max((long)0, Olap1), Olap2);

          // penalize off diagonal matches
        const long int Pen = Olap + std::abs ( (A [i] . Start2 - A [i] . Start1) -
                                               (A [j] . Start2 - A [j] . Start1) );

        if  (A [j] . Simple_Score + A [i] . Len - Pen > A [i] . Simple_Score) {
          A [i] . Simple_From = j;
          A [i] . Simple_Score = A [j] . Simple_Score + A [i] . Len - Pen;
          A [i] . Simple_Adj = Olap;
        }
      }
    }

    int best = 0;
    for  (int i = 1;  i < N;  i ++)
      if  (A [i] . Simple_Score > A [best] . Simple_Score)
        best = i;
    long int total = 0;
    long int hi    = LONG_MIN;
    long int lo    = LONG_MAX;
    for  (int i = best;  i >= 0;  i = A [i] . Simple_From) {
      A [i] . Good = true;
      total += A [i] . Len;
      hi = std::max(hi, A[i].Start1 + A[i].Len);
      lo = std::min(lo, A[i].Start1);
    }
    const long int score = Use_Extents ? hi - lo : total;

    if  (score >= Min_Output_Score) {
      count ++;
      std::vector<Match_t>* cluster = nullptr;
      for  (int i = 0;  i < N;  i ++)
        if  (A [i] . Good) {
          if(!cluster) {
            out.push_back(std::vector<Match_t>());
            cluster = &out.back();
          }
          assert(cluster != nullptr);
          cluster->push_back(A[i]);
        }
    }

    int k = 0;
    for  (int i = 0;  i < N;  i ++) {
      if  (! A [i] . Good) {
        if  (i != k)
          A [k] = A [i];
        k ++;
      }
    }
    N = k;
  }  while  (N > 0);

  return count;
}


static void Print_Cluster(const std::vector<std::vector<Match_t>>& clusters, const char* label) {
  for(const auto cl : clusters) {
    std::cout << label << '\n'
              << std::setw(8) << cl[0].Start1 << ' '
              << std::setw(8) << cl[0].Start2 << ' '
              << std::setw(6) << cl[0].Len << ' '
              << "   none      -      -\n";
    for(size_t i = 1; i < cl.size(); ++i) {
      const int long adj = cl[i] . Simple_Adj;
      std::cout << std::setw(8) << (cl[i].Start1 + adj) << ' '
                << std::setw(8) << (cl[i].Start2 + adj) << ' '
                << std::setw(6) << (cl[i].Len - adj) << ' ';
      if  (adj == 0)
        std::cout << std::setw(7) << "none" << ' ';
      else
        std::cout << std::setw(7) << (-adj) << ' ';
      std::cout << std::setw(6) << (cl[i].Start1 + adj - cl[i-1].Start1 - cl[i-1].Len) << ' '
                << std::setw(6) << (cl[i].Start2 + adj - cl[i-1].Start2 - cl[i-1].Len) << '\n';
    }
    label = "#";
  }
}

static int Process_Matches(Match_t * A, UnionFind& UF, int N, std::vector<std::vector<Match_t>>& clusters) {
  //  Process matches  A [1 .. N]  and output them after
  //  a line containing  label .

  //  Use Union-Find to create connected-components based on
  //  separation and similar diagonals between matches
  UF.reset(N);

  std::sort(A + 1, A + N + 1, By_Start2);
  Filter_Matches (A + 1, N);

  for  (int i = 1;  i < N;  i ++) {
    long int i_end  = A [i] . Start2 + A [i] . Len;
    long int i_diag = A [i] . Start2 - A [i] . Start1;

    for  (int j = i + 1;  j <= N;  j ++) {
      long int sep = A [j] . Start2 - i_end;
      if  (sep > Max_Separation)
        break;

      long int diag_diff = std::abs ((A [j] . Start2 - A [j] . Start1) - i_diag);
      if  (diag_diff <= std::max(Fixed_Separation, (int)(Separation_Factor * sep)))
        UF.union_sets(UF.find(i), UF.find(j));
    }
  }

  //  Set the cluster id of each match and reset Good flag
  for  (int i = 1;  i <= N;  i ++) {
    A [i] . cluster_id = UF.find (i);
    assert(A[i].cluster_id > 0);
    A[i].Good = false;
  }
  std::sort(A + 1, A + N + 1, By_Cluster);

  // Determine and process clusters
  int cluster_size, print_ct = 0;
  for (int i = 1;  i <= N;  i += cluster_size) {
    int j;
    for  (j = i + 1;  j <= N && A [i] . cluster_id == A [j] . cluster_id;  j ++)
      ;
    cluster_size = j - i;
    print_ct += Process_Cluster (A + i, cluster_size, clusters);
  }
  return print_ct;
}

} // namespace mummer_mgaps

using namespace mummer_mgaps;

static void  Usage
    (char * command)

//  Print to stderr description of options and command line for
//  this program.   command  is the command that was used to
//  invoke it.

  {
    std::cerr <<
      "USAGE:  " << command << " [-d <DiagDiff>] [-f <DiagFactor>] [-l <MatchLen>]\n"
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
      "-s num   Maximum separation between matches in cluster\n";

   return;
  }

static void  Parse_Command_Line
    (int argc, char * argv [])

//  Get options and parameters from command line with  argc
//  arguments in  argv [0 .. (argc - 1)] .

  {
   int  ch, errflg = false;
   char  * p;

   optarg = NULL;

   while  (! errflg
             && ((ch = getopt (argc, argv, "Cd:ef:l:s:")) != EOF))
     switch  (ch)
       {
        case  'C' :
          Check_Labels = true;
          break;

        case  'd' :
          Fixed_Separation = strtol (optarg, & p, 10);
          break;

        case  'e' :
          Use_Extents = true;
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
          std::cerr << "Unrecognized option -" << optopt << '\n';

        default :
          errflg = true;
       }

   if  (errflg || optind != argc)
       {
        Usage (argv [0]);
        exit (EXIT_FAILURE);
       }

   return;
  }


int  main(int argc, char * argv []) {
  std::ios::sync_with_stdio(false);

  std::string line, header;
  int  header_line_ct = 0;
  long int  S1, S2, Len;

  Parse_Command_Line  (argc, argv);

  std::vector<Match_t> A(1);
  UnionFind UF;
  std::vector<std::vector<Match_t>> clusters;

  int c = std::cin.peek();
  while(c != '>' && c != EOF) // Skip to first header
    std::getline(std::cin, line);
  while(c != EOF) {
    std::getline(std::cin, header); // Get header
    if(Check_Labels && (++ header_line_ct % 2 == 0))
      assert (strstr (header.c_str(), "Reverse") != NULL);
    A.resize(1);
    for(c = std::cin.peek(); c != '>' && c != EOF; c = std::cin.peek()) {
      std::getline(std::cin, line);
      if  (sscanf (line.c_str(), "%ld %ld %ld", & S1, & S2, & Len) == 3)
        A.push_back(Match_t(S1, S2, Len));
    }
    clusters.clear();
    const int nb_elements = Process_Matches (A.data(), UF, A.size() - 1, clusters);
    if(nb_elements > 0)
      Print_Cluster(clusters, header.c_str());
    else
      std::cout << header << '\n';
  }

   return  0;
  }

