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

#include <mummer/tigrinc.hh>
#include <mummer/mgaps.hh>


using namespace mummer::mgaps;

static const int      DEFAULT_FIXED_SEPARATION  = 5;
static const long int DEFAULT_MAX_SEPARATION    = 1000;
static const long int DEFAULT_MIN_OUTPUT_SCORE  = 200;
static const double   DEFAULT_SEPARATION_FACTOR = 0.05;

static bool     Check_Labels      = false;
static int      Fixed_Separation  = DEFAULT_FIXED_SEPARATION;
static long int Max_Separation    = DEFAULT_MAX_SEPARATION;
static long int Min_Output_Score  = DEFAULT_MIN_OUTPUT_SCORE;
static double   Separation_Factor = DEFAULT_SEPARATION_FACTOR;
static bool     Use_Extents       = false;


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

  std::vector<Match_t>          A(1);
  UnionFind                     UF;

  ClusterMatches clusterer(Fixed_Separation, Max_Separation, Min_Output_Score, Separation_Factor, Use_Extents);

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
    const char* label = header.c_str();
    clusterer.Cluster_each(A.data(), UF, A.size() - 1, [&](const cluster_type&& cl) {
        clusterer.Print_Cluster(cl, label, std::cout);
        label = "#";
      });
    if(label == header.c_str()) // Empty cluster, output empty header
      std::cout << label << '\n';
  }

   return  0;
  }

