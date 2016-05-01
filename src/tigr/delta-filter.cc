//------------------------------------------------------------------------------
//   Programmer: Adam M Phillippy, The Institute for Genomic Research
//         File: delta-filter.cc
//         Date: 09 / 22 / 2004
//
//        Usage: delta-filter  [options]  <deltafile>
//               Try 'show-coords -h' for more information
//
//  Description: For use in conjunction with the MUMmer package.
//              "delta-filter" cleans delta alignment files by filtering
//             alignments that fail to meet the specifications required by the
//            command line switches. Produces a new, filtered delta file on
//           stdout, and works for both nucleotide and amino-acid alignments.
//
//------------------------------------------------------------------------------

#include <mummer/delta.hh>
#include <mummer/tigrinc.hh>
#include <map>
#include <vector>
#include <string>
#include <sstream>
#include <iostream>
#include <algorithm>
using namespace std;


//=============================================================== Options ====//
string         OPT_AlignName;                // alignment name parameter

bool           OPT_QLIS         = false;     // do query based LIS
bool           OPT_RLIS         = false;     // do reference based LIS
bool           OPT_GLIS         = false;     // do global LIS
bool           OPT_1to1         = false;     // do 1-to-1 alignment
bool           OPT_MtoM         = false;     // do M-to-M alignment
long int       OPT_MinLength    = 0;         // minimum alignment length
float          OPT_MinIdentity  = 0.0;       // minimum %identity
float          OPT_MinUnique    = 0.0;       // minimum %unique
float          OPT_MaxOverlap   = 100.0;     // maximum olap as % of align len
float          OPT_Epsilon      = -1.0;      // negligible alignment score


//========================================================== Fuction Decs ====//

//------------------------------------------------------------- ParseArgs ----//
void ParseArgs(int argc, char ** argv);


//------------------------------------------------------------- PrintHelp ----//
void PrintHelp(const char * s);


//------------------------------------------------------------ PrintUsage ----//
void PrintUsage(const char * s);

//----------------------------------------------------------- FilterDelta ----//
int FilterDelta(std::istream& is, float min_len, float min_idy);

//========================================================= Function Defs ====//
int main(int argc, char ** argv)
{
  srand(1);

  //-- Command line parsing
  ParseArgs(argc, argv);

  // Optimization: if only filtering by identity or length, stream
  // through, do not create graph.
  const bool stream = !(OPT_MinUnique > 0) && !OPT_QLIS && !OPT_RLIS && !OPT_GLIS && !OPT_MtoM && !OPT_1to1;
  const bool noop = stream && !(OPT_MinIdentity > 0) && !(OPT_MinLength);
  if(stream) {
    std::ifstream is(OPT_AlignName);
    if(!is.good()) {
      std::cerr << "Error opening delta file '" << OPT_AlignName << "'" << std::endl;
      return EXIT_FAILURE;
    }

    if(noop) { // No operation! Copy input to output
      std::cout << is.rdbuf();
      return EXIT_SUCCESS;
    }
    return FilterDelta(is, OPT_MinLength, OPT_MinIdentity);
  }

  //-- Build the alignment graph from the delta file
  DeltaGraph_t graph;
  graph.build(OPT_AlignName, true);

  //-- Identity requirements
  if ( OPT_MinIdentity > 0  ||  OPT_MinLength > 0 )
    graph.flagScore(OPT_MinLength, OPT_MinIdentity);

  //-- Uniqueness requirements
  if ( OPT_MinUnique > 0 )
    graph.flagUNIQ(OPT_MinUnique);

  //-- Query-based LIS
  if ( OPT_QLIS )
    graph.flagQLIS(OPT_Epsilon, OPT_MaxOverlap);

  //-- Reference-based LIS
  if ( OPT_RLIS )
    graph.flagRLIS(OPT_Epsilon, OPT_MaxOverlap);

  //-- Global LIS
  if ( OPT_GLIS )
    graph.flagGLIS(OPT_Epsilon);

  //-- M-to-M
  if ( OPT_MtoM )
    graph.flagMtoM(OPT_Epsilon, OPT_MaxOverlap);

  //-- 1-to-1
  if ( OPT_1to1 )
    graph.flag1to1(OPT_Epsilon, OPT_MaxOverlap);

  //-- Output the filtered delta file
  graph.outputDelta(cout);

  return EXIT_SUCCESS;
}




//------------------------------------------------------------- ParseArgs ----//
void ParseArgs(int argc, char ** argv)
{
  int ch, errflg = 0;
  optarg = NULL;
  
  while ( !errflg  &&
         ((ch = getopt(argc, argv, "e:ghi:l:o:qru:m1")) != EOF) )
    switch (ch)
      {
      case 'e':
        OPT_Epsilon = atof(optarg);
        break;

      case 'g':
        OPT_GLIS = true;
        break;

      case 'h':
        PrintHelp(argv[0]);
        exit(EXIT_SUCCESS);
        break;

      case 'i':
        OPT_MinIdentity = atof(optarg);
        break;

      case 'l':
        OPT_MinLength = atol(optarg);
        break;

      case 'o':
	OPT_MaxOverlap = atof(optarg);
	break;

      case 'q':
        OPT_QLIS = true;
        break;

      case 'r':
        OPT_RLIS = true;
        break;

      case 'u':
        OPT_MinUnique = atof(optarg);
        break;

      case 'm':
        OPT_MtoM = true;
        break;

      case '1':
        OPT_1to1 = true;
        break;

      default:
        errflg ++;
      }

  if ( OPT_MinIdentity < 0.0  ||  OPT_MinIdentity > 100.0 )
    {
      cerr << "ERROR: Minimum identity must be within the range [0, 100]\n";
      errflg ++;
    }

  if ( OPT_MinLength < 0 )
    {
      cerr << "ERROR: Minimum length must be greater than or equal to zero\n";
      errflg ++;
    }

  if ( OPT_MinUnique < 0.0  ||  OPT_MinUnique > 100.0 )
    {
      cerr << "ERROR: Minimum uniqueness must be within the range [0, 100]\n";
      errflg ++;
    }

  if ( OPT_MaxOverlap < 0.0  ||  OPT_MaxOverlap > 100.0 )
    {
      cerr << "ERROR: Maximum overlap must be within the range [0, 100]\n";
      errflg ++;
    }

  if ( errflg > 0  ||  optind != argc - 1 )
    {
      PrintUsage(argv[0]);
      cerr << "Try '" << argv[0] << " -h' for more information.\n";
      exit(EXIT_FAILURE);
    }

  OPT_AlignName = argv [optind ++];
}




//------------------------------------------------------------- PrintHelp ----//
void PrintHelp(const char * s)
{
  PrintUsage(s);
  cerr
    << "-1            1-to-1 alignment allowing for rearrangements\n"
    << "              (intersection of -r and -q alignments)\n"
    << "-g            1-to-1 global alignment not allowing rearrangements\n"
    << "-h            Display help information\n"
    << "-i float      Set the minimum alignment identity [0, 100], default "
    << OPT_MinIdentity << endl
    << "-l int        Set the minimum alignment length, default "
    << OPT_MinLength << endl
    << "-m            Many-to-many alignment allowing for rearrangements\n"
    << "              (union of -r and -q alignments)\n"
    << "-q            Maps each position of each query to its best hit in\n"
    << "              the reference, allowing for reference overlaps\n"
    << "-r            Maps each position of each reference to its best hit\n"
    << "              in the query, allowing for query overlaps\n"
    << "-u float      Set the minimum alignment uniqueness, i.e. percent of\n"
    << "              the alignment matching to unique reference AND query\n"
    << "              sequence [0, 100], default "
    << OPT_MinUnique << endl
    << "-o float      Set the maximum alignment overlap for -r and -q options\n"
    << "              as a percent of the alignment length [0, 100], default "
    << OPT_MaxOverlap << endl
    << endl;

  cerr
    << "  Reads a delta alignment file from either nucmer or promer and\n"
    << "filters the alignments based on the command-line switches, leaving\n"
    << "only the desired alignments which are output to stdout in the same\n"
    << "delta format as the input. For multiple switches, order of operations\n"
    << "is as follows: -i -l -u -q -r -g -m -1. If an alignment is excluded\n"
    << "by a preceding operation, it will be ignored by the succeeding\n"
    << "operations.\n"
    << "  An important distinction between the -g option and the -1 and -m\n"
    << "options is that -g requires the alignments to be mutually consistent\n"
    << "in their order, while the -1 and -m options are not required to be\n"
    << "mutually consistent and therefore tolerate translocations,\n"
    << "inversions, etc. In general cases, the -m option is the best choice,\n"
    << "however -1 can be handy for applications such as SNP finding which\n"
    << "require a 1-to-1 mapping. Finally, for mapping query contigs, or\n"
    << "sequencing reads, to a reference genome, use -q.\n"
    << endl;

  return;
}


//------------------------------------------------------------ PrintUsage ----//
void PrintUsage(const char * s)
{
  cerr
    << "\nUSAGE: " << s << "  [options]  <deltafile>\n\n";
  return;
}

//----------------------------------------------------------- FilterDelta ----//
int FilterDelta(std::istream& is, float min_len, float min_idy) {
  std::string line;
  // Print header unmodified
  std::getline(is, line);
  std::cout << line << '\n';
  std::getline(is, line);
  const bool promer = line == PROMER_STRING;
  if(!promer && line != NUCMER_STRING) {
    std::cerr << "Unsupported format '" << line << "'\n";
    return EXIT_FAILURE;
  }
  std::cout << line << '\n';

  int c = is.peek();
  if(c != '>') {
    std::cerr << "Invalid format. Expected '>' got '" << c << "'\n";
    return EXIT_FAILURE;
  }

  DeltaRecord_t record;
  DeltaAlignment_t alignment;

  while(c != EOF) {
    is >> record;
    bool first = true;
    for(c = is.peek(); c != '>' && c != EOF; c = is.peek()) {
      alignment.read(is, promer, true);
      if(alignment.idy < min_idy ||
         std::abs(alignment.eR - alignment.sR) + 1 < min_len ||
         std::abs(alignment.eQ - alignment.sQ) + 1 < min_len)
        continue;
      if(first) {
        std::cout << record << '\n';
        first = false;
      }
      std::cout << alignment;
    }
  }
  return EXIT_SUCCESS;
}
