//------------------------------------------------------------------------------
//   Programmer: Adam M Phillippy, The Institute for Genomic Research
//         File: show-aligns.cc
//         Date: 10 / 18 / 2002
//
//        Usage: show-aligns [options] <deltafile>
//               Try 'show-aligns -h' for more information
//
//  Description: For use in conjunction with the MUMmer package.
//              "show-aligns" displays human readable information from the
//             .delta output of the "nucmer" and "promer" programs. Outputs
//            pairwise alignments to stdout. Works for both nucleotide and
//           amino-acid alignments.
//
//------------------------------------------------------------------------------

#include <sys/ioctl.h>

#include <mummer/delta.hh>
#include <mummer/tigrinc.hh>
#include <mummer/translate.hh>
#include <mummer/sw_alignscore.hh>
#include <mummer/redirect_to_pager.hpp>
#include <jellyfish/stream_manager.hpp>
#include <jellyfish/whole_sequence_parser.hpp>
#include <string>
#include <vector>
#include <algorithm>
using namespace std;

//------------------------------------------------------------- Jellyfish parser ----//
typedef std::vector<std::string>::const_iterator         path_iterator;
typedef jellyfish::stream_manager<path_iterator>         stream_manager;
typedef jellyfish::whole_sequence_parser<stream_manager> sequence_parser;


//------------------------------------------------------------- Constants ----//

const char NUCMER_MISMATCH_CHAR = '^';
const char NUCMER_MATCH_CHAR = ' ';
const char PROMER_SIM_CHAR = '+';
const char PROMER_MISMATCH_CHAR = ' ';

// ANSI terminal codes
static const char* const ANSI_RESET     = "\033[0m";
static const char* const ANSI_BOLD    = "\033[1m";
static const char* const ANSI_UNDERLINE = "\033[4m";

//-- Note: if coord exceeds LINE_PREFIX_LEN - 1 digits,
//         increase these accordingly
#define LINE_PREFIX_LEN 11
#define PREFIX_FORMAT "%-10ld "

#define DEFAULT_SCREEN_WIDTH 60
int Screen_Width = 0;

#define DEFAULT_MARKER_WIDTH 10
int Marker_Width = DEFAULT_MARKER_WIDTH;

bool Colorize = false;

class ColoredBuffer {
  std::string m_buf;
  bool        m_bright;
  int         m_color;
  size_t      m_bases;

public:
  ColoredBuffer()
    : m_bright(false)
    , m_color(-1)
    , m_bases(0)
  { }

  // Behave somewhat like a std::string
  void clear() { m_buf.clear(); m_bases = 0; }
  bool empty() const { return m_buf.empty(); }
  size_t size() const { return m_buf.size(); }
  char& back() { return m_buf.back(); }
  const char* c_str() { return m_buf.c_str(); }

  // Number of bases. Ignore formatting characters
  size_t bases() const { return m_bases; }

  void reset() {
    if(m_bright || m_color != -1)
      m_buf += ANSI_RESET;
    m_bright = false;
    m_color = -1;
  }
  template<typename T>
  ColoredBuffer& operator+=(const T& x) {
    const size_t before = size();
    m_buf += x;
    m_bases += size() - before;
    return *this;
  }
  enum colors { Black, Red, Green, Yello, Blue, Magenta, Cyan, White };

  void color(colors c, bool b) {
    if(c == m_color && b == m_bright) return;
    if(m_bright && !b)
      reset();
    m_buf    += ansi_color(c, b);
    m_color   = c;
    m_bright  = b;
  }

  std::string ansi_color(colors c, bool b) {
    std::string res("\033[3");
    res += std::to_string(c);
    if(b) res += ";7"; // 1 for bold, 7 for reverse video
    res += 'm';
    // std::cerr << "ansi_color " << c << ',' << b << ' ' << res << '\n';
    return res;
  }

  bool bright() const { return m_bright; }
};

//-- Get screen width from system. If not available, returns the
//-- default.
int get_screen_width() {
  struct winsize w;
  if(ioctl(1, TIOCGWINSZ, &w) == -1)
    return DEFAULT_SCREEN_WIDTH;
  return w.ws_col;
}



//------------------------------------------------------ Type Definitions ----//
struct AlignStats
     //-- Alignment statistics data structure
{
  long int sQ, eQ, sR, eR;              // start and end in Query and Reference
                                        // relative to the directional strand
  vector<long int> Delta;               // delta information
};



struct sR_Sort
//-- For sorting alignments by their sR coordinate
{
  bool operator( ) (const AlignStats & pA, const AlignStats & pB)
  {
    //-- sort sR
    if ( pA.sR < pB.sR )
      return true;
    else
      return false;
  }
};



struct sQ_Sort
//-- For sorting alignments by their sQ coordinate
{
  bool operator( ) (const AlignStats & pA, const AlignStats & pB)
  {
    //-- sort sQ
    if ( pA.sQ < pB.sQ )
      return true;
    else
      return false;
  }
};




//------------------------------------------------------ Global Variables ----//
bool isSortByQuery = false;              // -q option
bool isSortByReference = false;          // -r option

int DATA_TYPE = NUCMER_DATA;
int MATRIX_TYPE = mummer::sw_align::BLOSUM62;

char InputFileName [MAX_LINE];
// Today, only 1 ref and qry file supported in delta format. It may change one day!
std::vector<std::string> RefFileNames, QryFileNames;


//------------------------------------------------- Function Declarations ----//
long int toFwd
     (long int coord, long int len, int frame);

void parseDelta
(vector<AlignStats> & Aligns, const std::string& IdR, const std::string& IdQ);

void printAlignments
(vector<AlignStats> Aligns, const std::string& R, const std::string& Q);

void add_prefix(ColoredBuffer& Buff, long int pos, long int seqlen, int frame);

void append(ColoredBuffer& Buff1, ColoredBuffer& Buff2, std::string &Buff3,
            char c1, char c2, char c3);

void print_buffers(ColoredBuffer& b1, ColoredBuffer& b2, std::string& b3);

void print_markers(int max_len = Screen_Width);

void printHelp
     (const char * s);

void printUsage
     (const char * s);

long int revC
     (long int coord, long int len);

bool find_sequence(const std::vector<string>& paths, const std::string& IdS, std::string& seq);

//-------------------------------------------------- Function Definitions ----//
int main
     (int argc, char ** argv)
{
  long int i;

  vector<AlignStats> Aligns;

  std::string R, Q; // Reference & Query sequence

  std::string IdR, IdQ;

  //-- Parse the command line arguments
  {
    int ch, errflg = 0;
    optarg = NULL;

    while ( !errflg  &&  ((ch = getopt
                           (argc, argv, "hqrw:x:m:c")) != EOF) )
      switch (ch)
        {
        case 'h' :
	  printHelp (argv[0]);
	  exit (EXIT_SUCCESS);
          break;

	case 'q' :
	  isSortByQuery = true;
	  break;

	case 'r' :
	  isSortByReference = true;
	  break;

	case 'w' :
	  Screen_Width = atoi (optarg);
	  if ( Screen_Width <= LINE_PREFIX_LEN )
	    {
	      fprintf(stderr,
		      "WARNING: invalid screen width %d, using default\n",
		      Screen_Width);
	      Screen_Width = DEFAULT_SCREEN_WIDTH;
	    }
	  break;

        case 'm':
          Marker_Width = atoi(optarg);
          break;

	case 'x' :
	  MATRIX_TYPE = atoi (optarg);
	  if ( MATRIX_TYPE < 1 || MATRIX_TYPE > 3 )
	    {
	      fprintf(stderr,
		      "WARNING: invalid matrix type %d, using default\n",
		      MATRIX_TYPE);
	      MATRIX_TYPE = mummer::sw_align::BLOSUM62;
	    }
	  break;

        case 'c' :
          Colorize = true;
          break;

        default :
          errflg ++;
        }

    if ( errflg > 0  ||  argc - optind != 3 )
      {
        printUsage (argv[0]);
        exit (EXIT_FAILURE);
      }

    if ( isSortByQuery  &&  isSortByReference )
      fprintf (stderr,
               "WARNING: both -r and -q were passed, -q ignored\n");
    if(Screen_Width == 0)
      Screen_Width = get_screen_width();
  }

  strcpy (InputFileName, argv[optind ++]);
  IdR = argv[optind++];
  IdQ = argv[optind++];

  //-- Read in the alignment data
  parseDelta (Aligns, IdR.c_str(), IdQ.c_str());

  //-- Find, and read in the reference sequence
  if(!find_sequence(RefFileNames, IdR, R))
    {
      fprintf(stderr,"ERROR: Could not find %s in the reference file\n", IdR.c_str());
      exit (EXIT_FAILURE);
    }

  if(!find_sequence(QryFileNames, IdQ, Q))
  {
    fprintf(stderr,"ERROR: Could not find %s in the query file\n", IdQ.c_str());
    exit (EXIT_FAILURE);
  }

  //-- Sort the alignment regions if user passed -r or -q option
  if ( isSortByReference )
    sort (Aligns.begin( ), Aligns.end( ), sR_Sort( ));
  else if ( isSortByQuery )
    sort (Aligns.begin( ), Aligns.end( ), sQ_Sort( ));


  //-- Output the alignments to stdout
  stdio_launch_pager redirect_to_pager;
  printf("%s %s\n\n", RefFileNames[0].c_str(), QryFileNames[0].c_str());
  for ( i = 0; i < Screen_Width; i ++ ) printf("=");
  printf("\n-- Alignments between %s and %s\n\n", IdR.c_str(), IdQ.c_str());
  printAlignments (Aligns, R, Q);
  printf("\n");
  for ( i = 0; i < Screen_Width; i ++ ) printf("=");
  printf("\n");
  fclose(stdout);

  return EXIT_SUCCESS;
}




long int toFwd
     (long int coord, long int len, int frame)

     // Switch relative coordinate to reference forward DNA strand

{
  long int newc = coord;

  if ( DATA_TYPE == PROMER_DATA )
    newc = newc * 3 - (3 - labs(frame));

  if ( frame < 0 )
    return revC ( newc, len );
  else
    return newc;
}




void parseDelta
(vector<AlignStats> & Aligns, const std::string& IdR, const std::string& IdQ)

     // Read in the alignments from the desired region

{
  AlignStats aStats;                     //  single alignment region
  bool found = false;

  DeltaReader_t dr;
  dr.open (InputFileName);
  DATA_TYPE = dr.getDataType( ) == NUCMER_STRING ?
    NUCMER_DATA : PROMER_DATA;
  RefFileNames.push_back(dr.getReferencePath());
  QryFileNames.push_back(dr.getQueryPath());

  while ( dr.readNext( ) )
    {
      if ( dr.getRecord( ).idR == IdR  &&
	   dr.getRecord( ).idQ == IdQ )
	{
	  found = true;
	  break;
	}
    }
  if ( !found )
    {
      fprintf(stderr, "ERROR: Could not find any alignments for %s and %s\n",
	      IdR.c_str(), IdQ.c_str());
      exit (EXIT_FAILURE);
    }

  for ( unsigned int i = 0; i < dr.getRecord( ).aligns.size( ); i ++ )
    {
      aStats.sR = dr.getRecord( ).aligns[i].sR;
      aStats.eR = dr.getRecord( ).aligns[i].eR;
      aStats.sQ = dr.getRecord( ).aligns[i].sQ;
      aStats.eQ = dr.getRecord( ).aligns[i].eQ;

      aStats.Delta = dr.getRecord( ).aligns[i].deltas;

      //-- Add the new alignment
      Aligns.push_back (aStats);
    }
  dr.close( );

  return;
}




void printAlignments
(vector<AlignStats> Aligns, const std::string& R, const std::string& Q)

     // Print the alignments to the screen

{
  vector<AlignStats>::iterator Ap;
  vector<long int>::iterator Dp;

  const char * A[7] = {NULL, NULL, NULL, NULL, NULL, NULL, NULL};
  const char * B[7] = {NULL, NULL, NULL, NULL, NULL, NULL, NULL};
  int Ai, Bi, i;

  ColoredBuffer Buff1, Buff2;
  std::string Buff3;

  int Sign;
  long int Delta;
  long int Total, Remain;
  // long int Errors;
  long int Pos;
  char c; // Character to add to Buff3

  long int sR, eR, sQ, eQ;
  long int Apos, Bpos;
  long int SeqLenR, SeqLenQ;
  int frameR, frameQ;

  SeqLenR = R.size() - 1; // strlen (R + 1);
  SeqLenQ = Q.size() - 1; // strlen (Q + 1);

  if ( DATA_TYPE == NUCMER_DATA )
    {
      A[1] = R.c_str();
      char* rcR = (char *) Safe_malloc ( sizeof(char) * (SeqLenR + 2) );
      strcpy ( rcR + 1, A[1] + 1 );
      rcR[0] = '\0';
      Reverse_Complement ( rcR, 1, SeqLenR );
      A[4] = rcR;

      B[1] = Q.c_str();
      char* rcQ = (char *) Safe_malloc ( sizeof(char) * (SeqLenQ + 2) );
      strcpy ( rcQ + 1, B[1] + 1 );
      rcQ[0] = '\0';
      Reverse_Complement ( rcQ, 1, SeqLenQ );
      B[4] = rcQ;
    }

  for ( Ap = Aligns.begin( ); Ap < Aligns.end( ); Ap ++ )
    {
      sR = Ap->sR;
      eR = Ap->eR;
      sQ = Ap->sQ;
      eQ = Ap->eQ;

      //-- Get the coords and frame right
      frameR = 1;
      if ( sR > eR )
	{
	  sR = revC (sR, SeqLenR);
	  eR = revC (eR, SeqLenR);
	  frameR += 3;
	}
      frameQ = 1;
      if ( sQ > eQ )
	{
	  sQ = revC (sQ, SeqLenQ);
	  eQ = revC (eQ, SeqLenQ);
	  frameQ += 3;
	}

      if ( DATA_TYPE == PROMER_DATA )
	{
	  frameR += (sR + 2) % 3;
	  frameQ += (sQ + 2) % 3;

	  //-- Translate the coordinates from DNA to Amino Acid
	  //   remeber that eR and eQ point to the last nucleotide in the codon
	  sR = (sR + 2) / 3;
	  eR = eR / 3;
	  sQ = (sQ + 2) / 3;
	  eQ = eQ / 3;
	}
      Ai = frameR;
      Bi = frameQ;
      if ( frameR > 3 )
	frameR = -(frameR - 3);
      if ( frameQ > 3 )
	frameQ = -(frameQ - 3);

      //-- Get the sequence
      if ( A[Ai] == NULL )
	{
	  assert ( DATA_TYPE == PROMER_DATA );
	  char* translated = (char *) Safe_malloc ( sizeof(char) * ( SeqLenR / 3 + 2 ) );
	  translated[0] = '\0';
	  Translate_DNA ( R.c_str(), R.size() - 1, translated, Ai );
          A[Ai] = translated;
	}
      if ( B[Bi] == NULL )
	{
	  assert ( DATA_TYPE == PROMER_DATA );
          char* translated = (char *) Safe_malloc ( sizeof(char) * ( SeqLenQ / 3 + 2 ) );
	  translated[0] = '\0';
	  Translate_DNA ( Q.c_str(), Q.size() - 1, translated, Bi );
	  B[Bi] = translated;

	}


      //-- Generate the alignment
      printf("-- BEGIN alignment [ %s%d %ld - %ld | %s%d %ld - %ld ]\n\n",
	     frameR > 0 ? "+" : "-", abs(frameR), Ap->sR, Ap->eR,
	     frameQ > 0 ? "+" : "-", abs(frameQ), Ap->sQ, Ap->eQ);

      Apos = sR;
      Bpos = sQ;

      //      Errors = 0;
      Total = 0;
      Remain = eR - sR + 1;

      add_prefix(Buff1, Apos, SeqLenR, frameR);
      add_prefix(Buff2, Bpos, SeqLenQ, frameQ);
      Pos = LINE_PREFIX_LEN;

      for ( Dp = Ap->Delta.begin( );
	    Dp < Ap->Delta.end( ) &&
	    *Dp != 0; Dp ++ )
	{
	  Delta = *Dp;
	  Sign = Delta > 0 ? 1 : -1;
	  Delta = labs ( Delta );


	  //-- For all the bases before the next indel
	  for ( i = 1; i < Delta; i ++ )
	    {
	      if ( Pos >= Screen_Width ) {
                print_buffers(Buff1, Buff2, Buff3);
                add_prefix(Buff1, Apos, SeqLenR, frameR);
                add_prefix(Buff2, Bpos, SeqLenQ, frameQ);
                Pos = LINE_PREFIX_LEN;
              }

	      if ( DATA_TYPE == NUCMER_DATA )
                c = A[Ai][Apos] == B[Bi][Bpos] ? NUCMER_MATCH_CHAR : NUCMER_MISMATCH_CHAR;
              else if ( A[Ai][Apos] == B[Bi][Bpos] )
                c = A[Ai][Apos];
	      else
		c = mummer::sw_align::MATCH_SCORE
		  [MATRIX_TYPE]
		  [toupper(A[Ai][Apos]) - 'A']
		  [toupper(B[Bi][Bpos]) - 'A'] > 0 ?
		  PROMER_SIM_CHAR : PROMER_MISMATCH_CHAR;
              append(Buff1, Buff2, Buff3, A[Ai][Apos++], B[Bi][Bpos++], c);
              ++Pos;
	    }


	  //-- For the indel
	  Remain -= i - 1;

	  if ( Pos >= Screen_Width ) {
            print_buffers(Buff1, Buff2, Buff3);
            add_prefix(Buff1, Apos, SeqLenR, frameR);
            add_prefix(Buff2, Bpos, SeqLenQ, frameQ);
            Pos = LINE_PREFIX_LEN;
          }

	  if ( Sign == 1 )
	    {
	      if ( DATA_TYPE == NUCMER_DATA )
		c = NUCMER_MISMATCH_CHAR;
	      else
		c = PROMER_MISMATCH_CHAR;
              append(Buff1, Buff2, Buff3, A[Ai][Apos++], '.', c);
	      Remain --;
              ++Pos;
	    }
	  else
	    {
	      if ( DATA_TYPE == NUCMER_DATA )
		c = NUCMER_MISMATCH_CHAR;
	      else
		c = PROMER_MISMATCH_CHAR;
              append(Buff1, Buff2, Buff3, '.', B[Bi][Bpos++], c);
	      Total ++;
              ++Pos;
	    }
	}


      //-- For all the bases remaining after the last indel
      for ( i = 0; i < Remain; i ++ )
	{
	  if ( Pos >= Screen_Width ) {
            print_buffers(Buff1, Buff2, Buff3);
            add_prefix(Buff1, Apos, SeqLenR, frameR);
            add_prefix(Buff2, Bpos, SeqLenQ, frameQ);
            Pos = LINE_PREFIX_LEN;
          }

	  if ( DATA_TYPE == NUCMER_DATA )
	    c = A[Ai][Apos] == B[Bi][Bpos] ?
	      NUCMER_MATCH_CHAR : NUCMER_MISMATCH_CHAR;
	  else if ( A[Ai][Apos] == B[Bi][Bpos] )
	    c = A[Ai][Apos];
	  else
	    c = mummer::sw_align::MATCH_SCORE
	      [MATRIX_TYPE]
	      [toupper(A[Ai][Apos]) - 'A']
	      [toupper(B[Bi][Bpos]) - 'A'] > 0 ?
	      PROMER_SIM_CHAR : PROMER_MISMATCH_CHAR;
          append(Buff1, Buff2, Buff3, A[Ai][Apos++], B[Bi][Bpos++], c);
          ++Pos;
	}


      //-- For the remaining buffered output
      if ( Pos > LINE_PREFIX_LEN ) {
        print_buffers(Buff1, Buff2, Buff3);
        add_prefix(Buff1, Apos, SeqLenR, frameR);
        add_prefix(Buff2, Bpos, SeqLenQ, frameQ);
        Pos = LINE_PREFIX_LEN;
      }

      printf("\n--   END alignment [ %s%d %ld - %ld | %s%d %ld - %ld ]\n",
	     frameR > 0 ? "+" : "-", abs(frameR), Ap->sR, Ap->eR,
	     frameQ > 0 ? "+" : "-", abs(frameQ), Ap->sQ, Ap->eQ);
    }

  //-- Free the sequences, except for the originals
  for ( i = 0; i < 7; i ++ )
    {
      if ( (DATA_TYPE != NUCMER_DATA || i != 1)  &&  A[i] != NULL )
	free ( (void*)A[i] );
      if ( (DATA_TYPE != NUCMER_DATA || i != 1)  &&  B[i] != NULL )
	free ( (void*)B[i] );
    }

  return;
}

void base_color(ColoredBuffer& b, char c, bool m) {
  // std::cerr << "base_color " << c << ',' << m << '\n';
  switch(c) {
  case 'a': case 'A': b.color(ColoredBuffer::Green, m); break;
  case 'c': case 'C': b.color(ColoredBuffer::Blue, m);  break;
  case 'g': case 'G': b.color(ColoredBuffer::Black, m); break;
  case 't': case 'T': b.color(ColoredBuffer::Red, m);   break;
  default: b.reset(); break;
  }
  b += c;
}

void append(ColoredBuffer& Buff1, ColoredBuffer& Buff2, std::string &Buff3,
            char c1, char c2, char c3) {
  if(!Colorize || DATA_TYPE == PROMER_DATA) {
    Buff1 += c1;
    Buff2 += c2;
  } else {
    const bool mismatch =
      (DATA_TYPE == NUCMER_DATA && c3 == NUCMER_MISMATCH_CHAR) ||
      (DATA_TYPE == PROMER_DATA && c3 == PROMER_MISMATCH_CHAR);
    // std::cerr << "append " << c3 << ',' << NUCMER_MISMATCH_CHAR << ',' << PROMER_MISMATCH_CHAR << ' ' << mismatch << '\n';
    base_color(Buff1, c1, mismatch);
    base_color(Buff2, c2, mismatch);
    if(Buff3.empty())
      Buff3 += ANSI_BOLD;
  }
  Buff3 += c3;
}

void add_prefix(ColoredBuffer& Buff, long int pos, long int seqlen, int frame) {
  char b[LINE_PREFIX_LEN + 1];
  sprintf(b, PREFIX_FORMAT, toFwd(pos, seqlen, frame));
  Buff.clear();
  Buff += b;
}

void print_buffers(ColoredBuffer& Buff1, ColoredBuffer& Buff2, std::string& Buff3) {
  print_markers(Buff1.bases());
  if(Colorize) { // Make sure that we reset any color/property setting
    Buff1.reset();
    Buff2.reset();
    Buff3 += ANSI_RESET;
  }
  const char* b2 = Buff2.c_str();
  const char* b3 = Buff3.c_str();
  if(DATA_TYPE != NUCMER_DATA)
    std::swap(b2, b3);
  printf("%s\n"
         "%s\n"
         "%*s%s\n",
         Buff1.c_str(),
         Buff2.c_str(),
         LINE_PREFIX_LEN, "", Buff3.c_str());
  Buff1.clear(); Buff2.clear(); Buff3.clear();
}

void print_markers(int max_len) {
  static const int maximums[7] = {1, 10, 100, 1000, 10000, 100000, 1000000};
  if(Marker_Width <= 0) {
    printf("\n");
    return;
  }

  const int max = maximums[std::min(6, Marker_Width - 1)];
  printf("%*s", LINE_PREFIX_LEN, "");
  if(Colorize)
    printf("%s", ANSI_UNDERLINE);
  int i = Marker_Width;
  for( ; i <= max_len - LINE_PREFIX_LEN; i += Marker_Width) {
    if(i < max)
      printf("%*d|", Marker_Width - 1, i);
    else
      printf("%*s|", Marker_Width - 1, "");
  }
  if(Colorize) {
    i -= Marker_Width;
    if(i < max_len - LINE_PREFIX_LEN)
      printf("%*s", max_len - LINE_PREFIX_LEN - i, "");
    printf("%s", ANSI_RESET);
  }
  printf("\n");
}


void printHelp
     (const char * s)

      //  Display the program's help information to stderr

{
  fprintf (stderr,
           "\nUSAGE: %s  [options]  <deltafile>  <ref ID>  <qry ID>\n\n", s);
  fprintf (stderr,
           "-h            Display help information\n"
           "-q            Sort alignments by the query start coordinate\n"
           "-r            Sort alignments by the reference start coordinate\n"
           "-w int        Set the screen width - default is terminal width\n"
           "-c            Colorize bases on output\n"
           "-x int        Set the matrix type - default is 2 (BLOSUM 62)\n"
           "-m int        Space between markers - default is 10, disable with 0\n"
           "              other options include 1 (BLOSUM 45) and 3 (BLOSUM 80)\n"
           "              note: only has effect on amino acid alignments\n\n");
  fprintf (stderr,
           "  Input is the .delta output of either the \"nucmer\" or the\n"
           "\"promer\" program passed on the command line.\n"
           "  Output is to stdout, and consists of all the alignments between the\n"
           "query and reference sequences identified on the command line.\n"
           "  NOTE: No sorting is done by default, therefore the alignments\n"
           "will be ordered as found in the <deltafile> input.\n\n");
  return;
}




void printUsage
     (const char * s)

      //  Display the program's usage information to stderr.

{
  fprintf (stderr,
           "\nUSAGE: %s  [options]  <deltafile>  <ref ID>  <qry ID>\n\n", s);
  fprintf (stderr, "Try '%s -h' for more information.\n", s);
  return;
}




long int revC
     (long int coord, long int len)
{
  return len - coord + 1;
}

bool find_sequence(const std::vector<string>& paths, const std::string& Id, std::string& seq)
{
  //-- Find, and read in sequences. Return if find one with name Id, and store it in seq.
  stream_manager streams(paths.cbegin(), paths.cend());
  sequence_parser parser(16, 10, 1, streams);
  bool found = false;
  while(!found) {
    sequence_parser::job j(parser);
    if(j.is_empty()) break;
    for(size_t i = 0; i < j->nb_filled; ++i) {
      // Compare up to first white space or tab
      const auto n = j->data[i].header.find_first_of(" \t");
      if(j->data[i].header.compare(0, n, Id) == 0) {
        found = true;
        seq = std::string(1, '\0') + j->data[i].seq; // seq is 1-based
        for(auto& base : seq)
          base = tolower(base);
        break;
      }
    }
  }
  return found;
}
