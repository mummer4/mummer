//------------------------------------------------------------------------------
//   Programmer: Adam M Phillippy, The Institute for Genomic Research
//         File: show-tiling.cc
//         Date: 09 / 13 / 2002
//
//        Usage: show-tiling [options] <deltafile>
//               Try 'show-tiling -h' for more information
//
//  Description: For use in conjunction with the MUMmer package.
//              "show-tiling" predicts and displays a query to reference tiling
//             path. Note that this program is most effective if the query
//            sequences are small fragments of a large refence.
//
//------------------------------------------------------------------------------

#include <mummer/delta.hh>
#include <mummer/tigrinc.hh>
#include <mummer/redirect_to_pager.hpp>
#include <vector>
#include <algorithm>
using namespace std;

//------------------------------------------------------------- Constants ----//
const int CHARS_PER_LINE = 60;

const long int DEFAULT_MIN_CONTIG_LENGTH        =     1;

const long int DEFAULT_NUCMER_MAX_GAP_SIZE      =  1000;
const long int DEFAULT_PROMER_MAX_GAP_SIZE      =    -1;

const float DEFAULT_NUCMER_MIN_COVERAGE         =  95.0;
const float DEFAULT_PROMER_MIN_COVERAGE         =  50.0;

const float DEFAULT_NUCMER_MIN_COVERAGE_DIFF    =  10.0;
const float DEFAULT_PROMER_MIN_COVERAGE_DIFF    =  30.0;

const float DEFAULT_NUCMER_MIN_PIDY             =  90.0;
const float DEFAULT_PROMER_MIN_PIDY             =  55.0;

const char FORWARD_CHAR  =  '+';
const char REVERSE_CHAR  =  '-';

const float REPEAT_PIDY_DIFF = 0.25;

const int UNUSABLE_TILE_LEVEL  =  -2;
const int IGNORE_TILE_LEVEL    =  -1;
const int UNUSED_TILE_LEVEL    =   0;
const int USED_TILE_LEVEL      =   1;

char NULL_STRING[1] = "";



//------------------------------------------------------ Type Definitions ----//
struct AlignStats
     //-- Alignment statistics data structure
{
  float Idy;                            // percent identity   [0.0, 100.0]
  float Sim;                            // percent similarity [0.0, 100.0]
  float Stp;                            // precent stop codon [0.0, 100.0]

  char DirQ;                            // contig orientation (relative to ref)

  long int sQ, eQ, sR, eR;              // start and end in Query and Reference
                                        // relative to the directional strand

  char * IdR;                           // FASTA Id of the mapping reference
  long int SeqLenR;                     // length of the reference

  bool isTiled;                          // is the alignment be tiled?
};



struct QueryContig
     //-- Query sequence tiling contig data structure
{
  char * IdQ;                              // FASTA Id of the query
  long int SeqLenQ;                        // length of the query
  char * SeqQ;                             // query sequence

  vector<AlignStats> Aligns;               // alignments for this contig

  int TileLevel;                      // describes tiling status of query contig

  //-- Things to be filled in once the contig is tiled
  char DirQ;                               // orientation of the contig
  char * IdR;                              // Id of the mapping reference
  long int SeqLenR;                        // sequence length of the reference
  long int StartR, EndR;                   // contig -> reference mapping coords

  long int LoTrim, HiTrim;                 // lo and hi trim lengths

  vector<QueryContig>::iterator linksTo;   // who does this contig link to?
  bool isLinkHead;                         // is the head of the linked list?
};



struct IdR_StartR_Sort
//-- For sorting QueryContigs by IdR, then StartR coordinates
{
  bool operator( ) (const QueryContig & pA, const QueryContig & pB)
  {
    int cmp = strcmp (pA.IdR, pB.IdR);

    //-- sort IdR
    if ( cmp < 0 )
      return true;
    else if ( cmp > 0 )
      return false;

    //-- sort StartR
    else if ( pA.StartR < pB.StartR )
      return true;
    else if ( pA.StartR > pB.StartR )
      return false;

    //-- sort SeqLenQ
    else if ( pA.SeqLenQ > pB.SeqLenQ )
      return true;
    else
      return false;
  }
};



struct IdR_StartRTrimmed_Sort
//-- For sorting contig output by trim offset
{
  bool operator( ) (const QueryContig & pA, const QueryContig & pB)
  {
    int cmp = strcmp (pA.IdR, pB.IdR);

    //-- sort IdR
    if ( cmp < 0 )
      return true;
    else if ( cmp > 0 )
      return false;

    //-- sort StartR
    else if ( pA.StartR + pA.LoTrim < pB.StartR + pA.LoTrim )
      return true;
    else
      return false;
  }
};



struct IdQ_Sort
//-- For sorting QueryContigs by IdQ
{
  bool operator( ) (const QueryContig & pA, const QueryContig & pB)
  {
    //-- sort IdQ
    if ( strcmp (pA.IdQ, pB.IdQ) < 0 )
      return true;
    else
      return false;
  }
};



struct IdR_sQ_Sort
//-- For sorting AlignStats by IdR, then sQ
{
  bool operator( ) (const AlignStats & pA, const AlignStats & pB)
  {
    int cmp = strcmp (pA.IdR, pB.IdR);

    //-- sort IdR
    if ( cmp < 0 )
      return true;
    else if ( cmp > 0 )
      return false;

    //-- sort sQ
    else if ( pA.sQ < pB.sQ )
      return true;
    else
      return false;
  }
};




//------------------------------------------------------ Global Variables ----//
bool isOutputContigs = false;            // set by -t option
bool isOutputPseudoMolecule = false;     // set by -p option
bool isOutputUnusable = false;           // set by -u option
bool isPrintAlignments = false;          // set by -a option
bool isPrintXML = false;                 // set by -x option
bool isCircularReference = false;        // set by -c option
bool isRandomRepeats = false;            // set by -R option

long int MIN_CONTIG_LENGTH = DEFAULT_MIN_CONTIG_LENGTH;

long int MAX_GAP_SIZE    = 0;
float MIN_COVERAGE       = 0;
float MIN_COVERAGE_DIFF  = 0;
float MIN_PIDY           = 0;

bool isdef_MAX_GAP_SIZE       = false;
bool isdef_MIN_COVERAGE       = false;
bool isdef_MIN_COVERAGE_DIFF  = false;
bool isdef_MIN_PIDY           = false;

int   DATA_TYPE = NUCMER_DATA;         // set by .delta header

char InputFileName [MAX_LINE];
char RefFileName [MAX_LINE], QryFileName [MAX_LINE];



//------------------------------------------------- Function Declarations ----//
float getAlignmentQryCoverage
     (vector<AlignStats>::iterator Ap, long int SeqLenQ);

float getSubsetQryCoverage
     (vector<AlignStats>::iterator begin,
      vector<AlignStats>::iterator end, long int SeqLenQ);

float getSubsetQrySyntenyCoverage
     (vector<AlignStats>::iterator begin,
      vector<AlignStats>::iterator end, long int SeqLenQ);

float getSubsetIdentity
     (vector<AlignStats>::iterator begin,
      vector<AlignStats>::iterator end);

void linkContigs
     (vector<QueryContig> & Contigs);

long int longestConsistentSubset
     (vector<AlignStats>::iterator begin,
      vector<AlignStats>::iterator end);

void outputContigs
     (vector<QueryContig> Contigs, FILE * Output);

void outputPseudoMolecule
     (vector<QueryContig> & Contigs, FILE * QryFile, FILE * Output);

void outputUnusable
     (vector<QueryContig> Contigs, FILE * Output);

void parseDelta
     (vector<QueryContig> & Contigs);

void placeContig
     (vector<QueryContig>::iterator Cp);

void printAlignment
     (vector<QueryContig>::iterator Cp, vector<AlignStats>::iterator Ap,
      FILE * Output);

void printTilingAlignments
     (vector<QueryContig> & Contigs);

void printTilingPath
     (vector<QueryContig> & Contigs);

void printTilingXML
     (vector<QueryContig> & Contigs, char * QryFileName,
      int argc, char ** argv);

void printHelp
     (const char * s);

void printUsage
     (const char * s);

inline long int revC
     (long int Coord, long int Len);

void tileContigs
     (vector<QueryContig> & Contigs);




//-------------------------------------------------- Function Definitions ----//
int main
     (int argc, char ** argv)
{
  char ContigsFileName [MAX_LINE];
  char PseudoMoleculeFileName [MAX_LINE];
  char UnusableFileName [MAX_LINE];

  FILE * ContigsFile = NULL;
  FILE * PseudoMoleculeFile = NULL;
  FILE * UnusableFile = NULL;
  FILE * QryFile = NULL;

  vector<QueryContig> Contigs;

  //-- Parse the command line arguments
  {
    int ch, errflg = 0;
    optarg = NULL;

    while ( !errflg  &&  ((ch = getopt
			   (argc, argv, "achg:i:l:p:Rt:u:v:V:x")) != EOF) )
      switch (ch)
        {
	case 'a' :
	  isPrintAlignments = true;
	  break;

	case 'c' :
	  isCircularReference = true;
	  break;

        case 'h' :
          printHelp (argv[0]);
          exit (EXIT_SUCCESS);
          break;

	case 'g' :
	  MAX_GAP_SIZE = atoi (optarg);
	  isdef_MAX_GAP_SIZE = true;
	  break;

	case 'i' :
	  MIN_PIDY = atof (optarg);
	  isdef_MIN_PIDY = true;
	  break;

	case 'l' :
	  MIN_CONTIG_LENGTH = atoi (optarg);
	  break;

	case 'p' :
	  strcpy ( PseudoMoleculeFileName, optarg );
	  isOutputPseudoMolecule = true;
	  break;

	case 'R' :
	  isRandomRepeats = true;
	  break;

	case 't' :
	  strcpy ( ContigsFileName, optarg );
	  isOutputContigs = true;
	  break;

	case 'u' :
	  strcpy ( UnusableFileName, optarg );
	  isOutputUnusable = true;
	  break;

	case 'v' :
	  MIN_COVERAGE = atof (optarg);
	  isdef_MIN_COVERAGE = true;
	  break;

	case 'V' :
	  MIN_COVERAGE_DIFF = atof (optarg);
	  isdef_MIN_COVERAGE_DIFF = true;
	  break;

	case 'x' :
	  isPrintXML = true;
	  break;

        default :
          errflg ++;
        }

    if ( errflg > 0  ||  argc - optind != 1 )
      {
        printUsage (argv[0]);
        exit (EXIT_FAILURE);
      }

    if ( MIN_PIDY < 0.0 || MIN_PIDY > 100.0 ||
	 MIN_COVERAGE < 0.0 || MIN_COVERAGE > 100.0 ||
	 MIN_COVERAGE_DIFF < 0.0 || MIN_COVERAGE_DIFF > 100.0 )
      {
	fprintf(stderr,
		"\nERROR: Percentages must be within the range [0.0, 100.0]\n");
	exit (EXIT_FAILURE);
      }

    if ( MIN_CONTIG_LENGTH < -1 || MAX_GAP_SIZE < -1 )
      {
	fprintf(stderr,
		"\nERROR: Size values must be within the range [-1, oo]\n");
	exit (EXIT_FAILURE);
      }

    if ( isRandomRepeats )
      {
	srand ( time(NULL) );
	MIN_COVERAGE_DIFF = 0;
	isdef_MIN_COVERAGE_DIFF = true;
      }
  }


  //-- Parse the delta file
  strcpy (InputFileName, argv[optind ++]);
  parseDelta (Contigs);


  //-- Try and open all the output files
  if ( isOutputContigs )
    ContigsFile = File_Open ( ContigsFileName, "w" );
  if ( isOutputPseudoMolecule )
    PseudoMoleculeFile = File_Open ( PseudoMoleculeFileName, "w" );
  if ( isOutputPseudoMolecule )
    QryFile = File_Open ( QryFileName, "r" );
  if ( isOutputUnusable )
    UnusableFile = File_Open ( UnusableFileName, "w" );


  //-- Find each contig's mapping and sort by reference
  tileContigs (Contigs);


  if ( isOutputContigs )
    {
      outputContigs (Contigs, ContigsFile);
      fclose (ContigsFile);
    }


  //-- Print the unusable contigs w/ alignments to the user specified file
  if ( isOutputUnusable )
    {
      outputUnusable (Contigs, UnusableFile);
      fclose (UnusableFile);
    }


  //-- Link the contigs to form the tiling path
  linkContigs (Contigs);


  //-- Print the pseudo molecule to the user specifed file
  if ( isOutputPseudoMolecule )
    {
      outputPseudoMolecule (Contigs, QryFile, PseudoMoleculeFile);
      fclose (QryFile);
      fclose (PseudoMoleculeFile);
    }


  //-- Print the tiling path to stdout
  stdio_launch_pager redirect_to_pager;
  if ( isPrintAlignments )
    printTilingAlignments (Contigs);
  else if ( isPrintXML )
    printTilingXML (Contigs, QryFileName, argc, argv);
  else
    printTilingPath (Contigs);


  return EXIT_SUCCESS;
}




float getAlignmentQryCoverage
     (vector<AlignStats>::iterator Ap, long int SeqLenQ)
{
  //-- Return query coverage for a single alignment
  return (float)(Ap->eQ - Ap->sQ + 1) / (float)SeqLenQ * 100.0;
}




float getSubsetQryCoverage
     (vector<AlignStats>::iterator begin,
      vector<AlignStats>::iterator end, long int SeqLenQ)
{
  vector<AlignStats>::iterator Ap;
  vector<AlignStats>::iterator preAp = end;

  long int cov = 0;
  long int olap;

  //-- Get non-redundant query coverage for the subset of alignments
  for ( Ap = begin; Ap < end; Ap ++ )
    {
      if ( ! Ap->isTiled )
	continue;
      
      if ( preAp == end )
	cov += Ap->eQ - Ap->sQ + 1;
      else
	{
	  //-- Alignments must be pre-sorted by sQ
	  assert ( Ap->sQ >= preAp->sQ );
	  olap = preAp->eQ - Ap->sQ + 1;
	  olap = olap > 0 ? olap : 0;
	  cov += Ap->eQ - Ap->sQ + 1 - olap;
	}
      preAp = Ap;
    }

  return (float)(cov) / (float)SeqLenQ * 100.0;
}




float getSubsetQrySyntenyCoverage
     (vector<AlignStats>::iterator begin,
      vector<AlignStats>::iterator end, long int SeqLenQ)
{
  vector<AlignStats>::iterator Ap;
  long int cov = 0;
  long int low, high;

  low = LONG_MAX;
  high = -LONG_MAX;
  for ( Ap = begin; Ap < end; Ap ++ )
    if ( Ap->isTiled )
      {
	if ( Ap->sQ < low )
	  low = Ap->sQ;
	if ( Ap->eQ > high )
	  high = Ap->eQ;
      }
  cov = high - low + 1;

  return (float)(cov) / (float)SeqLenQ * 100.0;
}




float getSubsetIdentity
     (vector<AlignStats>::iterator begin,
      vector<AlignStats>::iterator end)
{
  vector<AlignStats>::iterator Ap;

  long int len;
  long int N = 0;
  float tot = 0;

  //-- Get weighted (by length) average identity for the subset of alignments
  for ( Ap = begin; Ap < end; Ap ++ )
    {
      if ( ! Ap->isTiled )
	continue;

      len = Ap->eQ - Ap->sQ + 1;
      N += len;
      tot += Ap->Idy * len;
    }

  return tot / (float)N;
}




void linkContigs
     (vector<QueryContig> & Contigs)
{
  vector<QueryContig>::iterator Cp;
  vector<QueryContig>::iterator Cip;
  vector<QueryContig>::iterator nxtCp;
  vector<QueryContig>::iterator firstCp;
  long int best_end;

  //-- Ignore shadowed contigs and generate the tiling path that uses
  //   the minimum number of contigs while covering the maximum area
  firstCp = Contigs.end( );
  for ( Cp = Contigs.begin( ); Cp < Contigs.end( ); Cp ++ )
    {
      if ( Cp->TileLevel != USED_TILE_LEVEL )
	continue;

      if ( firstCp == Contigs.end( ) )
	{
	  firstCp = Cp;
	  firstCp->isLinkHead = true;
	}

      //-- If several contigs overlap Cp, pick the best one...
      //   or if none overlap, pick the next used tile
      nxtCp = Contigs.end( );
      best_end = -(LONG_MAX);
      for ( Cip = Cp + 1; Cip < Contigs.end( ); Cip ++ )
	{
	  if ( Cip->TileLevel != USED_TILE_LEVEL )
	    continue;

	  if ( strcmp ( Cip->IdR, Cp->IdR ) != 0 )
	    break;

	  //-- If no overlap
	  if ( Cip->StartR > Cp->EndR )
	    {
	      if ( nxtCp == Contigs.end( ) )
		nxtCp = Cip;
	      break;
	    }

	  //-- If overlap
	  Cip->TileLevel = IGNORE_TILE_LEVEL;
	  if ( Cip->EndR > Cp->EndR  &&  Cip->EndR > best_end )
	    {
	      best_end = Cip->EndR;
	      nxtCp = Cip;
	    }
	}


      //-- Link Cp to nxtCp
      if ( nxtCp != Contigs.end( ) )
	{
	  //-- Link to the next valid tile
	  nxtCp->TileLevel = USED_TILE_LEVEL;
	  Cp->linksTo = nxtCp;
	}
      else
	{
	  if ( isCircularReference )
	    {
	      if ( firstCp->StartR <= 0 )
		{
		  nxtCp = Cp;
		  for (; Cp > firstCp; Cp -- )
		    if ( Cp->TileLevel == USED_TILE_LEVEL )
		      {
			if ( firstCp->SeqLenR + firstCp->StartR <= Cp->EndR )
			  {
			    nxtCp = Cp;
			    Cp->TileLevel = IGNORE_TILE_LEVEL;
			  }
			else
			  break;
		      }

		  if ( firstCp->SeqLenR + firstCp->StartR > nxtCp->StartR )
		    {
		      Cp = nxtCp;
		      Cp->TileLevel = USED_TILE_LEVEL;
		    }
		}

	      //-- Circular link (note: this creates a circular linked list)
	      Cp->linksTo = firstCp;
	    }
	  else
	    Cp->linksTo = Contigs.end( );

	  firstCp = Contigs.end( );
	}
    }

  return;
}




long int longestConsistentSubset
     (vector<AlignStats>::iterator begin, vector<AlignStats>::iterator end)
{
  //-- For dynamic programming scoring and backtracking
  struct node
  {
    vector<AlignStats>::iterator Ap;
    long int Score;
    long int From;
  };

  vector<AlignStats>::iterator Ap;
  long int i, j, N, OlapR, OlapQ, Olap, Best;
  node * A;

  //-- Build the data array
  N = 0;
  A = (node *) Safe_malloc (sizeof(node) * (end - begin));
  for ( Ap = begin; Ap < end; Ap ++ )
    {
      assert ( Ap->isTiled == false );
      A[N].Ap = Ap;
      A[N].Score = Ap->eQ - Ap->sQ + 1;
      A[N].From = -1;
      N ++;
    }
  assert ( N == (end - begin) );

  //-- Isn't it dynamic?
  for ( i = 0; i < N; i ++ )
    for ( j = 0; j < i; j ++ )
      {
	assert ( strcmp (A[i].Ap->IdR, A[j].Ap->IdR) == 0 );
	if ( A[i].Ap->DirQ != A[j].Ap->DirQ )
	  continue;

	OlapR = A[j].Ap->eR - A[i].Ap->sR + 1;
	if ( isCircularReference  &&  OlapR > A[j].Ap->SeqLenR / 2 )
	  OlapR = (A[j].Ap->eR - A[j].Ap->SeqLenR) + (1 - A[i].Ap->sR);

	if ( MAX_GAP_SIZE >= 0  &&  OlapR < -(MAX_GAP_SIZE) )
	  continue;
	Olap = OlapR < 0 ? 0 : OlapR;
	OlapQ = A[j].Ap->eQ - A[i].Ap->sQ + 1;
	if ( MAX_GAP_SIZE >= 0  &&  OlapQ < -(MAX_GAP_SIZE) )
	  continue;
	Olap = Olap > OlapQ ? Olap : OlapQ;

	if ( A[j].Score + (A[i].Ap->eQ - A[i].Ap->sQ + 1) - Olap > A[i].Score )
	  {
	    A[i].From = j;
	    A[i].Score = A[j].Score + (A[i].Ap->eQ - A[i].Ap->sQ + 1) - Olap;
	  }
      }

  //-- Find the highest score and backtrack to extract the subset
  Best = 0;
  for ( i = 1; i < N; i ++ )
    {
      if ( A[i].Score > A[Best].Score )
	Best = i;
      else if ( A[i].Score == A[Best].Score &&
		A[i].Ap->Idy > A[Best].Ap->Idy )
	Best = i;
    }

  if ( isRandomRepeats )
    {
      vector<long int> ties;
      for ( i = 0; i < N; i ++ )
	if ( A[i].Score == A[Best].Score &&
	     A[i].Ap->Idy >= A[Best].Ap->Idy - REPEAT_PIDY_DIFF )
	  ties . push_back (i);

      if ( ties.size( ) > 1 )
	Best = ties
	  [(long int)((double)ties.size( ) * rand( ) / (RAND_MAX + 1.0))];
    }

  for ( i = Best; i >= 0; i = A[i].From )
    A[i].Ap->isTiled = true;

  Best = A[Best].Score;

  free ( A );

  return Best;
}




void outputContigs
     (vector<QueryContig> Contigs, FILE * Output)
{
  vector<QueryContig>::iterator Cp;
  vector<QueryContig>::iterator endCp;
 
  //-- Sort by reference trim offset
  sort ( Contigs.begin( ), Contigs.end( ), IdR_StartRTrimmed_Sort( ) );

  long int seqs;

  Cp = Contigs.begin( );
  while ( Cp < Contigs.end( ) )
    {
      if ( Cp->TileLevel != USED_TILE_LEVEL )
	{
	  Cp ++;
	  continue;
	}

      seqs = 0;
      for ( endCp = Cp; endCp < Contigs.end( ); endCp ++ )
	{
	  if ( endCp->TileLevel != USED_TILE_LEVEL )
	    continue;

	  if ( strcmp (Cp->IdR, endCp->IdR) != 0 )
	    break;
	  seqs ++;
	}
	  
      fprintf(Output,
      "##%s %ld %ld bases, 00000000 checksum.\n",
	      Cp->IdR, seqs, Cp->SeqLenR);

      for ( ; Cp < endCp; Cp ++ )
	{
	  if ( Cp->TileLevel != USED_TILE_LEVEL )
	    continue;

	  if ( Cp->DirQ == FORWARD_CHAR )
	    fprintf(Output,
            "#%s(%ld) [%s] %ld bases, 00000000 checksum. {%ld %ld} <%ld %ld>\n",
		    Cp->IdQ, Cp->StartR + Cp->LoTrim - 1, "", Cp->SeqLenQ,
		    Cp->LoTrim + 1, Cp->SeqLenQ - Cp->HiTrim,
		    Cp->StartR + Cp->LoTrim, Cp->EndR - Cp->HiTrim);
	  else
	    fprintf(Output,
            "#%s(%ld) [%s] %ld bases, 00000000 checksum. {%ld %ld} <%ld %ld>\n",
		    Cp->IdQ, Cp->StartR + Cp->LoTrim - 1, "RC", Cp->SeqLenQ,
		    Cp->SeqLenQ - Cp->LoTrim, Cp->HiTrim + 1,
		    Cp->StartR + Cp->LoTrim, Cp->EndR - Cp->HiTrim);
	}
    }

  return;
}




void outputPseudoMolecule
     (vector<QueryContig> & Contigs, FILE * QryFile, FILE * Output)
{
  vector<QueryContig>::iterator Cp;
  vector<QueryContig>::iterator beginCp;

  vector<AlignStats>::iterator Apb;
  vector<AlignStats>::reverse_iterator Ape;

  long int ct = 0;
  long int start, end, gap, i, endR, startR;

  char * A;
  long int InitSize = INIT_SIZE;
  char Line [MAX_LINE];


  //-- Read in the needed query contig sequences
  A = (char *) Safe_malloc ( sizeof(char) * InitSize );
  while ( Read_String (QryFile, A, InitSize, Line, false) )
    {
      for ( Cp = Contigs.begin( ); Cp < Contigs.end( ); Cp ++ )
	if ( Cp->TileLevel == USED_TILE_LEVEL )
	  if ( strcmp ( Line, Cp->IdQ ) == 0 )
	    break;
      
      if ( Cp < Contigs.end( ) )
	{
	  assert ( (long int)strlen(A+1) == Cp->SeqLenQ );
	  Cp->SeqQ = (char *) Safe_malloc
	    ( sizeof(char) * (Cp->SeqLenQ + 2) );
	  Cp->SeqQ[0] = '\0';
	  strcpy ( Cp->SeqQ + 1, A + 1 );
	  if ( Cp->DirQ == REVERSE_CHAR )
	    Reverse_Complement (Cp->SeqQ, 1, Cp->SeqLenQ);
	}
    }
  free ( A );

  //-- For all contigs, create pseudo
  for ( beginCp = Contigs.begin( ); beginCp < Contigs.end( ); beginCp ++ )
    {
      if ( ! beginCp->isLinkHead )
	continue;

      //-- New Reference sequence
      start = 1;
      Cp = beginCp;
      fprintf (Output, ">pseudo_used_%s\n", Cp->IdR);

      //-- For all contigs mapping to this Reference
      while ( Cp != Contigs.end( ) )
	{
	  gap = 0;
	  end = Cp->SeqLenQ;
	  if ( Cp->linksTo != Contigs.end( ) )
	    {
	      if ( ! Cp->linksTo->isLinkHead )
		{
		  //-- Internal sequential link, gap N's needed if necessary
		  startR = Cp->linksTo->StartR;
		  endR = Cp->EndR;

		  gap = startR - endR - 1;
		  if ( gap < 0 )
		    {
		      for ( Apb = Cp->linksTo->Aligns.begin( );
			    Apb < Cp->linksTo->Aligns.end( )  &&  !Apb->isTiled;
			    Apb ++ )
			{}
		      assert ( Apb < Cp->linksTo->Aligns.end( ) );

		      for ( Ape = Cp->Aligns.rbegin( );
			    Ape < Cp->Aligns.rend( )  &&  !Ape->isTiled;
			    Ape ++ )
			{}
		      assert ( Ape < Cp->Aligns.rend( ) );
		      
		      //-- Use the contig with the least 'junk' sequence
		      if ( Cp->SeqLenQ - Ape->eQ > Apb->sQ - 1  ||
			   ( Cp->SeqLenQ - Ape->eQ == Apb->sQ - 1  &&
			     getSubsetIdentity (Cp->Aligns.begin( ),
						Cp->Aligns.end( ))  <
			     getSubsetIdentity (Cp->linksTo->Aligns.begin( ),
						Cp->linksTo->Aligns.end( )) ) )
			end += gap;
		    }
		}
	    }

	  if ( Cp->SeqQ == NULL )
	    {
	      fprintf (stderr,
		 "\nERROR: Sequence \"%s\" was not found in the query file.\n"
		 "       Please check the validity of the query file listed\n"
		 "       at the top of the .delta input file and rerun.\n",
		       Cp->IdQ);
	      exit (EXIT_FAILURE);
	    }

	  //-- Print the sequence
	  for ( i = start; i <= end; i ++ )
	    {
	      fputc (toupper(Cp->SeqQ[i]), Output);
	      if ( ++ ct == CHARS_PER_LINE )
		{
		  ct = 0;
		  fputc ('\n', Output);
		}
	    }
	  free (Cp->SeqQ);

	  //-- Print the gap
	  for ( i = 1; i <= gap; i ++ )
	    {
	      fputc ('N', Output);
	      if ( ++ ct == CHARS_PER_LINE )
		{
		  ct = 0;
		  fputc ('\n', Output);
		}
	    }

	  //-- Adjust the next start if this contig was used for overlap
	  if ( gap < 0  &&  end == Cp->SeqLenQ )
	    start = 1 + -(gap);

	  //-- Walk the pointer down the list of contigs
	  if ( Cp->linksTo == Contigs.end( )  ||  Cp->linksTo->isLinkHead )
	    Cp = Contigs.end( );
	  else
	    Cp = Cp->linksTo;
	}

      if ( ct != 0 )
	fputc ('\n', Output);
      ct = 0;
    }

  return;
}




void outputUnusable
     (vector<QueryContig> Contigs, FILE * Output)
{
  vector<QueryContig>::iterator Cp;
  vector<AlignStats>::iterator Ap;

  sort ( Contigs.begin( ), Contigs.end( ), IdQ_Sort( ) );

  for ( Cp = Contigs.begin( ); Cp < Contigs.end( ); Cp ++ )
    {
      if ( Cp->TileLevel != UNUSABLE_TILE_LEVEL  &&
	   Cp->TileLevel != UNUSED_TILE_LEVEL )
	continue;

      for ( Ap = Cp->Aligns.begin( ); Ap < Cp->Aligns.end( ); Ap ++ )
	printAlignment (Cp, Ap, Output);
    }
  
  return;
}




void parseDelta
     (vector<QueryContig> & Contigs)

     //  parse the delta file

{
  vector<QueryContig>::iterator Cp;

  char * CurrIdQ;                        // the current contig Id
  long int temp;

  QueryContig aContig;                   //  single query contig
  AlignStats aStats;                     //  single alignment region

  DeltaReader_t dr;
  dr.open (InputFileName);
  DATA_TYPE = dr.getDataType( ) == NUCMER_STRING ?
    NUCMER_DATA : PROMER_DATA;
  strcpy (RefFileName, dr.getReferencePath( ).c_str( ));
  strcpy (QryFileName, dr.getQueryPath( ).c_str( ));


  Contigs.clear( );
  CurrIdQ = NULL_STRING;
  aContig.SeqQ = NULL;
  aContig.DirQ = '*';
  aContig.IdR = NULL_STRING;
  aContig.SeqLenR = 0;
  aContig.Aligns.clear( );
  aContig.StartR = aContig.EndR = 0;
  aContig.LoTrim = aContig.HiTrim = 0;
  aContig.isLinkHead = false;


  if ( DATA_TYPE == NUCMER_DATA )
    {
      if ( !isdef_MAX_GAP_SIZE )
	MAX_GAP_SIZE = DEFAULT_NUCMER_MAX_GAP_SIZE;
      if ( !isdef_MIN_COVERAGE )
	MIN_COVERAGE = DEFAULT_NUCMER_MIN_COVERAGE;
      if ( !isdef_MIN_COVERAGE_DIFF )
	MIN_COVERAGE_DIFF = DEFAULT_NUCMER_MIN_COVERAGE_DIFF;
      if ( !isdef_MIN_PIDY )
	MIN_PIDY = DEFAULT_NUCMER_MIN_PIDY;
    }
  else
    {
      if ( !isdef_MAX_GAP_SIZE )
	MAX_GAP_SIZE = DEFAULT_PROMER_MAX_GAP_SIZE;
      if ( !isdef_MIN_COVERAGE )
	MIN_COVERAGE = DEFAULT_PROMER_MIN_COVERAGE;
      if ( !isdef_MIN_COVERAGE_DIFF )
	MIN_COVERAGE_DIFF = DEFAULT_PROMER_MIN_COVERAGE_DIFF;
      if ( !isdef_MIN_PIDY )
	MIN_PIDY = DEFAULT_PROMER_MIN_PIDY;
    }


  //-- Process the delta input file
  while ( dr.readNext( ) )
    {
      aStats.SeqLenR = dr.getRecord( ).lenR;
      aContig.SeqLenQ = dr.getRecord( ).lenQ;
      
      aStats.IdR = (char *) Safe_malloc (dr.getRecord( ).idR.length( ) + 1);
      aContig.IdQ = (char *) Safe_malloc (dr.getRecord( ).idQ.length( ) + 1);
      strcpy (aStats.IdR, dr.getRecord( ).idR.c_str( ));
      strcpy (aContig.IdQ, dr.getRecord( ).idQ.c_str( ));

      if ( strcmp (CurrIdQ, aContig.IdQ) )
	{
	  CurrIdQ = aContig.IdQ;
	  aContig.TileLevel = aContig.SeqLenQ < MIN_CONTIG_LENGTH ?
	    IGNORE_TILE_LEVEL : UNUSED_TILE_LEVEL;
	  Contigs.push_back (aContig);
	}

      for ( unsigned int i = 0; i < dr.getRecord( ).aligns.size( ); i ++ )
	{
	  aStats.sR = dr.getRecord( ).aligns[i].sR;
	  aStats.eR = dr.getRecord( ).aligns[i].eR;
	  aStats.sQ = dr.getRecord( ).aligns[i].sQ;
	  aStats.eQ = dr.getRecord( ).aligns[i].eQ;

          //-- Check match orientation
	  if ( (aStats.sR <= aStats.eR  &&  aStats.sQ <= aStats.eQ) ||
	       (aStats.sR > aStats.eR  &&  aStats.sQ > aStats.eQ) )
	    aStats.DirQ = FORWARD_CHAR;
          else
	    aStats.DirQ = REVERSE_CHAR;

	  //-- Force ascending coordinates
	  if ( aStats.sR > aStats.eR )
	    {
	      temp = aStats.sR;
	      aStats.sR = aStats.eR;
	      aStats.eR = temp;
	    }
	  if ( aStats.sQ > aStats.eQ )
	    {
	      temp = aStats.sQ;
	      aStats.sQ = aStats.eQ;
	      aStats.eQ = temp;
	    }

	  //-- If flipped orientation, reverse query coordinates
	  if ( aStats.DirQ == REVERSE_CHAR )
	    {
	      temp = aStats.sQ;
	      aStats.sQ = revC (aStats.eQ, aContig.SeqLenQ);
	      aStats.eQ = revC (temp, aContig.SeqLenQ);
	    }
	  
	  //-- Set the statistics for this alignment region
	  aStats.Idy = dr.getRecord( ).aligns[i].idy;
	  aStats.Sim = dr.getRecord( ).aligns[i].sim;
	  aStats.Stp = dr.getRecord( ).aligns[i].stp;
	  aStats.isTiled = false;

	  //-- Add the alignment region
	  Contigs.rbegin( )->Aligns.push_back (aStats);
        }
    }
  dr.close( );

  for ( Cp = Contigs.begin( ); Cp < Contigs.end( ); Cp ++ )
    sort ( Cp->Aligns.begin( ), Cp->Aligns.end( ), IdR_sQ_Sort( ) );

  return;
}




void placeContig
     (vector<QueryContig>::iterator Cp)
{
  vector<AlignStats>::iterator Aip, Ap;
  float max_cov, cov;
  long int start;

  assert ( Cp->TileLevel == USED_TILE_LEVEL );

  //-- Find the 'representative' alignment
  Ap = Cp->Aligns.end( );
  max_cov = -(FLT_MAX);
  for ( Aip = Cp->Aligns.begin( ); Aip < Cp->Aligns.end( ); Aip ++ )
    if ( Aip->isTiled )
      {
	//-- Set trim values
	if ( Ap == Cp->Aligns.end( ) )
	  Cp->LoTrim = Aip->sQ - 1;
	Cp->HiTrim = Cp->SeqLenQ - Aip->eQ;

	cov = getAlignmentQryCoverage (Aip, Cp->SeqLenQ);
	if ( cov > max_cov )
	  {
	    max_cov = cov;
	    Ap = Aip;
	  }
      }
  
  //-- Set mapping reference data
  assert ( Ap != Cp->Aligns.end( ) );
  Cp->DirQ = Ap->DirQ;
  Cp->IdR = Ap->IdR;
  Cp->SeqLenR = Ap->SeqLenR;

  //-- Position the contig
  start = Ap->sQ;
  Cp->StartR = Ap->sR - (start - 1);
  Cp->EndR = Cp->StartR + Cp->SeqLenQ - 1;

  //-- Force negative offset if circular, else just let it overhang 
  if ( isCircularReference  &&  Cp->EndR > Cp->SeqLenR )
    {
      Cp->StartR = Cp->StartR - Cp->SeqLenR;
      Cp->EndR = Cp->EndR - Cp->SeqLenR;
    }
  assert ( Cp->StartR <= Cp->EndR );
  assert ( Cp->LoTrim >= 0  &&  Cp->HiTrim >= 0 );

  return;
}




void printAlignment
     (vector<QueryContig>::iterator Cp, vector<AlignStats>::iterator Ap,
      FILE * Output)
{
  long int len1, len2;
  float covA, covB;

  len1 = Ap->eR - Ap->sR + 1;
  len2 = Ap->eQ - Ap->sQ + 1;
  covA = (float)len1 / (float)Ap->SeqLenR * 100.0;
  covB = getAlignmentQryCoverage ( Ap, Cp->SeqLenQ );

  //-- Output the statistics for this alignment region
  fprintf(Output,"%ld\t%ld\t", Ap->sR, Ap->eR);
  if ( Ap->DirQ == FORWARD_CHAR )
    fprintf(Output,"%ld\t%ld\t", Ap->sQ, Ap->eQ);
  else
    fprintf(Output,"%ld\t%ld\t", revC(Ap->sQ, Cp->SeqLenQ),
	                         revC(Ap->eQ, Cp->SeqLenQ));
  fprintf(Output,"%ld\t%ld\t", len1, len2);
  fprintf(Output,"%.2f\t", Ap->Idy);
  fprintf(Output,"%ld\t%ld\t", Ap->SeqLenR, Cp->SeqLenQ);
  fprintf(Output,"%.2f\t%.2f\t", covA, covB);
  fprintf(Output,"%s\t%s", Ap->IdR, Cp->IdQ);
  fprintf(Output,"\n"); 

  return;
}




void printTilingAlignments
     (vector<QueryContig> & Contigs)
{
  vector<QueryContig>::iterator beginCp;
  vector<QueryContig>::iterator Cp;
  vector<AlignStats>::iterator Ap;
  
  for ( beginCp = Contigs.begin( ); beginCp < Contigs.end( ); beginCp ++ )
    {
      if ( ! beginCp->isLinkHead )
	continue;

      //-- A new Reference sequence
      Cp = beginCp;

      //-- For all contigs mapping to this reference
      while ( Cp != Contigs.end( ) )
	{
	  for ( Ap = Cp->Aligns.begin( ); Ap < Cp->Aligns.end( ); Ap ++ )
	    if ( ! Ap->isTiled )
	      continue;
	    else
	      printAlignment(Cp,Ap,stdout);

	  //-- Walk the pointer down the list of contigs
	  if ( Cp->linksTo == Contigs.end( )  ||  Cp->linksTo->isLinkHead )
	    Cp = Contigs.end( );
	  else
	    Cp = Cp->linksTo;
	}
    }

  return;
}




void printTilingPath
     (vector<QueryContig> & Contigs)
{
  vector<QueryContig>::iterator beginCp;
  vector<QueryContig>::iterator Cp;

  long int len, gap;
  float pcov, pidy;


  for ( beginCp = Contigs.begin( ); beginCp < Contigs.end( ); beginCp ++ )
    {
      if ( ! beginCp->isLinkHead )
	continue;

      //-- A new Reference sequence
      Cp = beginCp;
      printf (">%s %ld bases\n", Cp->IdR, Cp->SeqLenR);

      //-- For all contigs mapping to this reference
      while ( Cp != Contigs.end( ) )
	{
	  len = Cp->SeqLenQ;
	  pcov = getSubsetQryCoverage
	    ( Cp->Aligns.begin( ), Cp->Aligns.end( ), Cp->SeqLenQ );
	  pidy = getSubsetIdentity
	    ( Cp->Aligns.begin( ), Cp->Aligns.end( ) );

	  //-- Calculate the gap size between this and the next contig
	  if ( Cp->linksTo == Contigs.end( ) )
	    gap = Cp->SeqLenR - Cp->EndR;
	  else
	    {
	      if ( Cp->linksTo->isLinkHead )
		gap = (Cp->SeqLenR + Cp->linksTo->StartR) - Cp->EndR - 1;
	      else
		gap = Cp->linksTo->StartR - Cp->EndR - 1;
	    }

	  //-- Print the data
	  printf ("%ld\t%ld\t%ld\t%ld\t%.2f\t%.2f\t%c\t%s\n",
		  Cp->StartR, Cp->EndR, gap, len, pcov,
		  pidy, Cp->DirQ, Cp->IdQ);

	  //-- Walk the pointer down the list of contigs
	  if ( Cp->linksTo == Contigs.end( )  ||  Cp->linksTo->isLinkHead )
	    Cp = Contigs.end( );
	  else
	    Cp = Cp->linksTo;
	}
    }

  return;
}




void printTilingXML
     (vector<QueryContig> & Contigs, char * QryFileName, int argc, char ** argv)
{
  vector<AlignStats>::iterator Ap;
  vector<QueryContig>::iterator Cp;
  vector<QueryContig>::iterator beginCp;

  int i;
  long int gap;
  long int ct = 0;

  time_t tt = time(NULL);
  char * time_str;
  time_str = ctime(&tt);
  time_str[strlen(time_str)-1] = '\0';

  printf ("<?xml version = \"1.0\" ?>\n");

  printf ("<EVIDENCE ID = \"project_1\"\n");
  printf ("        DATE = \"%s\"\n", time_str);
  printf ("     PROJECT = \"%s\"\n", QryFileName);
  printf ("  PARAMETERS = \"");
  for ( i = 0; i < argc; i ++ ) printf ("%s%s", i == 0 ? "" : " ", argv[i]);
  printf ("\">\n");

  //-- Print the used contig information
  for ( beginCp = Contigs.begin( ); beginCp < Contigs.end( ); beginCp ++ )
    {
      if ( ! beginCp->isLinkHead )
	continue;

      //-- A new Reference sequence
      Cp = beginCp;

      //-- For all contigs mapping to this reference
      while ( Cp != Contigs.end( ) )
	{
	  printf
	    ("\t<CONTIG ID = \"contig_%s\" NAME = \"%s\" LEN = \"%ld\"/>\n",
	     Cp->IdQ, Cp->IdQ, Cp->SeqLenQ);

	  //-- Walk the pointer down the list of contigs
	  if ( Cp->linksTo == Contigs.end( )  ||  Cp->linksTo->isLinkHead )
	    Cp = Contigs.end( );
	  else
	    Cp = Cp->linksTo;
	}
    }
  printf("\n");


  //-- Print the linking information
  for ( beginCp = Contigs.begin( ); beginCp < Contigs.end( ); beginCp ++ )
    {
      if ( ! beginCp->isLinkHead )
	continue;

      //-- A new Reference sequence
      Cp = beginCp;

      //-- For all contigs mapping to this reference
      while ( Cp != Contigs.end( ) )
	{
	  if ( Cp->linksTo == Contigs.end( ) )
	    break;

	  //-- Calculate the gap size between this and the next contig
	  if ( Cp->linksTo->isLinkHead )
	    gap = (Cp->SeqLenR + Cp->linksTo->StartR) - Cp->EndR - 1;
	  else
	    gap = Cp->linksTo->StartR - Cp->EndR - 1;

	  printf
	    ("\t<LINK ID = \"link_%ld\" SIZE = \"%ld\" TYPE = \"MUMmer\">\n",
	     ++ ct, gap);
	  printf("\t\t<CONTIG ID = \"contig_%s\" ORI = \"%s\">\n",
		 Cp->IdQ, Cp->DirQ == FORWARD_CHAR ? "BE" : "EB");

	  for ( Ap = Cp->Aligns.begin( );
		Ap < Cp->Aligns.end( ); Ap ++ )
	    if ( ! Ap->isTiled )
	      continue;
	    else
	      {
		printf("\t\t");
		printAlignment(Cp,Ap,stdout);
	      }

	  printf("\t\t</CONTIG>\n");
	  printf("\t\t<CONTIG ID = \"contig_%s\" ORI = \"%s\">\n",
		 Cp->linksTo->IdQ,
		 Cp->linksTo->DirQ == FORWARD_CHAR ? "BE" : "EB");

	  for ( Ap = Cp->linksTo->Aligns.begin( );
		Ap < Cp->linksTo->Aligns.end( ); Ap ++ )
	    if ( ! Ap->isTiled )
	      continue;
	    else
	      {
		printf("\t\t");
		printAlignment(Cp->linksTo,Ap,stdout);
	      }

	  printf("\t\t</CONTIG>\n");
	  printf("\t</LINK>\n");

	  //-- Walk the pointer down the list of contigs
	  if ( Cp->linksTo == Contigs.end( )  ||  Cp->linksTo->isLinkHead )
	    Cp = Contigs.end( );
	  else
	    Cp = Cp->linksTo;
	}
    }

  printf ("</EVIDENCE>\n");

  return;
}




void printHelp
     (const char * s)

     //  Display the program's help information to stderr.

{
  fprintf (stderr,
           "\nUSAGE: %s  [options]  <deltafile>\n\n", s);
  fprintf (stderr,
    "-a            Describe the tiling path by printing the tab-delimited\n"
    "              alignment region coordinates to stdout\n"
    "-c            Assume the reference sequences are circular, and allow\n"
    "              tiled contigs to span the origin\n"
    "-h            Display help information\n"
    "-g int        Set maximum gap between clustered alignments [-1, INT_MAX]\n"
    "              A value of -1 will represent infinity\n"
    "              (nucmer default = %ld)\n"
    "              (promer default = %ld)\n"
    "-i float      Set minimum percent identity to tile [0.0, 100.0]\n"
    "              (nucmer default = %.1f)\n"
    "              (promer default = %.1f)\n"
    "-l int        Set minimum length contig to report [-1, INT_MAX]\n"
    "              A value of -1 will represent infinity\n"
    "              (common default = %ld)\n"
    "-p file       Output a pseudo molecule of the query contigs to 'file'\n"
    "-R            Deal with repetitive contigs by randomly placing them\n"
    "              in one of their copy locations (implies -V 0)\n"
    "-t file       Output a TIGR style contig list of each query sequence\n"
    "              that sufficiently matches the reference (non-circular)\n"
    "-u file       Output the tab-delimited alignment region coordinates\n"
    "              of the unusable contigs to 'file'\n"
    "-v float      Set minimum contig coverage to tile [0.0, 100.0]\n"
    "              (nucmer default = %.1f) sum of individual alignments\n"
    "              (promer default = %.1f) extent of syntenic region\n"
    "-V float      Set minimum contig coverage difference [0.0, 100.0]\n"
    "              i.e. the difference needed to determine one alignment\n"
    "              is 'better' than another alignment\n"
    "              (nucmer default = %.1f) sum of individual alignments\n"
    "              (promer default = %.1f) extent of syntenic region\n"
    "-x            Describe the tiling path by printing the XML contig\n"
    "              linking information to stdout\n\n",
	   DEFAULT_NUCMER_MAX_GAP_SIZE, DEFAULT_PROMER_MAX_GAP_SIZE,
	   DEFAULT_NUCMER_MIN_PIDY, DEFAULT_PROMER_MIN_PIDY,
	   DEFAULT_MIN_CONTIG_LENGTH,
	   DEFAULT_NUCMER_MIN_COVERAGE, DEFAULT_PROMER_MIN_COVERAGE,
	   DEFAULT_NUCMER_MIN_COVERAGE_DIFF, DEFAULT_PROMER_MIN_COVERAGE_DIFF);
  fprintf (stderr,
       "  Input is the .delta output of the nucmer program, run on very\n"
       "similar sequence data, or the .delta output of the promer program,\n"
       "run on divergent sequence data.\n"
       "  Output is to stdout, and consists of the predicted location of\n"
       "each aligning query contig as mapped to the reference sequences.\n"
       "These coordinates reference the extent of the entire query contig,\n"
       "even when only a certain percentage of the contig was actually\n"
       "aligned (unless the -a option is used). Columns are, start in ref,\n"
       "end in ref, distance to next contig, length of this contig, alignment\n"
       "coverage, identity, orientation, and ID respectively.\n\n");
  return;
}




void printUsage
     (const char * s)

     //  Display the program's usage information to stderr.

{
  fprintf (stderr,
           "\nUSAGE: %s  [options]  <deltafile>\n\n", s);
  fprintf (stderr, "Try '%s -h' for more information.\n", s);
  return;
}




inline long int revC
     (long int Coord, long int Len)

     //  Reverse complement the given coordinate for the given length.

{
  assert (Len - Coord + 1 > 0);
  return (Len - Coord + 1);
}




void tileContigs
     (vector<QueryContig> & Contigs)
{
  vector<QueryContig>::iterator Cp;
  vector<AlignStats>::iterator Ap, Aip, Atemp;

  long int start, end, hang;

  char * IdR = NULL;
  char * IdRhang = NULL;

  float cov, covhang;
  float max_cov, maxx_cov;
  float idy;
  float max_covhang, maxx_covhang;
  float idyhang, tmpidy;

  for ( Cp = Contigs.begin( ); Cp < Contigs.end( ); Cp ++ )
    {
      if ( Cp->TileLevel != UNUSED_TILE_LEVEL )
	continue;

      max_cov = maxx_cov = idy = -(FLT_MAX);
      max_covhang = maxx_covhang = idyhang = -(FLT_MAX);  
      Ap = Cp->Aligns.begin( );
      while ( Ap < Cp->Aligns.end( ) )
	{
	  //-- For a single reference
	  for ( Aip = Ap + 1; Aip < Cp->Aligns.end( ); Aip ++ )
	    if ( strcmp (Ap->IdR, Aip->IdR) != 0 )
	      break;

	  //-- Cluster the alignments
	  longestConsistentSubset (Ap, Aip);
	  
	  //-- Get the extents and overhang of the best subset
	  hang = 0;
	  end = start = -1;
	  for ( Atemp = Ap; Atemp < Aip; Atemp ++ )
	    if ( Atemp->isTiled )
	      {
		if ( Atemp->sQ - Atemp->sR > 0 )
		  hang += Atemp->sQ - Atemp->sR;
		start = Atemp->sR;
		break;
	      }
	  for ( Atemp = Aip - 1; Atemp >= Ap; Atemp -- )
	    if ( Atemp->isTiled )
	      {
		if ( Atemp->eR + (Cp->SeqLenQ - Atemp->eQ) > Atemp->SeqLenR )
		  hang += (Atemp->eR + (Cp->SeqLenQ - Atemp->eQ)) -
		    Atemp->SeqLenR;
		end = Atemp->eR;
		break;
	      }
	  assert ( end > 0 && start > 0 );

	  //-- Check the query coverage for mapping to this reference
	  if ( DATA_TYPE == NUCMER_DATA )
	    cov = getSubsetQryCoverage (Ap, Aip, Cp->SeqLenQ);
	  else // PROMER_DATA
	    cov = getSubsetQrySyntenyCoverage (Ap, Aip, Cp->SeqLenQ);
	  if ( cov >= max_cov )
	    {
	      tmpidy = getSubsetIdentity (Ap, Aip);
	      if ( cov > max_cov  ||  (cov == max_cov  &&  tmpidy > idy) )
		{
		  IdR = Ap->IdR;
		  maxx_cov = max_cov;
		  max_cov = cov;
		  idy = tmpidy;
		}
	    }
	  else if ( cov > maxx_cov )
	    maxx_cov = cov;

	  if ( !isCircularReference )
	    {
	      //-- Check the query coverage for mapping to this reference
	      if ( DATA_TYPE == NUCMER_DATA )
		covhang = getSubsetQryCoverage
		  (Ap, Aip, Cp->SeqLenQ - hang);
	      else // PROMER_DATA
		covhang = getSubsetQrySyntenyCoverage
		  (Ap, Aip, Cp->SeqLenQ - hang);
	      if ( covhang >= max_covhang )
		{
		  tmpidy = getSubsetIdentity (Ap, Aip);
		  if ( covhang > max_covhang  ||  (covhang == max_covhang  &&
						   tmpidy > idyhang) )
		    {
		      IdRhang = Ap->IdR;
		      maxx_covhang = max_covhang;
		      max_covhang = covhang;
		      idyhang = tmpidy;
		    }
		}
	      else if ( covhang > maxx_covhang )
		maxx_covhang = covhang;
	    }

	  for ( ; Ap < Aip; Ap ++ )
	    if ( !Ap->isTiled )
	      {
		//-- Don't worry about the single if it overlaps
		if (start <= end  &&  (Ap->sR <= end  &&  Ap->eR >= start))
		  continue;
		else if (start > end  &&  (Ap->eR >= start  ||  Ap->sR <= end))
		  continue;

		//-- Look for competing single alignments
		cov = getAlignmentQryCoverage (Ap, Cp->SeqLenQ);
		if ( cov > maxx_cov )
		  maxx_cov = cov;

		if ( !isCircularReference )
		  {
		    hang = 0;
		    if ( Ap->sQ - Ap->sR > 0 )
		      hang += Ap->sQ - Ap->sR;
		    if ( Ap->eR + (Cp->SeqLenQ - Ap->eQ) >
			 Ap->SeqLenR )
		      hang += (Ap->eR + (Cp->SeqLenQ - Ap->eQ)) -
			Ap->SeqLenR;

		    covhang = getAlignmentQryCoverage (Ap, Cp->SeqLenQ - hang);
		    if ( covhang > maxx_covhang )
		      maxx_covhang = covhang;
		  }
	      }
	  // Ap now equals Aip
	}

      //-- If clustered coverage is...
      if ( max_cov - maxx_cov >= MIN_COVERAGE_DIFF  &&
	   max_cov >= MIN_COVERAGE  &&  idy >= MIN_PIDY )
	{	  
	  for ( Aip = Cp->Aligns.begin( ); Aip < Cp->Aligns.end( ); Aip ++ )
	    if ( Aip->isTiled  &&  strcmp (Aip->IdR, IdR) != 0 )
	      Aip->isTiled = false;

	  //-- Tile the contig
	  Cp->TileLevel = USED_TILE_LEVEL;
	  placeContig ( Cp );
	}
      else if ( !isCircularReference  &&
		max_covhang - maxx_covhang >= MIN_COVERAGE_DIFF  &&
		max_covhang >= MIN_COVERAGE  &&  idyhang >= MIN_PIDY )
	{
	  for ( Aip = Cp->Aligns.begin( ); Aip < Cp->Aligns.end( ); Aip ++ )
	    if ( Aip->isTiled  &&  strcmp (Aip->IdR, IdRhang) != 0 )
	      Aip->isTiled = false;

	  //-- Tile the contig
	  Cp->TileLevel = USED_TILE_LEVEL;
	  placeContig ( Cp );
	}
      else
	Cp->TileLevel = UNUSABLE_TILE_LEVEL;
    }

  //-- Sort by reference mapping location
  sort ( Contigs.begin( ), Contigs.end( ), IdR_StartR_Sort( ) );

  return;
}
