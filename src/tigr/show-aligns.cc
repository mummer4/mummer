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

#include "delta.hh"
#include "tigrinc.hh"
#include "translate.hh"
#include "sw_alignscore.hh"
#include <vector>
#include <algorithm>
using namespace std;

//------------------------------------------------------------- Constants ----//

const char NUCMER_MISMATCH_CHAR = '^';
const char NUCMER_MATCH_CHAR = ' ';
const char PROMER_SIM_CHAR = '+';
const char PROMER_MISMATCH_CHAR = ' ';

//-- Note: if coord exceeds LINE_PREFIX_LEN - 1 digits,
//         increase these accordingly
#define LINE_PREFIX_LEN 11
#define PREFIX_FORMAT "%-10ld "

#define DEFAULT_SCREEN_WIDTH 60
int Screen_Width = DEFAULT_SCREEN_WIDTH;



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
int MATRIX_TYPE = BLOSUM62;

char InputFileName [MAX_LINE];
char RefFileName [MAX_LINE], QryFileName [MAX_LINE];



//------------------------------------------------- Function Declarations ----//
long int toFwd
     (long int coord, long int len, int frame);

void parseDelta
     (vector<AlignStats> & Aligns, char * IdR, char * IdQ);

void printAlignments
     (vector<AlignStats> Aligns, char * R, char * Q);

void printHelp
     (const char * s);

void printUsage
     (const char * s);

long int revC
     (long int coord, long int len);



//-------------------------------------------------- Function Definitions ----//
int main
     (int argc, char ** argv)
{
  long int i;

  FILE * RefFile = NULL;
  FILE * QryFile = NULL;

  vector<AlignStats> Aligns;

  char * R, * Q;

  long int InitSize = INIT_SIZE;
  char Id [MAX_LINE], IdR [MAX_LINE], IdQ [MAX_LINE];

  //-- Parse the command line arguments
  {
    int ch, errflg = 0;
    optarg = NULL;

    while ( !errflg  &&  ((ch = getopt
                           (argc, argv, "hqrw:x:")) != EOF) )
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
		      DEFAULT_SCREEN_WIDTH);
	      Screen_Width = DEFAULT_SCREEN_WIDTH;
	    }
	  break;

	case 'x' :
	  MATRIX_TYPE = atoi (optarg);
	  if ( MATRIX_TYPE < 1 || MATRIX_TYPE > 3 )
	    {
	      fprintf(stderr,
		      "WARNING: invalid matrix type %d, using default\n",
		      MATRIX_TYPE);
	      MATRIX_TYPE = BLOSUM62;
	    }
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
  }

  strcpy (InputFileName, argv[optind ++]);
  strcpy (IdR, argv[optind ++]);
  strcpy (IdQ, argv[optind ++]);

  //-- Read in the alignment data
  parseDelta (Aligns, IdR, IdQ);

  //-- Find, and read in the reference sequence
  RefFile = File_Open (RefFileName, "r");
  InitSize = INIT_SIZE;
  R = (char *) Safe_malloc ( sizeof(char) * InitSize );
  while ( Read_String (RefFile, R, InitSize, Id, FALSE) )
    if ( strcmp (Id, IdR) == 0 )
      break;
  fclose (RefFile);
  if ( strcmp (Id, IdR) != 0 )
    {
      fprintf(stderr,"ERROR: Could not find %s in the reference file\n", IdR);
      exit (EXIT_FAILURE);
    }


  //-- Find, and read in the query sequence
  QryFile = File_Open (QryFileName, "r");
  InitSize = INIT_SIZE;
  Q = (char *) Safe_malloc ( sizeof(char) * InitSize );
  while ( Read_String (QryFile, Q, InitSize, Id, FALSE) )
    if ( strcmp (Id, IdQ) == 0 )
      break;
  fclose (QryFile);
  if ( strcmp (Id, IdQ) != 0 )
    {
      fprintf(stderr,"ERROR: Could not find %s in the query file\n", IdQ);
      exit (EXIT_FAILURE);
    }

  //-- Sort the alignment regions if user passed -r or -q option
  if ( isSortByReference )
    sort (Aligns.begin( ), Aligns.end( ), sR_Sort( ));
  else if ( isSortByQuery )
    sort (Aligns.begin( ), Aligns.end( ), sQ_Sort( ));


  //-- Output the alignments to stdout
  printf("%s %s\n\n", RefFileName, QryFileName);
  for ( i = 0; i < Screen_Width; i ++ ) printf("=");
  printf("\n-- Alignments between %s and %s\n\n", IdR, IdQ);
  printAlignments (Aligns, R, Q);
  printf("\n");
  for ( i = 0; i < Screen_Width; i ++ ) printf("=");
  printf("\n");

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
     (vector<AlignStats> & Aligns, char * IdR, char * IdQ)

     // Read in the alignments from the desired region

{
  AlignStats aStats;                     //  single alignment region
  bool found = false;

  DeltaReader_t dr;
  dr.open (InputFileName);
  DATA_TYPE = dr.getDataType( ) == NUCMER_STRING ?
    NUCMER_DATA : PROMER_DATA;
  strcpy (RefFileName, dr.getReferencePath( ).c_str( ));
  strcpy (QryFileName, dr.getQueryPath( ).c_str( ));

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
	      IdR, IdQ);
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
     (vector<AlignStats> Aligns, char * R, char * Q)

     // Print the alignments to the screen

{
  vector<AlignStats>::iterator Ap;
  vector<long int>::iterator Dp;

  char * A[7] = {NULL, NULL, NULL, NULL, NULL, NULL, NULL};
  char * B[7] = {NULL, NULL, NULL, NULL, NULL, NULL, NULL};
  int Ai, Bi, i;

  char Buff1 [Screen_Width + 1],
       Buff2 [Screen_Width + 1],
       Buff3 [Screen_Width + 1];

  int Sign;
  long int Delta;
  long int Total, Errors, Remain;
  long int Pos;

  long int sR, eR, sQ, eQ;
  long int Apos, Bpos;
  long int SeqLenR, SeqLenQ;
  int frameR, frameQ;

  for ( i = 0; i < LINE_PREFIX_LEN; i ++ )
    Buff3[i] = ' ';

  SeqLenR = strlen (R + 1);
  SeqLenQ = strlen (Q + 1);

  if ( DATA_TYPE == NUCMER_DATA )
    {
      A[1] = R;
      A[4] = (char *) Safe_malloc ( sizeof(char) * (SeqLenR + 2) );
      strcpy ( A[4] + 1, A[1] + 1 );
      A[4][0] = '\0';
      Reverse_Complement ( A[4], 1, SeqLenR );

      B[1] = Q;
      B[4] = (char *) Safe_malloc ( sizeof(char) * (SeqLenQ + 2) );
      strcpy ( B[4] + 1, B[1] + 1 );
      B[4][0] = '\0';
      Reverse_Complement ( B[4], 1, SeqLenQ );
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
	  A[Ai] = (char *) Safe_malloc ( sizeof(char) * ( SeqLenR / 3 + 2 ) );
	  A[Ai][0] = '\0';
	  Translate_DNA ( R, A[Ai], Ai );
	}
      if ( B[Bi] == NULL )
	{
	  assert ( DATA_TYPE == PROMER_DATA );
	  B[Bi] = (char *) Safe_malloc ( sizeof(char) * ( SeqLenQ / 3 + 2 ) );
	  B[Bi][0] = '\0';
	  Translate_DNA ( Q, B[Bi], Bi );
	}


      //-- Generate the alignment
      printf("-- BEGIN alignment [ %s%d %ld - %ld | %s%d %ld - %ld ]\n\n\n",
	     frameR > 0 ? "+" : "-", abs(frameR), Ap->sR, Ap->eR,
	     frameQ > 0 ? "+" : "-", abs(frameQ), Ap->sQ, Ap->eQ);

      Apos = sR;
      Bpos = sQ;

      Errors = 0;
      Total = 0;
      Remain = eR - sR + 1;

      sprintf(Buff1, PREFIX_FORMAT, toFwd (Apos, SeqLenR, frameR));
      sprintf(Buff2, PREFIX_FORMAT, toFwd (Bpos, SeqLenQ, frameQ));
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
	      if ( Pos >= Screen_Width )
		{
		  Buff1[Pos] = Buff2[Pos] = Buff3[Pos] = '\0';
		  if ( DATA_TYPE == NUCMER_DATA )
		    printf("%s\n%s\n%s\n\n", Buff1, Buff2, Buff3);
		  else
		    printf("%s\n%s\n%s\n\n", Buff1, Buff3, Buff2);
		  sprintf(Buff1, PREFIX_FORMAT, toFwd (Apos, SeqLenR, frameR));
		  sprintf(Buff2, PREFIX_FORMAT, toFwd (Bpos, SeqLenQ, frameQ));
		  Pos = LINE_PREFIX_LEN;
		}

	      if ( DATA_TYPE == NUCMER_DATA )
		Buff3[Pos] = A[Ai][Apos] == B[Bi][Bpos] ?
		  NUCMER_MATCH_CHAR : NUCMER_MISMATCH_CHAR;
	      else if ( A[Ai][Apos] == B[Bi][Bpos] )
		Buff3[Pos] = A[Ai][Apos];
	      else
		Buff3[Pos] = MATCH_SCORE
		  [MATRIX_TYPE]
		  [toupper(A[Ai][Apos]) - 'A']
		  [toupper(B[Bi][Bpos]) - 'A'] > 0 ?
		  PROMER_SIM_CHAR : PROMER_MISMATCH_CHAR;
	      Buff1[Pos] = A[Ai][Apos ++];
	      Buff2[Pos ++] = B[Bi][Bpos ++];
	    }


	  //-- For the indel
	  Remain -= i - 1;

	  if ( Pos >= Screen_Width )
	    {
	      Buff1[Pos] = Buff2[Pos] = Buff3[Pos] = '\0';
	      if ( DATA_TYPE == NUCMER_DATA )
		printf("%s\n%s\n%s\n\n", Buff1, Buff2, Buff3);
	      else
		printf("%s\n%s\n%s\n\n", Buff1, Buff3, Buff2);
	      sprintf(Buff1, PREFIX_FORMAT, toFwd (Apos, SeqLenR, frameR));
	      sprintf(Buff2, PREFIX_FORMAT, toFwd (Bpos, SeqLenQ, frameQ));
	      Pos = LINE_PREFIX_LEN;
	    }

	  if ( Sign == 1 )
	    {
	      if ( DATA_TYPE == NUCMER_DATA )
		Buff3[Pos] = NUCMER_MISMATCH_CHAR;
	      else
		Buff3[Pos] = PROMER_MISMATCH_CHAR;
	      Buff1[Pos] = A[Ai][Apos ++];
	      Buff2[Pos ++] = '.';
	      Remain --;
	    }
	  else
	    {
	      if ( DATA_TYPE == NUCMER_DATA )
		Buff3[Pos] = NUCMER_MISMATCH_CHAR;
	      else
		Buff3[Pos] = PROMER_MISMATCH_CHAR;
	      Buff1[Pos] = '.';
	      Buff2[Pos ++] = B[Bi][Bpos ++];
	      Total ++;
	    }
	}


      //-- For all the bases remaining after the last indel
      for ( i = 0; i < Remain; i ++ )
	{
	  if ( Pos >= Screen_Width )
	    {
	      Buff1[Pos] = Buff2[Pos] = Buff3[Pos] = '\0';
	      if ( DATA_TYPE == NUCMER_DATA )
		printf("%s\n%s\n%s\n\n", Buff1, Buff2, Buff3);
	      else
		printf("%s\n%s\n%s\n\n", Buff1, Buff3, Buff2);
	      sprintf(Buff1, PREFIX_FORMAT, toFwd (Apos, SeqLenR, frameR));
	      sprintf(Buff2, PREFIX_FORMAT, toFwd (Bpos, SeqLenQ, frameQ));
	      Pos = LINE_PREFIX_LEN;
	    }
	  
	  if ( DATA_TYPE == NUCMER_DATA )
	    Buff3[Pos] = A[Ai][Apos] == B[Bi][Bpos] ?
	      NUCMER_MATCH_CHAR : NUCMER_MISMATCH_CHAR;
	  else if ( A[Ai][Apos] == B[Bi][Bpos] )
	    Buff3[Pos] = A[Ai][Apos];
	  else
	    Buff3[Pos] = MATCH_SCORE
	      [MATRIX_TYPE]
	      [toupper(A[Ai][Apos]) - 'A']
	      [toupper(B[Bi][Bpos]) - 'A'] > 0 ?
	      PROMER_SIM_CHAR : PROMER_MISMATCH_CHAR;
	  Buff1[Pos] = A[Ai][Apos ++];
	  Buff2[Pos ++] = B[Bi][Bpos ++];
	}


      //-- For the remaining buffered output
      if ( Pos > LINE_PREFIX_LEN )
	{
	  Buff1[Pos] = Buff2[Pos] = Buff3[Pos] = '\0';
	  if ( DATA_TYPE == NUCMER_DATA )
	    printf("%s\n%s\n%s\n\n", Buff1, Buff2, Buff3);
	  else
	    printf("%s\n%s\n%s\n\n", Buff1, Buff3, Buff2);
	  sprintf(Buff1, PREFIX_FORMAT, toFwd (Apos, SeqLenR, frameR));
	  sprintf(Buff2, PREFIX_FORMAT, toFwd (Bpos, SeqLenQ, frameQ));
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
	free ( A[i] );
      if ( (DATA_TYPE != NUCMER_DATA || i != 1)  &&  B[i] != NULL )
	free ( B[i] );
    }

  return;
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
       "-w int        Set the screen width - default is 60\n"
       "-x int        Set the matrix type - default is 2 (BLOSUM 62),\n"
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
