//------------------------------------------------------------------------------
//   Programmer: Adam M Phillippy, The Institute for Genomic Research
//         File: show-coords.cc
//         Date: 08 / 13 / 2002
//
//        Usage: show-coords [options] <deltafile>
//               Try 'show-coords -h' for more information
//
//  Description: For use in conjunction with the MUMmer package.
//              "show-coords" displays human readable information from the
//             .delta output of the "nucmer" and "promer" programs. Outputs
//            various useful information to stdout. Works for both nucleotide
//           and amino-acid alignments.
//
//------------------------------------------------------------------------------

#include <mummer/delta.hh>
#include <mummer/tigrinc.hh>
#include <mummer/redirect_to_pager.hpp>
#include <vector>
#include <algorithm>
using namespace std;

//------------------------------------------------------------- Constants ----//

//-- for the -k option, used in comparing two overlapping alignments
const float MIN_OVERLAP_P = 0.5; // min percent overlap to spawn knockout
const float MIN_LENGTH_P = 0.75; // winner must be atleast this big as loser

const float MAX_OHANG_P = 0.05; // max overlap hang as a percentage of overlap
const float MAX_PIDYDIFF = 0.01;


//------------------------------------------------------ Type Definitions ----//
struct AlignStats
     //-- Alignment statistics data structure
{
  bool isLAS;               // involved in a longest ascending subset
  float Idy;                              // percent identity (0.0 - 100.0)
  float Sim;                              // percent similarity (0.0 - 100.0)
  float Stp;                              // precent stop codon (0.0 - 100.0)
  int FrameA, FrameB;                     // reading frame
  long int sA, eA, sB, eB;                // start, end in A, start, end in B
  long int SeqLenA, SeqLenB;              // length of seq A, seq B
  char * IdA, * IdB;                      // Id of seq A, Id of seq B
  char annot [12];                        // annotation string
};


struct LASstats
{
  bool ori;
  long int score, wscore;
  long int from, wfrom;
  bool wpoint;
  vector<AlignStats>::iterator asi;
};


struct LAS_Sort
{
  bool operator() (const LASstats & pA, const LASstats & pB)
  {
    if ( pA . asi -> sA < pB . asi -> sA )
      return true;
    else
      return false;
  }
};


struct sA_Sort
     //-- For sorting AlignStats by sA
{
  bool operator() (const AlignStats & pA, const AlignStats & pB)
  {
    //-- sort sA
    if ( pA.sA < pB.sA )
      return true;
    else
      return false;
  }
};


struct IdA_sA_IdB_sB_Sort
     //-- For sorting AlignStats by IdA, sA, IdB, sB in that order (-r option)
{
  bool operator() (const AlignStats & pA, const AlignStats & pB)
  {
    int i = strcmp (pA.IdA, pB.IdA);
    int j = strcmp (pA.IdB, pB.IdB);

    //-- sort IdA
    if ( i < 0 )
      return true;
    else if ( i > 0 )
      return false;
    //-- sort sA
    else
      {
	if ( pA.sA < pB.sA )
	  return true;
	else if ( pA.sA > pB.sA )
	  return false;
	//-- sort IdB
	else
	  {
	    if ( j < 0 )
	      return true;
	    else if ( j > 0 )
	      return false;
	    //-- sort sB
	    else
	      {
		if ( pA.sB < pB.sB )
		  return true;
		else
		  return false;
	      }
	  }
      }
  }
};


struct IdB_sB_IdA_sA_Sort
     //-- For sorting AlignStats by IdB, sB, IdA, sA in that order (-q option)
{
  bool operator() (const AlignStats & pA, const AlignStats & pB)
  {
    int i = strcmp (pA.IdB, pB.IdB);
    int j = strcmp (pA.IdA, pB.IdA);

    //-- sort IdB
    if ( i < 0 )
      return true;
    else if ( i > 0 )
      return false;
    //-- sort sB
    else
      {
	if ( pA.sB < pB.sB )
	  return true;
	else if ( pA.sB > pB.sB )
	  return false;
	//-- sort IdA
	else
	  {
	    if ( j < 0 )
	      return true;
	    else if ( j > 0 )
	      return false;
	    //-- sort sA
	    else
	      {
		if ( pA.sA < pB.sA )
		  return true;
		else
		  return false;
	      }
	  }
      }
  }
};




//--------------------------------------------------- Global Option Flags ----//
bool isKnockout = false;                // -k option
bool isAnnotateOverlaps = false;        // -o option
bool isBrief = false;                   // -b option
bool isBtab = false;                    // -B option
bool isShowCoverage = false;            // -c option
bool isShowDir = false;                 // -d option
bool isShowSeqLens = false;             // -l option
bool isShowWarnings = false;            // -w option
bool isSortByQuery = false;             // -q option
bool isSortByReference = false;         // -r option
bool isLAS = false;                     // -g option
bool isWLAS = false;                    // -G option
bool isPrintTabular = false;            // -T option
bool isPrintHeader = true;              // -H option
bool isAnnotation = false;              // true if either -w or -o
float idyCutoff = 0;                    // -I option
long int lenCutoff = 0;                 // -L option
int  whichDataType = NUCMER_DATA;       // set by .delta header
char InputFileName [MAX_LINE];          //  I/O filenames
char RefFileName [MAX_LINE], QryFileName [MAX_LINE];



//------------------------------------------------- Function Declarations ----//
void generateWarnings
     (vector<AlignStats> & Stats);

inline long int  Max
     (long int A, long int B) { if  (A <  B) return  B; else return  A; }

inline long int  Min
     (long int A, long int B) { if  (A >= B) return  B; else return  A; }

void annotateOverlaps
     (vector<AlignStats> & Stats);

void flagLAS
     (vector<AlignStats> & Stats);

void parseDelta
     (vector<AlignStats> & Stats);

void printBtab
     (vector<AlignStats> Stats);

void printHuman
     (vector<AlignStats> Stats);

void printTabular
     (vector<AlignStats> Stats);

inline long int revC
     (long int Coord, long int Len);

void simplifyAlignments
     (vector<AlignStats> & Stats);

void printHelp
     (const char * s);

void printUsage
     (const char * s);

inline void swap
     (long int & a, long int & b)
{
  static long int t;
  t = a; a = b; b = t;
}





//-------------------------------------------------- Function Definitions ----//
int main
     (int argc, char ** argv)
{  
  vector<AlignStats> Stats;            //  all the alignment regions

  //-- Parse the command line arguments
  {
    int ch, errflg = 0;
    optarg = NULL;

    while ( !errflg &&
	    ((ch = getopt (argc, argv, "bkBdhTHqrgGclowI:L:")) != EOF) )
      switch (ch)
	{
	case 'b' :
	  isBrief = true;
	  break;

	case 'B' :
	  isBtab = true;
	  break;

	case 'c' :
	  isShowCoverage = true;
	  break;

	case 'd' :
	  isShowDir = true;
	  break;

	case 'h' :
	  printHelp (argv[0]);
	  exit (EXIT_SUCCESS);
	  break;

	case 'H' :
	  isPrintHeader = false;
	  break;

	case 'I' :
	  idyCutoff = atof (optarg);
	  break;

	case 'k' :
	  isKnockout = true;
	  break;

	case 'l' :
	  isShowSeqLens = true;
	  break;

	case 'L' :
	  lenCutoff = atol (optarg);
	  break;

	case 'o' :
	  isAnnotateOverlaps = true;
	  break;

	case 'q' :
	  isSortByQuery = true;
	  break;

	case 'r' :
	  isSortByReference = true;
	  break;

	case 'g' :
	  isLAS = true;
	  break;

	case 'G' :
	  isLAS = true;
	  isWLAS = true;
	  break;

	case 'T' :
	  isPrintTabular = true;
	  break;

	case 'w' :
	  isShowWarnings = true;
	  break;

	default :
	  errflg ++;
	}

    if ( errflg > 0  ||  argc - optind != 1 )
      {
	printUsage (argv[0]);
	exit (EXIT_FAILURE);
      }
    
    if ( isShowWarnings && isAnnotateOverlaps )
      {
	fprintf (stderr,
		 "WARNING: -o and -w are mutually exclusive, -w ignored\n");
	isShowWarnings = false;
      }

    if ( isWLAS || isShowWarnings )
      {
	fprintf (stderr,
	"WARNING: hidden option used, I hope you know what you are doing\n");
      }

    if ( isLAS && (!isSortByReference && !isSortByQuery) )
      {
	fprintf (stderr,
		 "WARNING: -g option requires -r or -q, using -r\n");
	isSortByReference = true;
      }

    if ( isBtab && (!isPrintHeader || isShowCoverage || isShowSeqLens ||
		    isPrintTabular || isShowWarnings || isBrief ||
		    isAnnotateOverlaps) )
      {
	fprintf (stderr,
		 "WARNING: the options -bdHcloTw are all ignored for btab\n");
	isPrintHeader = isShowCoverage = isShowSeqLens =
	  isPrintTabular = isShowWarnings = isBrief =
	  isAnnotateOverlaps = false;
      }

    if ( isBrief && (isKnockout || isLAS) )
      {
	fprintf (stderr,
   "WARNING: -k and -g have no effect when used with -b, -k and -g ignored\n");
	isKnockout = false;
	isLAS = false;
      }

    if ( isSortByReference && isSortByQuery )
      fprintf (stderr,
	       "WARNING: both -r and -q were passed, -q ignored\n");

    if ( idyCutoff < 0.0 || idyCutoff > 100.0 )
      {
	fprintf(stderr,
		"\nERROR: Percentages must be within the range [0.0. 100.0]\n");
	exit (EXIT_FAILURE);
      }

    if ( lenCutoff < 0 )
      {
	fprintf(stderr,
		"\nERROR: Length values must be positive\n");
	exit (EXIT_FAILURE);
      }

    if ( isShowWarnings || isAnnotateOverlaps )
      isAnnotation = true;
  }

  srand (time (NULL));

  //-- Open and parse the delta file
  strcpy (InputFileName, argv[optind ++]);
  parseDelta (Stats);

  //-- Can only pick best frame for promer data
  if ( isKnockout && whichDataType != PROMER_DATA )
    {
      fprintf(stderr,"WARNING: -k only applies to PROmer data, -k ignored\n");
      isKnockout = false;
    }

  if ( whichDataType == PROMER_DATA )
    isShowDir = true;

  //-- NOTE: simplifyAlignments assumes all alignments from two
  //   sequences are grouped together (as ouput be postpro/nuc)
  if ( isBrief || isKnockout )
    simplifyAlignments (Stats);

  //-- NOTE: flagLAS assumes all alignments from two
  //   sequences are grouped together (as ouput be postpro/nuc)
  if ( isLAS )
    flagLAS (Stats);

  //-- Sort the alignment regions if user passed -r or -q option
  if ( isSortByReference )
    sort (Stats.begin( ), Stats.end( ), IdA_sA_IdB_sB_Sort( ));
  else if ( isSortByQuery )
    sort (Stats.begin( ), Stats.end( ), IdB_sB_IdA_sA_Sort( ));


  //-- Generate overlap warnings if user passed -w option
  if ( isShowWarnings )
    generateWarnings (Stats);

  //-- Generate overlap annotations if user passed -o option
  if ( isAnnotateOverlaps )
    annotateOverlaps (Stats);

  //-- Output data to stdout, tabular if -T option was used
  stdio_launch_pager redirect_to_pager;
  if ( isBtab )
    printBtab (Stats);
  else if ( isPrintTabular )
    printTabular (Stats);
  else
    printHuman (Stats);
  fclose(stdout);

  return EXIT_SUCCESS;
}




void flagLAS
     (vector<AlignStats> & Stats)
{
  LASstats alas;
  vector<LASstats> las;

  vector<long int> ties;
  vector<AlignStats>::iterator Sp;
  vector<AlignStats>::iterator Sep;
  long int i, j, n, best, next;
  long int lenA, lenB, len, olapA, olapB, olap;
  bool wrap;

  Sp = Stats.begin( );
  while ( Sp < Stats.end( ) )
    {
      las . clear( );
      for ( Sep = Sp; Sep < Stats.end( ); Sep ++ )
	{
	  if ( strcmp (Sp->IdA, Sep->IdA) != 0  ||
	       strcmp (Sp->IdB, Sep->IdB) != 0 )
	    break;

	  Sep -> isLAS = false;

	  if ( Sep -> sA < 0 )
	    continue;

	  if ( Sep -> FrameA < 0 )
	    swap (Sep -> sA, Sep -> eA);
	  if ( Sep -> FrameB < 0 )
	    swap (Sep -> sB, Sep -> eB);

	  if ( (Sep -> FrameA < 0 && Sep -> FrameB < 0) ||
	       (Sep -> FrameA > 0 && Sep -> FrameB > 0) )
	    alas . ori = 0;
	  else
	    {
	      alas . ori = 1;
	      swap (Sep -> sB, Sep -> eB);
	      Sep -> sB = revC (Sep -> sB, Sep -> SeqLenB);
	      Sep -> eB = revC (Sep -> eB, Sep -> SeqLenB);
	    }

	  lenA = Sep ->  eA - Sep ->  sA + 1;
	  lenB = Sep ->  eB - Sep ->  sB + 1;

	  alas . score = alas . wscore = lenA < lenB ? lenA : lenB;
	  alas . from = alas . wfrom = -1;
	  alas . wpoint = false;
	  alas . asi = Sep;

	  las . push_back (alas);
	}

      sort (las . begin( ), las . end( ), LAS_Sort( ));

      n = las . size( );
      for ( i = 0; i < n; i ++ )
	for ( j = 0; j < i; j ++ )
	  {
	    if ( las [i] . ori != las [j] . ori )
	      continue;

	    lenA = las [i] . asi ->  eA - las [i] . asi ->  sA + 1;
	    lenB = las [i] . asi ->  eB - las [i] . asi ->  sB + 1;
	    len = lenA < lenB ? lenA : lenB;

	    olapA = las [j] . asi -> eA - las [i] . asi -> sA + 1;
	    olap = olapA > 0 ? olapA : 0;
	    
	    if ( las [j] . score + len - olap > las [i] . wscore )
	      {
		las [i] . wfrom = j;
		las [i] . wscore = las [j] . score + len - olap;
		las [i] . wpoint = true;
	      }
	    
	    olapB = las [j] . asi -> eB - las [i] . asi -> sB + 1;
	    olap = olapB > olap ? olapB : olap;
	    
	    if ( las [j] . score + len - olap > las [i] . score )
	      {
		las [i] . from = j;
		las [i] . score = las [j] . score + len - olap;
	      }
	    if ( las [j] . wscore + len - olap > las [i] . wscore )
	      {
		las [i] . wfrom = j;
		las [i] . wscore = las [j] . wscore + len - olap;
		las [i] . wpoint = false;
	      }
	  }


      best = 0;
      if ( isWLAS )
	{
	  for ( i = 0; i < n; i ++ )
	    {
	      if ( las [i] . wscore > las [best] . wscore  ||
		   (las [i] . wscore == las [best] . wscore &&
		    las [i] . asi -> Idy > las [best] . asi -> Idy) )
		best = i;
	    }
	}
      else
	{
	  for ( i = 0; i < n; i ++ )
	    {	      
	      if ( las [i] . score > las [best] . score  ||
		   (las [i] . score == las [best] . score &&
		    las [i] . asi -> Idy > las [best] . asi -> Idy) )
		best = i;
	    }
	}

      ties . clear( );
      for ( i = 0; i < n; i ++ )
	{
	  if ( las [i] . ori == 1 )
	    {
	      swap (las [i] . asi -> sB, las [i] . asi -> eB);
	      las[i].asi->sB = revC (las[i].asi->sB,las[i].asi->SeqLenB);
	      las[i].asi->eB = revC (las[i].asi->eB,las[i].asi->SeqLenB);
	    }
	  if ( las [i] . asi -> FrameA < 0 )
	    swap (las [i] . asi -> sA, las [i] . asi -> eA);
	  if ( las [i] . asi -> FrameB < 0 )
	    swap (las [i] . asi -> sB, las [i] . asi -> eB);

	  if ( ((isWLAS  &&  las [i] . wscore == las [best] . wscore)  ||
		(!isWLAS  &&  las [i] . score == las [best] . score))
	       &&
	       las [i] . asi -> Idy >= las [best] . asi -> Idy - MAX_PIDYDIFF)
	    ties . push_back (i);
	}
      if ( ties . size( ) > 1 )
	best = ties
	  [(long int)((double)ties . size( ) * rand( ) / (RAND_MAX + 1.0))];

      if ( isWLAS )
	{
	  wrap = false;
	  for ( i = best; i >= 0; i = next )
	    {
	      las [i] . asi -> isLAS = true;
	      if ( wrap )
		{
		  next = las [i] . from;
		  las [i] . wpoint = false;
		}
	      else
		next = las [i] . wfrom;
	      
	      if ( las [i] . wpoint )
		wrap = true;
	    }
	}
      else
	{
	  for ( i = best; i >= 0; i = las [i] . from )
	    las [i] . asi -> isLAS = true;
	}

      Sp = Sep;
    }


  for ( Sp = Stats . begin( ); Sp < Stats . end( ); Sp ++ )
    if ( Sp -> sA > 0  &&  !Sp -> isLAS )
      Sp -> sA *= -1;
}




void annotateOverlaps
     (vector<AlignStats> & Stats)
{
  vector<AlignStats>::iterator Sip;
  long int loR, hiR, loQ, hiQ;
  long int max_ohang;

  for ( Sip = Stats . begin( ); Sip < Stats . end( ); Sip ++ )
    {
      Sip -> annot [0] = '\0';

      if ( Sip -> sA < 0 )
	continue;

      loR = Sip -> sA;
      hiR = Sip -> eA;
      loQ = Sip -> sB;
      hiQ = Sip -> eB;

      if ( Sip -> FrameA < 0 )
	swap (loR, hiR);
      if ( Sip -> FrameB < 0 )
	swap (loQ, hiQ);
      
      if ( !(Sip -> FrameA < 0 && Sip -> FrameB < 0) &&
	   !(Sip -> FrameA > 0 && Sip -> FrameB > 0) )
	{
	  swap (loQ, hiQ);
	  loQ = revC (loQ, Sip -> SeqLenB);
	  hiQ = revC (hiQ, Sip -> SeqLenB);
	}

      max_ohang =
	(long int)(MAX_OHANG_P * (float)(hiR - loR > hiQ - loQ ?
					 hiR - loR + 1 : hiQ - loQ + 1));

      if ( (loR <= max_ohang  ||
	    loQ <= max_ohang)
	   &&
	   (hiR + max_ohang > Sip -> SeqLenA  ||
	    hiQ + max_ohang > Sip -> SeqLenB) )
	{
	  if      ( loR <= max_ohang  &&
		    loQ <= max_ohang  &&
		    hiR + max_ohang > Sip -> SeqLenA  &&
		    hiQ + max_ohang > Sip -> SeqLenB )
	    strcpy (Sip -> annot, "[IDENTITY]");
	  else if ( loR <= max_ohang  &&
		    hiR + max_ohang > Sip -> SeqLenA )
	    strcpy (Sip -> annot, "[CONTAINED]");
	  else if ( loQ <= max_ohang  &&
		    hiQ + max_ohang > Sip -> SeqLenB )
	    strcpy (Sip -> annot, "[CONTAINS]");
	  else if ( loQ <= max_ohang )
	    strcpy (Sip -> annot, "[END]");
	  else if ( loR <= max_ohang )
	    strcpy (Sip -> annot, "[BEGIN]");
	  else
	    assert (false);
	}
    }
}



void generateWarnings
     (vector<AlignStats> & Stats)

     //  Generate warnings for each alignment region if it overlaps or in any
     //  other way conflicts with the alignment preceding it. Set its
     //  warning field appropriately.

{  
  vector<AlignStats>::iterator Sip;
  vector<AlignStats>::iterator Sprev;
  long int s1, e1, s2, e2, ps1, pe1, ps2, pe2;
  long int tmp;

  if ( !Stats.empty( ) )
    Stats.begin()->annot [0] = '\0';
  for ( Sip = Stats.begin( ) + 1; Sip < Stats.end( ); Sip ++ )
    {
      Sip->annot [0] = '\0';

      if ( Sip->sA < 0 )
	continue;

      for ( Sprev = Sip - 1;
	    Sprev >= Stats.begin( )  &&  Sprev->sA < 0;
	    Sprev -- )
	{ };

      if ( Sprev < Stats.begin( ) ||
	   //	   Sip->FrameA != Sprev->FrameA  ||
	   //	   Sip->FrameB != Sprev->FrameB  ||
	   strcmp ( Sip->IdA , Sprev->IdA ) != 0 ||
	   strcmp ( Sip->IdB , Sprev->IdB ) != 0 )
	continue;

      ps1 = Sprev->sA;
      pe1 = Sprev->eA;
      ps2 = Sprev->sB;
      pe2 = Sprev->eB;
      s1  = Sip->sA;
      e1  = Sip->eA;
      s2  = Sip->sB;
      e2  = Sip->eB;

      if ( Sprev->FrameA < 0 )
	{
	  tmp = ps1;
	  ps1 = pe1;
	  pe1 = tmp;
	}
      if ( Sprev->FrameB < 0 )
	{
	  tmp = ps2;
	  ps2 = pe2;
	  pe2 = tmp;
	}
      if ( Sip->FrameA < 0 )
	{
	  tmp = s1;
	  s1 = e1;
	  e1 = tmp;
	}
      if ( Sip->FrameB < 0 )
	{
	  tmp = s2;
	  s2 = e2;
	  e2 = tmp;
	}

      //-- Duplicate region
      if ( s1 == ps1 && s2 == ps2 && e1 == pe1 && e2 == pe2 )
	strcpy ( Sip->annot , "[DUPLICATE]" );
      //-- Shadowing region
      else if ( s1 < ps1 && s2 < ps2 && e1 > pe1 && e2 > pe2 )
	strcpy ( Sip->annot , "[CONTAINS]" );
      //-- Shadowed region
      else if ( ps1 < s1 && ps2 < s2 && pe1 > e1 && pe2 > e2 )
	strcpy ( Sip->annot , "[SHADOWED]" );
      //-- Overlapping region in reference sequence
      else if ( (s1 >= ps1 && s1 <= pe1) || (e1 >= ps1 && e1 <= pe1) )
	strcpy ( Sip->annot , "[OVERLAPS]" );
      //-- Overlapping region in query sequence
      else if ( (s2 >= ps2 && s2 <= pe2) || (e2 >= ps2 && e2 <= pe2) )
	strcpy ( Sip->annot , "[OVERLAPS]" );
      //-- else No Conflict
    }
}




void printBtab
     (vector<AlignStats> Stats)
{
  time_t currtime;
  const char * type;
  char date[MAX_LINE];
  long int len;
  vector<AlignStats>::iterator Sip;

  currtime = time(NULL);
  strftime (date, MAX_LINE, "%b %d %Y", localtime(&currtime));
  if ( whichDataType == NUCMER_DATA )
    type = "NUCMER";
  else if ( whichDataType == PROMER_DATA )
    type = "PROMER";
  else
    type = "NULL";
  
  //-- Print the stats
  for ( Sip = Stats.begin( ); Sip < Stats.end( ); Sip ++ )
    {
      if ( Sip->sA < 0 )
	continue;

      len = labs(Sip->eB - Sip->sB) + 1;

      //-- Output the stats for this alignment in btab format
      printf("%s\t%s\t%ld\t%s\t%s\t%s\t",
	     Sip->IdB, date, Sip->SeqLenB, type, RefFileName, Sip->IdA);

      printf("%ld\t%ld\t%ld\t%ld\t%f\t%f\t%ld\t0\t0\tNULL\t",
	     Sip->sB, Sip->eB, Sip->sA, Sip->eA,
	     Sip->Idy, Sip->Sim, len);

      printf("%d\t%s\t%ld\t0\t0\n",
	     whichDataType == NUCMER_DATA ? 0 : Sip->FrameA,
	     Sip->FrameB < 0 ? "Minus" : "Plus", Sip->SeqLenA);
    }

  return;
}




void printHuman
     (vector<AlignStats> Stats)

     //  Print stats in a human readable format

{
  long int len1, len2;
  float covA, covB;
  vector<AlignStats>::iterator Sip;

  //-- Print the output header
  if ( isPrintHeader )
    {
      printf ("%s %s\n%s\n\n", RefFileName, QryFileName,
	      whichDataType == NUCMER_DATA ? "NUCMER" : "PROMER");
      printf("%8s %8s  | ", "[S1]", "[E1]");
      printf("%8s %8s  | ", "[S2]", "[E2]");
      printf("%8s %8s  | ", "[LEN 1]", "[LEN 2]");
      if ( !isBrief )
	{
	  printf("%8s ", "[% IDY]");
	  if ( whichDataType == PROMER_DATA )
	    printf("%8s %8s ", "[% SIM]", "[% STP]");
	  printf(" | ");
	}
      if ( isShowSeqLens )
	printf("%8s %8s  | ", "[LEN R]", "[LEN Q]");
      if ( isShowCoverage )
	printf("%8s %8s  | ", "[COV R]", "[COV Q]");
      if ( !isBrief  &&  isShowDir )
	printf("%5s  ", "[FRM]");
      printf("%s", "[TAGS]");
      printf("\n");
      if ( isShowSeqLens )
	printf("=====================");
      if ( isShowCoverage )
	printf("=====================");
      if ( !isBrief )
	{
	  printf("============");
	  if ( isShowDir )
	    printf("=======");
	  if ( whichDataType == PROMER_DATA )
	    printf("==================");
	}
      printf("===================================="
	     "=====================================\n");
    }

  //-- Print the stats
  for ( Sip = Stats.begin( ); Sip < Stats.end( ); Sip ++ )
    {
      if ( Sip->sA < 0 )
	continue;

      len1 = labs(Sip->eA - Sip->sA) + 1;
      len2 = labs(Sip->eB - Sip->sB) + 1;
      covA = (float)len1 / (float)Sip->SeqLenA * 100.0;
      covB = (float)len2 / (float)Sip->SeqLenB * 100.0;

      //-- Output the statistics for this alignment region
      printf("%8ld %8ld  | ", Sip->sA, Sip->eA);
      printf("%8ld %8ld  | ", Sip->sB, Sip->eB);
      printf("%8ld %8ld  | ", len1, len2);
      if ( !isBrief )
	{
	  printf("%8.2f ", Sip->Idy);
	  if ( whichDataType == PROMER_DATA )
	    printf("%8.2f %8.2f ", Sip->Sim, Sip->Stp);
	  printf(" | ");
	}
      if ( isShowSeqLens )
	printf("%8ld %8ld  | ", Sip->SeqLenA, Sip->SeqLenB);
      if ( isShowCoverage )
	printf("%8.2f %8.2f  | ", covA, covB);
      if ( !isBrief  &&  isShowDir )
	printf("%2d %2d  ", Sip->FrameA, Sip->FrameB);
      printf("%s\t%s", Sip->IdA, Sip->IdB);
      if ( isAnnotation )
	printf("\t%s", Sip->annot);

      printf("\n");
    }

  return;
}




void printTabular
     (vector<AlignStats> Stats)

     //  Print stats in a column delimited format

{
  long int len1, len2;
  float covA, covB;
  vector<AlignStats>::iterator Sip;

  //-- Print the output header
  if ( isPrintHeader )
    {
      printf ("%s %s\n%s\n\n", RefFileName, QryFileName,
	      whichDataType == NUCMER_DATA ? "NUCMER" : "PROMER");
      printf("%s\t%s\t", "[S1]", "[E1]");
      printf("%s\t%s\t", "[S2]", "[E2]");
      printf("%s\t%s\t", "[LEN 1]", "[LEN 2]");
      if ( !isBrief )
	{
	  printf("%s\t", "[% IDY]");
	  if ( whichDataType == PROMER_DATA )
	    printf("%s\t%s\t", "[% SIM]", "[% STP]");
	}
      if ( isShowSeqLens )
	printf("%s\t%s\t", "[LEN R]", "[LEN Q]");
      if ( isShowCoverage )
	printf("%s\t%s\t", "[COV R]", "[COV Q]");
      if ( !isBrief  &&  isShowDir )
	printf("%s\t", "[FRM]");
      printf("%s\n", "[TAGS]");
    }

  //-- Print the stats
  for ( Sip = Stats.begin( ); Sip < Stats.end( ); Sip ++ )
    {
      if ( Sip->sA < 0 )
	continue;

      len1 = labs(Sip->eA - Sip->sA) + 1;
      len2 = labs(Sip->eB - Sip->sB) + 1;
      covA = (float)len1 / (float)Sip->SeqLenA * 100.0;
      covB = (float)len2 / (float)Sip->SeqLenB * 100.0;

      //-- Output the statistics for this alignment region
      printf("%ld\t%ld\t", Sip->sA, Sip->eA);
      printf("%ld\t%ld\t", Sip->sB, Sip->eB);
      printf("%ld\t%ld\t", len1, len2);
      if ( !isBrief )
	{
	  printf("%.2f\t", Sip->Idy);
	  if ( whichDataType == PROMER_DATA )
	    printf("%.2f\t%.2f\t", Sip->Sim, Sip->Stp);
	}
      if ( isShowSeqLens )
	  printf("%ld\t%ld\t", Sip->SeqLenA, Sip->SeqLenB);
      if ( isShowCoverage )
	  printf("%.2f\t%.2f\t", covA, covB);
      if ( !isBrief  &&  isShowDir )
	printf("%d\t%d\t", Sip->FrameA, Sip->FrameB);
      printf("%s\t%s", Sip->IdA, Sip->IdB);
      if ( isAnnotation )
	printf("\t%s", Sip->annot);
      printf("\n");
    }

  return;
}




void parseDelta
     (vector<AlignStats> & Stats)

     //  parse the delta file

{
  int frameA, frameB;                  //  sequence frame
  long int sA, eA, sB, eB;
  AlignStats CurrStats;                //  single alignment region

  DeltaReader_t dr;
  dr.open (InputFileName);
  whichDataType = dr.getDataType( ) == NUCMER_STRING ?
    NUCMER_DATA : PROMER_DATA;
  strcpy (RefFileName, dr.getReferencePath( ).c_str( ));
  strcpy (QryFileName, dr.getQueryPath( ).c_str( ));

  //-- for each delta record
  while ( dr.readNext( ) )
    {
      CurrStats.SeqLenA = dr.getRecord( ).lenR;
      CurrStats.SeqLenB = dr.getRecord( ).lenQ;

      CurrStats.IdA = (char *) Safe_malloc (dr.getRecord( ).idR.length( ) + 1);
      CurrStats.IdB = (char *) Safe_malloc (dr.getRecord( ).idQ.length( ) + 1);
      strcpy (CurrStats.IdA, dr.getRecord( ).idR.c_str( ));
      strcpy (CurrStats.IdB, dr.getRecord( ).idQ.c_str( ));

      //-- for each alignment
      for ( unsigned int i = 0; i < dr.getRecord( ).aligns.size( ); i ++ )
	{
	  CurrStats.sA = dr.getRecord( ).aligns[i].sR;
	  CurrStats.eA = dr.getRecord( ).aligns[i].eR;
	  CurrStats.sB = dr.getRecord( ).aligns[i].sQ;
	  CurrStats.eB = dr.getRecord( ).aligns[i].eQ;
	  CurrStats.Idy = dr.getRecord( ).aligns[i].idy;
	  CurrStats.Sim = dr.getRecord( ).aligns[i].sim;
	  CurrStats.Stp = dr.getRecord( ).aligns[i].stp;

	  if ( CurrStats.Idy < idyCutoff )
	    continue;

	  if ( CurrStats.Idy > 99.99  &&  CurrStats.Idy != 100 )
	    CurrStats.Idy = 99.99;
	  if ( CurrStats.Sim > 99.99  &&  CurrStats.Sim != 100 )
	    CurrStats.Sim = 99.99;
	  if ( CurrStats.Stp > 99.99  &&  CurrStats.Stp != 100 )
	    CurrStats.Stp = 99.99;
	  if ( CurrStats.Stp < 0.01  &&  CurrStats.Stp != 0 )
	    CurrStats.Stp = 0.01;

	  if ( labs (CurrStats.sA - CurrStats.eA) + 1 < lenCutoff  ||
	       labs (CurrStats.sB - CurrStats.eB) + 1 < lenCutoff )
	    continue;

	  sA = CurrStats.sA;
	  eA = CurrStats.eA;
	  sB = CurrStats.sB;
	  eB = CurrStats.eB;

	  //-- Reset the coordinates to reference the appropriate strand
	  frameA = 1;
	  frameB = 1;
	  if ( sA > eA )
	    {
	      sA = revC (sA, CurrStats.SeqLenA);
	      eA = revC (eA, CurrStats.SeqLenA);
	      frameA += 3;
	    }
	  if ( sB > eB )
	    {
	      sB = revC (sB, CurrStats.SeqLenB);
	      eB = revC (eB, CurrStats.SeqLenB);
	      frameB += 3;
	    }

	  if ( isBrief )
	    {
	      CurrStats.FrameA = frameA > 3 ? (frameA - 3) * -1 : frameA;
	      CurrStats.FrameB = frameB > 3 ? (frameB - 3) * -1 : frameB;
	      Stats.push_back (CurrStats);
	      continue;
	    }

	  //-- Set/Generate the correct frame for sequences A and B
	  if ( whichDataType == NUCMER_DATA )
	    assert ( frameA == 1  &&  (frameB == 1 || frameB == 4));
	  else   // PROMER_DATA
	    {
	      //-- Set the reading frame
	      frameA += (sA + 2) % 3;
	      frameB += (sB + 2) % 3;

              //-- Translated the coordinates from DNA to Amino Acid
              //   remeber that eA and eB point to the last nucleotide in the
              //   end codon
              sA = (sA + 2) / 3;
              eA = (eA) / 3;
              sB = (sB + 2) / 3;
              eB = (eB) / 3;
	    }

	  //-- Set the statistics for this alignment region
	  CurrStats.FrameA = frameA > 3 ? (frameA - 3) * -1 : frameA;
	  CurrStats.FrameB = frameB > 3 ? (frameB - 3) * -1 : frameB;
	  
	  //-- Add the alignment region
	  Stats.push_back (CurrStats);
	}
    }
  dr.close ( );

  return;
}




inline long int revC
     (long int Coord, long int Len)

     //  Reverse complement the given coordinate for the given length.

{
  assert (Len - Coord + 1 > 0);
  return (Len - Coord + 1);
}




void simplifyAlignments
     (vector<AlignStats> & Stats)

     //  clean up alignments if the isBrief or isKnockout options were selected

{
  long int tmp;
  long int olapA, olapB, lenA, lenAi, lenB, lenBi;
  vector<AlignStats>::iterator Sp;
  vector<AlignStats>::iterator Sip;
  vector<AlignStats>::iterator Sep;
 
  Sp = Stats.begin( );
  while ( Sp < Stats.end( ) )
    {
      for ( Sep = Sp; Sep < Stats.end( ); Sep ++ )
	{
	  if ( strcmp (Sp->IdA, Sep->IdA) != 0  ||
	       strcmp (Sp->IdB, Sep->IdB) != 0 )
	    break;
	  
	  if ( Sep->FrameA < 0 )
	    {
	      tmp = Sep->sA;
	      Sep->sA = Sep->eA;
	      Sep->eA = tmp;
	    }
	  if ( Sep->FrameB < 0 )
	    {
	      tmp = Sep->sB;
	      Sep->sB = Sep->eB;
	      Sep->eB = tmp;
	    }

	  if ( isBrief )
	    {
	      Sep->FrameA = 1;
	      Sep->FrameB = 1;
	    }
	}
      
      sort (Sp, Sep, sA_Sort( ));
      
      for ( ; Sp < Sep; Sp ++ )
	{
	  if ( Sp->sA < 0 )
	    continue;
	  
	  for ( Sip = Sp + 1; Sip < Sep; Sip ++ )
	    {
	      if ( Sip->sA > Sp->eA )
		break;
	      
	      if ( isBrief )
		{
		  if ( Sip->sB <= Sp->eB  &&  Sip->eB >= Sp->sB )
		    {
		      Sp->sA = Min ( Sp->sA, Sip->sA );
		      Sp->sB = Min ( Sp->sB, Sip->sB );
		      Sp->eA = Max ( Sp->eA, Sip->eA );
		      Sp->eB = Max ( Sp->eB, Sip->eB );
		      Sip->sA *= -1;
		    }
		}
	      else if ( isKnockout )
		{
		  //-- Are they on the same strand?
		  if ( (Sp->FrameA * Sp->FrameB) *
		       (Sip->FrameA * Sip->FrameB) < 0 )
		    continue;

		  //-- If overlap by more than 50% of their length
		  //   pick the best frame and throw out the loser
		  lenA = Sp->eA - Sp->sA + 1;
		  lenAi = Sip->eA - Sip->sA + 1;
		  lenB = Sp->eB - Sp->sB + 1;
		  lenBi = Sip->eB - Sip->sB + 1;

		  olapA = Sp->eA - Sip->sA + 1;
		  olapB = Sp->eB - Sp->sB + 1;
		  tmp = Sp->eB - Sip->eB;
		  if ( tmp > 0 )
		    olapB -= tmp;
		  tmp = Sip->sB - Sp->sB;
		  if ( tmp > 0 )
		    olapB -= tmp;
		  if ( (float)olapA /
		       (float)(lenA < lenAi ? lenA : lenAi) >= MIN_OVERLAP_P &&
		       (float)olapB /
		       (float)(lenB < lenBi ? lenB : lenBi) >= MIN_OVERLAP_P )
		    {
		      if ( Sp->Sim > Sip->Sim &&
			   ((float)lenA / (float)lenAi >= MIN_LENGTH_P &&
			    (float)lenB / (float)lenBi >= MIN_LENGTH_P) )
			Sip->sA *= -1;
		      else
			{
			  Sp->sA *= -1;
			  break;
			}
		    }
		}
	      else
		assert (false);
	    }

	  if ( Sp->sA > 0 && Sp->FrameA < 0 )
	    {
	      tmp = Sp->sA;
	      Sp->sA = Sp->eA;
	      Sp->eA = tmp;
	    }
	  if ( Sp->sA > 0 && Sp->FrameB < 0 )
	    {
	      tmp = Sp->sB;
	      Sp->sB = Sp->eB;
	      Sp->eB = tmp;
	    }
	}
    }

  return;
}




void printHelp
     (const char * s)

     //  Display the program's help information to stderr.

{
  fprintf (stderr,
	   "\nUSAGE: %s  [options]  <deltafile>\n\n", s);
  fprintf (stderr,
       "-b          Merges overlapping alignments regardless of match dir\n"
       "            or frame and does not display any idenitity information.\n"
       "-B          Switch output to btab format\n"
       "-c          Include percent coverage information in the output\n"
       "-d          Display the alignment direction in the additional\n"
       "            FRM columns (default for promer)\n"
       "-g          Deprecated option. Please use 'delta-filter' instead\n"
       "-h          Display help information\n"
       "-H          Do not print the output header\n"
       "-I float    Set minimum percent identity to display\n"
       "-k          Knockout (do not display) alignments that overlap\n"
       "            another alignment in a different frame by more than 50%%\n"
       "            of their length, AND have a smaller percent similarity\n"
       "            or are less than 75%% of the size of the other alignment\n"
       "            (promer only)\n"
       "-l          Include the sequence length information in the output\n"
       "-L long     Set minimum alignment length to display\n"
       "-o          Annotate maximal alignments between two sequences, i.e.\n"
       "            overlaps between reference and query sequences\n"
       "-q          Sort output lines by query IDs and coordinates\n"
       "-r          Sort output lines by reference IDs and coordinates\n"
       "-T          Switch output to tab-delimited format\n\n");
  fprintf (stderr,
	   "  Input is the .delta output of either the \"nucmer\" or the\n"
	   "\"promer\" program passed on the command line.\n"
	   "  Output is to stdout, and consists of a list of coordinates,\n"
	   "percent identity, and other useful information regarding the\n"
	   "alignment data contained in the .delta file used as input.\n"
	   "  NOTE: No sorting is done by default, therefore the alignments\n"
	   "will be ordered as found in the <deltafile> input.\n\n");
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
