//------------------------------------------------------------------------------
//   Programmer: Adam M Phillippy, The Institute for Genomic Research
//         File: prepro.cc
//         Date: 08 / 16 / 2002
//
//        Input: Input is a single multi-FASTA sequence file on the command
//              line, the command line switch '-r' specifies that the input
//             is the reference sequence. This switch tells the program to
//            suppress all but the first header and use a different masking
//           character than the query sequence. If it is the query sequence
//          being input, use the '-q' option. Either -r or -q must be specified.
//
//       Output: Output is to stdout, and consists of each sequence in the input
//              translated into all six reading frames. All the translations
//             for a particular sequence are appended together and each
//            seperated by a masking character. If the '-r' flag is specified,
//           all but the top header are removed from the output, if the '-q'
//          flag is specified each sequence header is left in place but the
//         translated frames remain concatenated together without headers.
//               In addition to translating the sequences, prepro also does
//              a simple masking step that masks amino acids that are
//             surrounded on either side with stop codons. The '-m len' option
//            allows the user to specify the maximum length of these regions
//           to be masked. For example, if '-m 5' were set and the sequence
//          "A*AAAAA*A" appeared, the output would be "A*XXXXX*A". The 0 and
//         Len+1 indices (although non-existent) are considered as stop codons
//        so this masking with -m 5 could turn the beginning of a sequence
//       "\0AAAAA*" into "\0XXXXX*" or the end "*AAAAA\0" into "*XXXXX\0".
//
//        Usage: prepro [-m len] -r/-q <multi-FASTA>
//
//------------------------------------------------------------------------------

#include <mummer/tigrinc.hh>
#include <mummer/translate.hh>

//-- Output this many sequence characters per line
#define CHARS_PER_LINE 60

const long int DEFAULT_MASK_LEN = 10;

const char TRANSLATE_MASK = 'X';     // translator masking character
const char REFERENCE_MASK = 'X';     // masking character to use for reference
const char QUERY_MASK = 'O';         // masking character to use for query
const char STOP_MASK = 'J';          // alpha character for stop codons
const char STOP_CHAR = '*';          // stop codon character

inline void mask
     (char * A, char mask_ch, long int x, long int y);

void printHelp
     (const char *);

void printUsage
     (const char *);




int main
     (int argc, char * argv[])
{
  bool isReference = false;
  bool isQuery = false;

  int frame;
  int mask_len = DEFAULT_MASK_LEN;
 
  long int InitSize, LenA, LentA, ct, i;
  long int last_index;

  char * A, * tA;
  char mask_char = 0;
  char Id [MAX_LINE];
  char InputFileName [MAX_LINE];

  FILE * InputFile;

  //-- Parse the command line arguments
  {
    optarg = NULL;
    int ch, errflg = 0;
    while ( !errflg  &&  ((ch = getopt (argc, argv, "hm:q:r:")) != EOF) )
      switch (ch)
        {
        case 'h' :
          printHelp (argv[0]);
          exit (EXIT_SUCCESS);
          break;

	case 'm' :
	  mask_len = atoi (optarg);
	  break;

	case 'q' :
	  strcpy (InputFileName, optarg);
	  isQuery = true;
	  mask_char = QUERY_MASK;
	  break;

	case 'r' :
	  strcpy (InputFileName, optarg);
	  isReference = true;
	  mask_char = REFERENCE_MASK;
	  break;

        default :
          errflg ++;
        }
    
    if ( errflg > 0 || argc - optind != 0 ||
	 (isReference  &&  isQuery)  ||  (!isReference  &&  !isQuery) )
      {
        printUsage (argv[0]);
        exit (EXIT_FAILURE);
      }

    if ( mask_len < 0 )
      {
	fprintf (stderr,
	  "WARNING: Invalid maximum mask length %d, ignored\n", mask_len);
	mask_len = DEFAULT_MASK_LEN;
      }
  }

  InputFile = File_Open (InputFileName, "r");

  InitSize = INIT_SIZE;
  A = (char *) Safe_malloc ( sizeof(char) * InitSize );
  tA = (char *) Safe_malloc ( sizeof(char) );
  tA [0] = '\0';
  
  ct = 0;
  if ( isReference )
    printf (">allcontigs %s\n", InputFileName);
  while ( Read_String (InputFile, A, InitSize, Id, false) )
    {
      LenA = strlen(A + 1);

      for ( frame = 1; frame <= 6; frame ++ )
	{
	  if ( isQuery )
	    printf (">%s.%d\n", Id, frame);

	  //-- Translate the current frame
	  tA = (char *) Safe_realloc (tA, sizeof(char) * ( (LenA / 3) + 2) );
	  LentA = Translate_DNA (A, tA, frame);
	  tA[++ LentA] = mask_char;
	  
	  //-- Mask the current frame
	  last_index = 0;
	  for ( i = 1; i <= LentA; i ++ )
	    {
	      if ( mask_char != TRANSLATE_MASK && tA[i] == TRANSLATE_MASK )
		tA[i] = mask_char;
	      else if ( tA[i] == STOP_CHAR )
		{
		  tA[i] = STOP_MASK;
		  if ( i - last_index - 1 <= mask_len )
		    mask (tA, mask_char, last_index + 1, i - 1);
		  last_index = i;
		}
	    }
	  if ( LentA - last_index - 1 <= mask_len )
	    mask (tA, mask_char, last_index + 1, i - 1);
	  
	  //-- Print the current frame
	  for ( i = 1; i <= LentA; i ++ )
	    {
	      fputc (tA[i], stdout);
	      if ( ++ ct == CHARS_PER_LINE )
		{
		  ct = 0;
		  fputc ('\n', stdout);
		}
	    }

	  if ( isQuery )
	    {
	      if ( ct != 0 )
		fputc ('\n', stdout);
	      ct = 0;
	    }
	}
    }
  if ( ct != 0 )
    fputc ('\n', stdout);

  fclose(InputFile);

  free(A);
  free(tA);

  return EXIT_SUCCESS;
}




inline void mask
     (char * A, char mask_ch, long int x, long int y)

  //  Mask sequence 'A' with 'mask_ch' from A [x...y] (inclusive)

{
  for ( ; x <= y; x ++ )
    A[x] = mask_ch;
}




void printHelp
     (const char * s)
{
  fprintf(stderr,
      "\nUSAGE: %s  [options]  -r/-q <fasta>\n\n", s);
  fprintf(stderr,
      "-h            display help information\n"
      "-m len        set maximum book-end masking length to 'len\n" 
      "-q query      input is the multi-fasta query file 'query'\n"
      "-r reference  input is the multi-fasta reference file 'reference'\n\n"
      "  Input is one multi-fasta sequence file, EITHER '-r reference' OR\n"
      "'-q query'. Both are not allowed.\n"
      "  Output is to stdout, and it consists of each sequence in the\n"
      "FASTA file translated in all six reading frames. This output is\n"
      "different depending on whether the the input was the reference\n"
      "or query sequence, and it is now ready to be passed to 'mummer2'\n"
      "for the match finding step.\n\n");
  return;
}




void printUsage
     (const char * s)
{
  fprintf(stderr,
      "\nUSAGE: %s  [options]  -r/-q <fasta>\n\n", s);
  fprintf (stderr, "Try '%s -h' for more information.\n", s);
  return;
}
