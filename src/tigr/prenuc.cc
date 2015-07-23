//------------------------------------------------------------------------------
//   Programmer: Adam M Phillippy, The Institute for Genomic Research
//         File: prenuc.cc
//         Date: 07 / 15 / 2002
//
//        Input: Input is a single multi-fasta sequence file
//
//       Output: Output is to stdout, and it consists of each sequence in the
//              FASTA file appended together with all the headers removed. A
//             new generic header is inserted at the beginning of the file to
//            adhere to FASTA standards.
//
//       NOTICE: An `x' is placed at the end of all sequences so that no MUMs
//              will span across two different sequences. This extra character
//             is not included in the sequence lengths.
//
//        Usage: prenuc <reference>
//
//------------------------------------------------------------------------------

#include <mummer/tigrinc.hh>

//-- Output this many sequence characters per line
#define CHARS_PER_LINE 60

void printHelp
     (const char *);

void printUsage
     (const char *);

int main
     (int argc, char * argv[])
{
  long int InitSize, SizeA, ct, i;
  char * A, Id[MAX_LINE];
  FILE * RefFile;

  //-- Parse the command line arguments
  {
    optarg = NULL;
    int ch, errflg = 0;
    while ( !errflg  &&  ((ch = getopt (argc, argv, "hT:qrclw")) != EOF) )
      switch (ch)
	{
	case 'h' :
	  printHelp (argv[0]);
	  exit (EXIT_SUCCESS);
	  break;

	default :
	  errflg ++;
	}
    
    if ( errflg > 0 || argc - optind != 1 )
      {
	printUsage (argv[0]);
	exit (EXIT_FAILURE);
      }
  }

  RefFile = File_Open(argv[1],"r");

  InitSize = INIT_SIZE;
  A = (char *) Safe_malloc (sizeof(char) * INIT_SIZE);

  ct = 0;

  printf(">allcontigs %s\n", argv[1]);
  while ( Read_String (RefFile, A, InitSize, Id, false) )
    {
      SizeA = strlen(A + 1);

      //-- Output the concatenated sequence
      for ( i = 1; i <= SizeA; i ++ )
        {
          fputc (A[i], stdout);
          if ( ++ ct == CHARS_PER_LINE )
            {
              ct = 0;
              fputc ('\n', stdout);
            }
        }

      fputc ('x', stdout);
      if ( ++ ct == CHARS_PER_LINE )
	{
	  ct = 0;
	  fputc ('\n', stdout);
	}
    }
  if ( ct != 0 )
    fputc ('\n', stdout);

  fclose(RefFile);

  free(A);

  return EXIT_SUCCESS;
}



void printHelp
     (const char * s)
{
  fprintf(stderr,
      "\nUSAGE: %s  [options]  <reference>\n\n", s);
  fprintf(stderr,"-h     display help information\n\n");
  fprintf(stderr,
      "  Input is one multi-fasta sequence file.\n"
      "  Output is to stdout, and it consists of each sequence in the\n"
      "FASTA file appended together with all the headers removed. A\n"
      "new generic header is inserted at the beginning of the file to\n"
      "adhere to FASTA standards. An `x' is placed at the end of all\n"
      "sequences so that no MUMs will span two different sequences.\n\n");
  return;
}



void printUsage
     (const char * s)
{
  fprintf(stderr,
      "\nUSAGE: %s  [options]  <reference>\n\n", s);
  fprintf (stderr, "Try '%s -h' for more information.\n", s);
  return;
}
