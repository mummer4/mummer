#!@PERL@

#-------------------------------------------------------------------------------
#   Programmer: Adam M Phillippy, The Institute for Genomic Research
#         File: promer
#         Date: 04 / 09 / 03
#
#        Usage:
#    promer  [options]  <Reference>  <Query>
#
#                Try 'promer -h' for more information.
#
#      Purpose: To create alignments between two multi-FASTA inputs by using
#              the MUMmer matching and clustering algorithms.
#
#-------------------------------------------------------------------------------

use lib "@LIB_DIR@";
use Foundation;
use File::Spec::Functions;
use strict;

my $BIN_DIR = "@BIN_DIR@";
my $AUX_BIN_DIR = "@LIBEXEC_DIR@";
my $LIB_DIR = "@LIB_DIR@";



my $VERSION_INFO = q~
PROmer (PROtein MUMmer) version 3.07
    ~;



my $HELP_INFO = q~
  USAGE: promer  [options]  <Reference>  <Query>

  DESCRIPTION:
    promer generates amino acid alignments between two mutli-FASTA DNA input
    files. The out.delta output file lists the distance between insertions
    and deletions that produce maximal scoring alignments between each
    sequence. The show-* utilities know how to read this format. The DNA
    input is translated into all 6 reading frames in order to generate the
    output, but the output coordinates reference the original DNA input.

  MANDATORY:
    Reference       Set the input reference multi-FASTA DNA file
    Query           Set the input query multi-FASTA DNA file

  OPTIONS:
    --mum           Use anchor matches that are unique in both the reference
                    and query
    --mumcand       Same as --mumreference
    --mumreference  Use anchor matches that are unique in in the reference
                    but not necessarily unique in the query (default behavior)
    --maxmatch      Use all anchor matches regardless of their uniqueness

    -b|breaklen     Set the distance an alignment extension will attempt to
                    extend poor scoring regions before giving up, measured in
                    amino acids (default 60)
    -c|mincluster   Sets the minimum length of a cluster of matches, measured in
                    amino acids (default 20)
    --[no]delta     Toggle the creation of the delta file (default --delta)
    --depend        Print the dependency information and exit
    -d|diagfactor   Set the clustering diagonal difference separation factor
                    (default .11)
    --[no]extend    Toggle the cluster extension step (default --extend)
    -g|maxgap       Set the maximum gap between two adjacent matches in a
                    cluster, measured in amino acids (default 30)
    -h
    --help          Display help information and exit.
    -l|minmatch     Set the minimum length of a single match, measured in amino
                    acids (default 6)
    -m|masklen      Set the maximum bookend masking lenth, measured in amino
                    acids (default 8)
    -o
    --coords        Automatically generate the original PROmer1.1 ".coords"
                    output file using the "show-coords" program
    --[no]optimize  Toggle alignment score optimization, i.e. if an alignment
                    extension reaches the end of a sequence, it will backtrack
                    to optimize the alignment score instead of terminating the
                    alignment at the end of the sequence (default --optimize)

    -p|prefix       Set the prefix of the output files (default "out")
    -V
    --version       Display the version information and exit
    -x|matrix       Set the alignment matrix number to 1 [BLOSUM 45], 2 [BLOSUM
                    62] or 3 [BLOSUM 80] (default 2)
    ~;


my $USAGE_INFO = q~
  USAGE: promer  [options]  <Reference>  <Query>
    ~;


my @DEPEND_INFO =
    (
     "$BIN_DIR/mummer",
     "$AUX_BIN_DIR/mgaps",
     "$BIN_DIR/show-coords",
     "$AUX_BIN_DIR/postpro",
     "$AUX_BIN_DIR/prepro",
     "$LIB_DIR/Foundation.pm"
     );


my %DEFAULT_PARAMETERS =
    (
     "OUTPUT_PREFIX"     =>   "out",      # prefix for all output files
     "MATCH_ALGORITHM"   =>   "-mumreference", # match finding algo switch
     "MIN_MATCH"         =>   "6",        # minimum match size (aminos)
     "MAX_GAP"           =>   "30",       # maximum gap between matches (aminos)
     "MIN_CLUSTER"       =>   "20",       # minimum cluster size (aminos)
     "DIAG_FACTOR"       =>   ".11",      # diagonal difference fraction
     "BREAK_LEN"         =>   "60",       # extension break length
     "BLOSUM_NUMBER"     =>   "2",        # options are 1,2,3 (BLOSUM 45,62,80)
     "MASKING_LENGTH"    =>   "8",        # set bookend masking length
     "POST_SWITCHES"     =>   ""          # switches for the post processing
     );


sub main ( )
{
    my $tigr;             # TIGR::Foundation object
    my @err;              # Error variable
    
    my $ref_file;         # path of the reference input file
    my $qry_file;         # path of the query input file

    #-- The command line options for the various programs
    my $pfx = $DEFAULT_PARAMETERS { "OUTPUT_PREFIX" };
    my $algo = $DEFAULT_PARAMETERS { "MATCH_ALGORITHM" };
    my $size = $DEFAULT_PARAMETERS { "MIN_MATCH" };
    my $gap = $DEFAULT_PARAMETERS { "MAX_GAP" };
    my $clus = $DEFAULT_PARAMETERS { "MIN_CLUSTER" };
    my $diff = $DEFAULT_PARAMETERS { "DIAG_FACTOR" };
    my $blen = $DEFAULT_PARAMETERS { "BREAK_LEN" };
    my $blsm = $DEFAULT_PARAMETERS { "BLOSUM_NUMBER" };
    my $mask = $DEFAULT_PARAMETERS { "MASKING_LENGTH" };
    my $psw = $DEFAULT_PARAMETERS { "POST_SWITCHES" };

    my $maxmatch;         # matching algorithm switches
    my $mumreference;
    my $mum;
    my $extend = 1;       # if true, extend clusters
    my $delta = 1;        # if true, create the delta file
    my $optimize = 1;     # if true, optimize alignment scores

    my $generate_coords;

    #-- Initialize TIGR::Foundation
    $tigr = new TIGR::Foundation;
    if ( !defined ($tigr) ) {
	print (STDERR "ERROR: TIGR::Foundation could not be initialized");
	exit (1);
    }
    
    #-- Set help and usage information
    $tigr->setHelpInfo ($HELP_INFO);
    $tigr->setUsageInfo ($USAGE_INFO);
    $tigr->setVersionInfo ($VERSION_INFO);
    $tigr->addDependInfo (@DEPEND_INFO);

    #-- Get command line parameters
    $err[0] = $tigr->TIGR_GetOptions
	(
         "maxmatch" => \$maxmatch,
	 "mumcand" => \$mumreference,
         "mumreference" => \$mumreference,
         "mum" => \$mum,
	 "b|breaklen=i" => \$blen,
	 "c|mincluster=i" => \$clus,
	 "delta!" => \$delta,
	 "d|diagfactor=f" => \$diff,
	 "extend!" => \$extend,
	 "g|maxgap=i" => \$gap,
	 "l|minmatch=i" => \$size,
	 "m|masklen=i" => \$mask,
	 "o|coords"   => \$generate_coords,
	 "optimize!" => \$optimize,
	 "p|prefix=s" => \$pfx,
	 "x|matrix=i" => \$blsm
	 );

    #-- Check if the parsing was successful
    if ( $err[0] == 0  ||  $#ARGV != 1 ) {
	$tigr->printUsageInfo( );
	print (STDERR "Try '$0 -h' for more information.\n");
	exit (1);
    }

    $ref_file = File::Spec->rel2abs ($ARGV[0]);
    $qry_file = File::Spec->rel2abs ($ARGV[1]);
    
    #-- Set up the program parameters
    if ( ! $extend ) {
	$psw .= "-e ";
    }
    if ( ! $delta ) {
	$psw .= "-d ";
    }
    if ( ! $optimize ) {
	$psw .= "-t ";
    }

    undef (@err);
    $err[0] = 0;
    if ( $mum ) {
	$err[0] ++;
	$algo = "-mum";
    }
    if ( $mumreference ) {
	$err[0] ++;
	$algo = "-mumreference";
    }
    if ( $maxmatch ) {
	$err[0] ++;
	$algo = "-maxmatch";
    }
    if ( $err[0] > 1 ) {
	$tigr->printUsageInfo( );
	print (STDERR "ERROR: Multiple matching algorithms selected\n");
	print (STDERR "Try '$0 -h' for more information.\n");
	exit (1);
    }

    #-- Set up the program path names
    my $algo_path = "$BIN_DIR/mummer";
    my $mgaps_path = "$AUX_BIN_DIR/mgaps";
    my $prepro_path = "$AUX_BIN_DIR/prepro";
    my $postpro_path = "$AUX_BIN_DIR/postpro";
    my $showcoords_path = "$BIN_DIR/show-coords";
		     
    #-- Check that the files needed are all there and readable/writable
    {
	undef (@err);
	if ( !$tigr->isExecutableFile ($algo_path) ) {
	    push (@err, $algo_path);
	}
	
	if ( !$tigr->isExecutableFile ($mgaps_path) ) {
	    push (@err, $mgaps_path);
	}
	
	if ( !$tigr->isExecutableFile ($prepro_path) ) {
	    push (@err, $prepro_path);
	}
	
	if ( !$tigr->isExecutableFile ($postpro_path) ) {
	    push (@err, $postpro_path);
	}
	
	if ( !$tigr->isReadableFile ($ref_file) ) {
	    push (@err, $ref_file);
	}
	
	if ( !$tigr->isReadableFile ($qry_file) ) {
	    push (@err, $qry_file);
	}
	
	if ( !$tigr->isCreatableFile ("$pfx.aaref") ) {
	    if ( !$tigr->isWritableFile ("$pfx.aaref") ) {
		push (@err, "$pfx.aaref");
	    }
	}

	if ( !$tigr->isCreatableFile ("$pfx.aaqry") ) {
	    if ( !$tigr->isWritableFile ("$pfx.aaqry") ) {
		push (@err, "$pfx.aaqry");
	    }
	}
	
	if ( !$tigr->isCreatableFile ("$pfx.mgaps") ) {
	    if ( !$tigr->isWritableFile ("$pfx.mgaps") ) {
		push (@err, "$pfx.mgaps");
	    }
	}
	
	if ( !$tigr->isCreatableFile ("$pfx.delta") ) {
	    if ( !$tigr->isWritableFile ("$pfx.delta") ) {
		push (@err, "$pfx.delta");
	    }
	}

	if ( $generate_coords ) {
	    if ( !$tigr->isExecutableFile ($showcoords_path) ) {
		push (@err, $showcoords_path);
	    }
	    if ( !$tigr->isCreatableFile ("$pfx.coords") ) {
		if ( !$tigr->isWritableFile ("$pfx.coords") ) {
		    push (@err, "$pfx.coords");
		}
	    }
	}
    
	#-- If 1 or more files could not be processed, terminate script
	if ( $#err >= 0 ) {
	    $tigr->logError
		("ERROR: The following critical files could not be used", 1);
	    while ( $#err >= 0 ) {
		$tigr->logError (pop(@err), 1);
	    }
	    $tigr->logError
		("Check your paths and file permissions and try again", 1);
	    $tigr->bail( );
	}
    }
    

    #-- Run prepro -r and -q and assert return value is zero
    print (STDERR "1: PREPARING DATA\n");
    $err[0] = $tigr->runCommand
	("$prepro_path -m $mask -r $ref_file > $pfx.aaref");

    if ( $err[0] != 0 ) {
	$tigr->bail
	    ("ERROR: prepro -r returned non-zero\n");
    }

    $err[0] = $tigr->runCommand
	("$prepro_path -m $mask -q $qry_file > $pfx.aaqry");

    if ( $err[0] != 0 ) {
	$tigr->bail ("ERROR: prepro -q returned non-zero\n");
    }


    #-- Run mummer | mgaps and assert return value is zero
    print (STDERR "2,3: RUNNING mummer AND CREATING CLUSTERS\n");
    open(ALGO_PIPE, "$algo_path $algo -l $size $pfx.aaref $pfx.aaqry |")
	or $tigr->bail ("ERROR: could not open $algo_path output pipe $!");
    open(CLUS_PIPE, "| $mgaps_path -l $clus -s $gap -f $diff > $pfx.mgaps")
	or $tigr->bail ("ERROR: could not open $mgaps_path input pipe $!");
    while ( <ALGO_PIPE> ) {
	print CLUS_PIPE
	or $tigr->bail ("ERROR: could not write to $mgaps_path pipe $!");
    }
    $err[0] = close(ALGO_PIPE);
    $err[1] = close(CLUS_PIPE);

    if ( $err[0] == 0  ||  $err[1] == 0 ) {
	$tigr->bail ("ERROR: mummer and/or mgaps returned non-zero\n");
    }


    #-- Run postpro and assert return value is zero
    print (STDERR "4: FINISHING DATA\n");
    $err[0] = $tigr->runCommand
	("$postpro_path $psw -x $blsm -b $blen ".
	 "$ref_file $qry_file $pfx < $pfx.mgaps");

    if ( $err[0] != 0 ) {
	$tigr->bail ("ERROR: postpro returned non-zero\n");
    }

    #-- If the -o flag was set, run show-coords using PROmer1.1 settings
    if ( $generate_coords ) {
	print (STDERR "5: GENERATING COORDS FILE\n");
	$err[0] = $tigr->runCommand
	    ("$showcoords_path -r $pfx.delta > $pfx.coords");
	
	if ( $err[0] != 0 ) {
	    $tigr->bail ("ERROR: show-coords returned non-zero\n");
	}
    }
 
    #-- Remove the temporary output
    $err[0] = unlink ("$pfx.aaref", "$pfx.aaqry", "$pfx.mgaps");
 
    if ( $err[0] != 3 ) {
 	$tigr->logError ("WARNING: there was a problem deleting".
			 " the temporary output files", 1);
    }
 
    #-- Return success
    return (0);
}

exit ( main ( ) );

#-- END OF SCRIPT
