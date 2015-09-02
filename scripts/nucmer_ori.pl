#!__PERL_PATH

#-------------------------------------------------------------------------------
#   Programmer: Adam M Phillippy, The Institute for Genomic Research
#         File: nucmer
#         Date: 04 / 09 / 03
#
#        Usage:
#    nucmer  [options]  <Reference>  <Query>
#
#                Try 'nucmer -h' for more information.
#
#      Purpose: To create alignments between two multi-FASTA inputs by using
#              the MUMmer matching and clustering algorithms.
#
#-------------------------------------------------------------------------------

use lib "__SCRIPT_DIR";
use Foundation;
use File::Spec::Functions;
use strict;

my $AUX_BIN_DIR = "__AUX_BIN_DIR";
my $BIN_DIR = "__BIN_DIR";
my $SCRIPT_DIR = "__SCRIPT_DIR";


my $VERSION_INFO = q~
NUCmer (NUCleotide MUMmer) version 3.1
    ~;


my $HELP_INFO = q~
  USAGE: nucmer  [options]  <Reference>  <Query>

  DESCRIPTION:
    nucmer generates nucleotide alignments between two mutli-FASTA input
    files. The out.delta output file lists the distance between insertions
    and deletions that produce maximal scoring alignments between each
    sequence. The show-* utilities know how to read this format.

  MANDATORY:
    Reference       Set the input reference multi-FASTA filename
    Query           Set the input query multi-FASTA filename

  OPTIONS:
    --mum           Use anchor matches that are unique in both the reference
                    and query
    --mumcand       Same as --mumreference
    --mumreference  Use anchor matches that are unique in in the reference
                    but not necessarily unique in the query (default behavior)
    --maxmatch      Use all anchor matches regardless of their uniqueness

    -b|breaklen     Set the distance an alignment extension will attempt to
                    extend poor scoring regions before giving up (default 200)
    --[no]banded    Enforce absolute banding of dynamic programming matrix
                    based on diagdiff parameter EXPERIMENTAL (default no)
    -c|mincluster   Sets the minimum length of a cluster of matches (default 65)
    --[no]delta     Toggle the creation of the delta file (default --delta)
    --depend        Print the dependency information and exit
    -D|diagdiff     Set the maximum diagonal difference between two adjacent
                    anchors in a cluster (default 5)
    -d|diagfactor   Set the maximum diagonal difference between two adjacent
                    anchors in a cluster as a differential fraction of the gap
                    length (default 0.12)
    --[no]extend    Toggle the cluster extension step (default --extend)
    -f
    --forward       Use only the forward strand of the Query sequences
    -g|maxgap       Set the maximum gap between two adjacent matches in a
                    cluster (default 90)
    -h
    --help          Display help information and exit
    -l|minmatch     Set the minimum length of a single match (default 20)
    -o
    --coords        Automatically generate the original NUCmer1.1 coords
                    output file using the 'show-coords' program
    --[no]optimize  Toggle alignment score optimization, i.e. if an alignment
                    extension reaches the end of a sequence, it will backtrack
                    to optimize the alignment score instead of terminating the
                    alignment at the end of the sequence (default --optimize)
    -p|prefix       Set the prefix of the output files (default "out")
    -r
    --reverse       Use only the reverse complement of the Query sequences
    --[no]simplify  Simplify alignments by removing shadowed clusters. Turn
                    this option off if aligning a sequence to itself to look
                    for repeats (default --simplify)
    -V
    --version       Display the version information and exit
    ~;


my $USAGE_INFO = q~
  USAGE: nucmer  [options]  <Reference>  <Query>
    ~;


my @DEPEND_INFO =
    (
     "$BIN_DIR/mummer",
     "$BIN_DIR/mgaps",
     "$BIN_DIR/show-coords",
     "$AUX_BIN_DIR/postnuc",
     "$AUX_BIN_DIR/prenuc",
     "$SCRIPT_DIR/Foundation.pm"
     );


my %DEFAULT_PARAMETERS =
    (
     "OUTPUT_PREFIX"     =>   "out",        # prefix for all output files
     "MATCH_ALGORITHM"   =>   "-mumreference", # match finding algo switch
     "MATCH_DIRECTION"   =>   "-b",         # match direction switch
     "MIN_MATCH"         =>   "20",         # minimum match size
     "MAX_GAP"           =>   "90",         # maximum gap between matches
     "MIN_CLUSTER"       =>   "65",         # minimum cluster size
     "DIAG_DIFF"         =>   "5",          # diagonal difference absolute
     "DIAG_FACTOR"       =>   ".12",        # diagonal difference fraction
     "BREAK_LEN"         =>   "200",        # extension break length
     "POST_SWITCHES"     =>   ""            # switches for the post processing
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
    my $mdir = $DEFAULT_PARAMETERS { "MATCH_DIRECTION" };
    my $size = $DEFAULT_PARAMETERS { "MIN_MATCH" };
    my $gap = $DEFAULT_PARAMETERS { "MAX_GAP" };
    my $clus = $DEFAULT_PARAMETERS { "MIN_CLUSTER" };
    my $ddiff = $DEFAULT_PARAMETERS { "DIAG_DIFF" };
    my $dfrac = $DEFAULT_PARAMETERS { "DIAG_FACTOR" };
    my $blen = $DEFAULT_PARAMETERS { "BREAK_LEN" };
    my $psw = $DEFAULT_PARAMETERS { "POST_SWITCHES" };

    my $fwd;              # if true, use forward strand
    my $rev;              # if true, use reverse strand
    my $maxmatch;         # matching algorithm switches
    my $mumreference;
    my $mum;
    my $banded = 0;       # if true, enforce absolute dp banding
    my $extend = 1;       # if true, extend clusters
    my $delta = 1;        # if true, create the delta file
    my $optimize = 1;     # if true, optimize alignment scores
    my $simplify = 1;     # if true, simplify shadowed alignments

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
         "banded!" => \$banded,
	 "c|mincluster=i" => \$clus,
	 "delta!" => \$delta,
         "D|diagdiff=i" => \$ddiff,
	 "d|diagfactor=f" => \$dfrac,
	 "extend!" => \$extend,
	 "f|forward"   => \$fwd,
	 "g|maxgap=i" => \$gap,
	 "l|minmatch=i" => \$size,
	 "o|coords"   => \$generate_coords,
	 "optimize!" => \$optimize,
	 "p|prefix=s" => \$pfx,
	 "r|reverse"   => \$rev,
	 "simplify!" => \$simplify
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
    if ( $fwd  &&  $rev ) {
	$mdir = "-b";
    } elsif ( $fwd ) {
	$mdir = "";
    } elsif ( $rev ) {
	$mdir = "-r";
    }
    if ( ! $extend ) {
	$psw .= "-e ";
    }
    if ( ! $delta ) {
	$psw .= "-d ";
    }
    if ( ! $optimize ) {
	$psw .= "-t ";
    }
    if ( ! $simplify ) {
	$psw .= "-s ";
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
    my $mgaps_path = "$BIN_DIR/mgaps";
    my $prenuc_path = "$AUX_BIN_DIR/prenuc";
    my $postnuc_path = "$AUX_BIN_DIR/postnuc";
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
	
	if ( !$tigr->isExecutableFile ($prenuc_path) ) {
	    push (@err, $prenuc_path);
	}
	
	if ( !$tigr->isExecutableFile ($postnuc_path) ) {
	    push (@err, $postnuc_path);
	}
	
	if ( !$tigr->isReadableFile ($ref_file) ) {
	    push (@err, $ref_file);
	}
	
	if ( !$tigr->isReadableFile ($qry_file) ) {
	    push (@err, $qry_file);
	}
	
	if ( !$tigr->isCreatableFile ("$pfx.ntref") ) {
	    if ( !$tigr->isWritableFile ("$pfx.ntref") ) {
		push (@err, "$pfx.ntref");
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
    

    #-- Run prenuc and assert return value is zero
    print (STDERR "1: PREPARING DATA\n");
    $err[0] = $tigr->runCommand
	("$prenuc_path $ref_file > $pfx.ntref");

    if ( $err[0] != 0 ) {
	$tigr->bail
	    ("ERROR: prenuc returned non-zero\n");
    }


    #-- Run mummer | mgaps and assert return value is zero
    print (STDERR "2,3: RUNNING mummer AND CREATING CLUSTERS\n");
    print("$algo_path $algo $mdir -l $size -n $pfx.ntref $qry_file |\n",
          "| $mgaps_path -l $clus -s $gap -d $ddiff -f $dfrac > $pfx.mgaps\n");
    open(ALGO_PIPE, "$algo_path $algo $mdir -l $size -n $pfx.ntref $qry_file |")
	or $tigr->bail ("ERROR: could not open $algo_path output pipe $!");
    open(CLUS_PIPE, "| $mgaps_path -l $clus -s $gap -d $ddiff -f $dfrac > $pfx.mgaps")
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


    #-- Run postnuc and assert return value is zero
    print (STDERR "4: FINISHING DATA\n");
    if ( $banded )
      {
        $err[0] = $tigr->runCommand
          ("$postnuc_path $psw -b $blen -B $ddiff $ref_file $qry_file $pfx < $pfx.mgaps");
      }
    else
      {
        $err[0] = $tigr->runCommand
          ("$postnuc_path $psw -b $blen $ref_file $qry_file $pfx < $pfx.mgaps");
      }

    if ( $err[0] != 0 ) {
	$tigr->bail ("ERROR: postnuc returned non-zero\n");
    }

    #-- If the -o flag was set, run show-coords using NUCmer1.1 settings
    if ( $generate_coords ) {
	print (STDERR "5: GENERATING COORDS FILE\n");
	$err[0] = $tigr->runCommand
	    ("$showcoords_path -r $pfx.delta > $pfx.coords");
	
	if ( $err[0] != 0 ) {
	    $tigr->bail ("ERROR: show-coords returned non-zero\n");
	}
    }

    #-- Remove the temporary output
#    $err[0] = unlink ("$pfx.ntref", "$pfx.mgaps");

    if ( $err[0] != 2 ) {
	$tigr->logError ("WARNING: there was a problem deleting".
			 " the temporary output files", 1);
    }

    #-- Return success
    return (0);
}

exit ( main ( ) );

#-- END OF SCRIPT
