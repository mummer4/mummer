#!@PERL@

################################################################################
#   Programmer: Adam M Phillippy, The Institute for Genomic Research
#         File: mummerplot
#         Date: 01 / 08 / 03
#               01 / 06 / 05 rewritten (v3.0)
#  
#        Usage:
#    mummerplot  [options]  <match file>
# 
#                Try 'mummerplot -h' for more information.
# 
#      Purpose: To generate a gnuplot plot for the display of mummer, nucmer,
#               promer, and show-tiling alignments.
# 
################################################################################

use lib "@LIB_DIR@";
use Foundation;
use strict;
use IO::Socket;

my $BIN_DIR     = "@BIN_DIR@";
my $LIB_DIR     = "@LIB_DIR@";
my $GNUPLOT_EXE = "@GNUPLOT_EXE@";


#================================================================= Globals ====#
#-- terminal types
my $X11    = "x11";
my $PS     = "postscript";
my $PNG    = "png";

#-- terminal sizes
my $SMALL  = "small";
my $MEDIUM = "medium";
my $LARGE  = "large";

my %TERMSIZE =
    (
     $X11 => { $SMALL => 500, $MEDIUM => 700,  $LARGE => 900  }, # screen pix
     $PS  => { $SMALL => 1,   $MEDIUM => 2,    $LARGE => 3    }, # pages
     $PNG => { $SMALL => 800, $MEDIUM => 1024, $LARGE => 1400 }  # image pix
     );

#-- terminal format
my $FFACE    = "Courier";
my $FSIZE    = "8";
my $TFORMAT  = "%.0f";
my $MFORMAT  = "[%.0f, %.0f]";

#-- output suffixes
my $FILTER  = "filter";
my $FWDPLOT = "fplot";
my $REVPLOT = "rplot";
my $HLTPLOT = "hplot";
my $GNUPLOT = "gnuplot";

my %SUFFIX =
    (
     $FILTER  => ".filter",
     $FWDPLOT => ".fplot",
     $REVPLOT => ".rplot",
     $HLTPLOT => ".hplot",
     $GNUPLOT => ".gp",
     $PS      => ".ps",
     $PNG     => ".png"
     );


#================================================================= Options ====#
my $OPT_breaklen;                  # -b option
my $OPT_color;                     # --[no]color option
my $OPT_coverage;                  # --[no]coverage option
my $OPT_filter;                    # -f option
my $OPT_layout;                    # -l option
my $OPT_prefix    = "out";         # -p option
my $OPT_rv;                        # --rv option
my $OPT_terminal  = $X11;          # -t option
my $OPT_IdR;                       # -r option
my $OPT_IdQ;                       # -q option
my $OPT_IDRfile;                   # -R option
my $OPT_IDQfile;                   # -Q option
my $OPT_rport;                     # -rport option
my $OPT_qport;                     # -qport option
my $OPT_size      = $SMALL;        # -small, -medium, -large
my $OPT_SNP;                       # -S option
my $OPT_xrange;                    # -x option
my $OPT_yrange;                    # -y option
my $OPT_title;                     # -title option

my $OPT_Mfile;                     # match file
my $OPT_Dfile;                     # delta filter file
my $OPT_Ffile;                     # .fplot output
my $OPT_Rfile;                     # .rplot output
my $OPT_Hfile;                     # .hplot output
my $OPT_Gfile;                     # .gp output
my $OPT_Pfile;                     # .ps .png output

my $OPT_gpstatus;                  # gnuplot status

my $OPT_ONLY_USE_FATTEST;          # Only use fattest alignment for layout


#============================================================== Foundation ====#
my $VERSION = '3.5';

my $USAGE = qq~
  USAGE: mummerplot  [options]  <match file>
    ~;

my $HELP = qq~
  USAGE: mummerplot  [options]  <match file>

  DESCRIPTION:
    mummerplot generates plots of alignment data produced by mummer, nucmer,
    promer or show-tiling by using the GNU gnuplot utility. After generating
    the appropriate scripts and datafiles, mummerplot will attempt to run
    gnuplot to generate the plot. If this attempt fails, a warning will be
    output and the resulting .gp and .[frh]plot files will remain so that the
    user may run gnuplot independently. If the attempt succeeds, either an x11
    window will be spawned or an additional output file will be generated
    (.ps or .png depending on the selected terminal). Feel free to edit the
    resulting gnuplot script (.gp) and rerun gnuplot to change line thinkness,
    labels, colors, plot size etc.

  MANDATORY:
    match file      Set the alignment input to 'match file'
                    Valid inputs are from mummer, nucmer, promer and
                    show-tiling (.out, .cluster, .delta and .tiling)

  OPTIONS:
    -b|breaklen     Highlight alignments with breakpoints further than
                    breaklen nucleotides from the nearest sequence end
    --[no]color     Color plot lines with a percent similarity gradient or
                    turn off all plot color (default color by match dir)
                    If the plot is very sparse, edit the .gp script to plot
                    with 'linespoints' instead of 'lines'
    -c
    --[no]coverage  Generate a reference coverage plot (default for .tiling)
    --depend        Print the dependency information and exit
    -f
    --filter        Only display .delta alignments which represent the "best"
                    hit to any particular spot on either sequence, i.e. a
                    one-to-one mapping of reference and query subsequences
    -h
    --help          Display help information and exit
    -l
    --layout        Layout a .delta multiplot in an intelligible fashion,
                    this option requires the -R -Q options
    --fat           Layout sequences using fattest alignment only
    -p|prefix       Set the prefix of the output files (default '$OPT_prefix')
    -rv             Reverse video for x11 plots
    -r|IdR          Plot a particular reference sequence ID on the X-axis
    -q|IdQ          Plot a particular query sequence ID on the Y-axis
    -R|Rfile        Plot an ordered set of reference sequences from Rfile
    -Q|Qfile        Plot an ordered set of query sequences from Qfile
                    Rfile/Qfile Can either be the original DNA multi-FastA
                    files or lists of sequence IDs, lens and dirs [ /+/-]
    -r|rport        Specify the port to send reference ID and position on
                    mouse double click in X11 plot window
    -q|qport        Specify the port to send query IDs and position on mouse
                    double click in X11 plot window
    -s|size         Set the output size to small, medium or large
                    --small --medium --large (default '$OPT_size')
    -S
    --SNP           Highlight SNP locations in each alignment
    -t|terminal     Set the output terminal to x11, postscript or png
                    --x11 --postscript --png (default '$OPT_terminal')
    -t|title        Specify the gnuplot plot title (default none)
    -x|xrange       Set the xrange for the plot '[min:max]'
    -y|yrange       Set the yrange for the plot '[min:max]'
    -V
    --version       Display the version information and exit
    ~;

my @DEPEND =
    (
     "$LIB_DIR/Foundation.pm",
     "$BIN_DIR/delta-filter",
     "$BIN_DIR/show-coords",
     "$BIN_DIR/show-snps",
     "gnuplot"
     );

my $tigr = new TIGR::Foundation
    or die "ERROR: TIGR::Foundation could not be initialized\n";

$tigr -> setVersionInfo ($VERSION);
$tigr -> setUsageInfo ($USAGE);
$tigr -> setHelpInfo ($HELP);
$tigr -> addDependInfo (@DEPEND);


#=========================================================== Function Decs ====#
sub GetParseFunc( );

sub ParseIDs($$);

sub ParseDelta($);
sub ParseCluster($);
sub ParseMummer($);
sub ParseTiling($);

sub LayoutIDs($$);
sub SpanXwY ($$$$$);

sub PlotData($$$);
sub WriteGP($$);
sub RunGP( );
sub ListenGP($$);

sub ParseOptions( );


#=========================================================== Function Defs ====#
MAIN:
{
    my @aligns;                # (sR eR sQ eQ sim lenR lenQ idR idQ)
    my %refs;                  # (id => (off, len, [1/-1]))
    my %qrys;                  # (id => (off, len, [1/-1]))

    #-- Get the command line options (sets OPT_ global vars)
    ParseOptions( );


    #-- Get the alignment type
    my $parsefunc = GetParseFunc( );

    if ( $parsefunc != \&ParseDelta &&
         ($OPT_filter || $OPT_layout || $OPT_SNP) ) {
        print STDERR "WARNING: -f -l -S only work with delta input\n";
        undef $OPT_filter;
        undef $OPT_layout;
        undef $OPT_SNP;
    }

    #-- Parse the reference and query IDs
    if    ( defined $OPT_IdR ) { $refs{$OPT_IdR} = [ 0, 0, 1 ]; }
    elsif ( defined $OPT_IDRfile ) {
        ParseIDs ($OPT_IDRfile, \%refs);
    }

    if    ( defined $OPT_IdQ ) { $qrys{$OPT_IdQ} = [ 0, 0, 1 ]; }
    elsif ( defined $OPT_IDQfile ) {
        ParseIDs ($OPT_IDQfile, \%qrys);
    }


    #-- Filter the alignments
    if ( $OPT_filter || $OPT_layout ) {
        print STDERR "Writing filtered delta file $OPT_Dfile\n";
        system ("$BIN_DIR/delta-filter -r -q $OPT_Mfile > $OPT_Dfile")
            and die "ERROR: Could not run delta-filter, $!\n";
        if ( $OPT_filter ) { $OPT_Mfile = $OPT_Dfile; }
    }


    #-- Parse the alignment data
    $parsefunc->(\@aligns);

    
    #-- Layout the alignment data if requested
    if ( $OPT_layout ) {
        if ( scalar (keys %refs) || scalar (keys %qrys) ) {
            LayoutIDs (\%refs, \%qrys);
        }
        else {
            print STDERR "WARNING: --layout option only works with -R or -Q\n";
            undef $OPT_layout;
        }
    }


    #-- Plot the alignment data
    PlotData (\@aligns, \%refs, \%qrys);


    #-- Write the gnuplot script
    WriteGP (\%refs, \%qrys);


    #-- Run gnuplot script and fork a clipboard listener
    unless ( $OPT_gpstatus == -1 ) {

        my $child = 1;
        if ( $OPT_gpstatus == 0 && $OPT_terminal eq $X11 ) {
            print STDERR "Forking mouse listener\n";
            $child = fork;
        }

        #-- parent runs gnuplot
        if ( $child ) {
            RunGP( );
            kill 1, $child;
        }
        #-- child listens to clipboard
        elsif ( defined $child ) {
            ListenGP(\%refs, \%qrys);
        }
        else {
            print STDERR "WARNING: Could not fork mouse listener\n";
        }
    }

    exit (0);
}


#------------------------------------------------------------ GetParseFunc ----#
sub GetParseFunc ( )
{
    my $fref;

    open (MFILE, "<$OPT_Mfile")
        or die "ERROR: Could not open $OPT_Mfile, $!\n";

    $_ = <MFILE>;
    if ( !defined ) { die "ERROR: Could not read $OPT_Mfile, File is empty\n" }

  SWITCH: {
      #-- tiling
      if ( /^>\S+ \d+ bases/ ) {
          $fref = \&ParseTiling;
          last SWITCH;
      }

      #-- mummer
      if ( /^> \S+/ ) {
          $fref = \&ParseMummer;
          last SWITCH;
      }

      #-- nucmer/promer
      if ( /^(\S+) (\S+)/ ) {
          if ( ! defined $OPT_IDRfile ) {
              $OPT_IDRfile = $1;
          }
          if ( ! defined $OPT_IDQfile ) {
              $OPT_IDQfile = $2;
          }

          $_ = <MFILE>;
          if ( (defined)  &&  (/^NUCMER$/  ||  /^PROMER$/) ) {
              $_ = <MFILE>;   # sequence header
              $_ = <MFILE>;   # alignment header
              if ( !defined ) {
                  $fref = \&ParseDelta;
                  last SWITCH;
              }
              elsif ( /^\d+ \d+ \d+ \d+ \d+ \d+ \d+$/ ) {
                  $fref = \&ParseDelta;
                  last SWITCH;
              }
              elsif ( /^[ \-][1-3] [ \-][1-3]$/ ) {
                  $fref = \&ParseCluster;
                  last SWITCH;
              }
          }
      }
      
      #-- default
      die "ERROR: Could not read $OPT_Mfile, Unrecognized file type\n";
  }

    close (MFILE)
        or print STDERR "WARNING: Trouble closing $OPT_Mfile, $!\n";

    return $fref;
}


#---------------------------------------------------------------- ParseIDs ----#
sub ParseIDs ($$)
{
    my $file = shift;
    my $href = shift;

    open (IDFILE, "<$file")
        or print STDERR "WARNING: Could not open $file, $!\n";

    my $dir;
    my $aref;
    my $isfasta;
    my $offset = 0;
    while ( <IDFILE> ) {
        #-- Ignore blank lines
        if ( /^\s*$/ ) { next; }

        #-- FastA header
        if ( /^>(\S+)/ ) {
            if ( exists $href->{$1} ) {
                print STDERR "WARNING: Duplicate sequence '$1' ignored\n";
                undef $aref;
                next;
            }

            if ( !$isfasta ) { $isfasta = 1; }
            if ( defined $aref ) { $offset += $aref->[1] - 1; }

            $aref = [ $offset, 0, 1 ];
            $href->{$1} = $aref;
            next;
        }
        
        #-- FastA sequence
        if ( $isfasta  &&  /^\S+$/ ) {
            if ( defined $aref ) { $aref->[1] += (length) - 1; }
            next;
        }

        #-- ID len dir
        if ( !$isfasta  &&  /^(\S+)\s+(\d+)\s+([+-]?)$/ ) {
            if ( exists $href->{$1} ) {
                print STDERR "WARNING: Duplicate sequence '$1' ignored\n";
                undef $aref;
                next;
            }

            $dir = (defined $3 && $3 eq "-") ? -1 : 1;
            $aref = [ $offset, $2, $dir ];
            $offset += $2 - 1;
            $href->{$1} = $aref;
            next;
        }

        #-- default
        print STDERR "WARNING: Could not parse $file\n$_";
        undef %$href;
        last;
    }

    close (IDFILE)
        or print STDERR "WARNING: Trouble closing $file, $!\n";
}


#-------------------------------------------------------------- ParseDelta ----#
sub ParseDelta ($)
{
    my $aref = shift;

    print STDERR "Reading delta file $OPT_Mfile\n";

    open (MFILE, "<$OPT_Mfile")
        or die "ERROR: Could not open $OPT_Mfile, $!\n";

    my @align;
    my $ispromer;
    my ($sim, $tot);
    my ($lenR, $lenQ, $idR, $idQ);

    $_ = <MFILE>;
    $_ = <MFILE>;
    $ispromer = /^PROMER/;

    while ( <MFILE> ) {
        #-- delta int
        if ( /^([-]?\d+)$/ ) {
            if ( $1 < 0 ) {
                $tot ++;
            }
            elsif ( $1 == 0 ) {
                $align[4] = ($tot - $sim) / $tot * 100.0;
                push @$aref, [ @align ];
                $tot = 0;
            }
            next;
        }

        #-- alignment header
        if ( /^(\d+) (\d+) (\d+) (\d+) \d+ (\d+) \d+$/ ) {
            if ( $tot == 0 ) {
                @align = ($1, $2, $3, $4, 0, $lenR, $lenQ, $idR, $idQ);
                $tot = abs($1 - $2) + 1;
                $sim = $5;
                if ( $ispromer ) { $tot /= 3.0; }
                next;
            }
            #-- drop to default
        }

        #-- sequence header
        if ( /^>(\S+) (\S+) (\d+) (\d+)$/ ) {
            ($idR, $idQ, $lenR, $lenQ) = ($1, $2, $3, $4);
            $tot = 0;
            next;
        }

        #-- default
        die "ERROR: Could not parse $OPT_Mfile\n$_";
    }

    close (MFILE)
        or print STDERR "WARNING: Trouble closing $OPT_Mfile, $!\n";
}


#------------------------------------------------------------ ParseCluster ----#
sub ParseCluster ($)
{
    my $aref = shift;

    print STDERR "Reading cluster file $OPT_Mfile\n";

    open (MFILE, "<$OPT_Mfile")
        or die "ERROR: Could not open $OPT_Mfile, $!\n";

    my @align;
    my ($dR, $dQ, $len);
    my ($lenR, $lenQ, $idR, $idQ);

    $_ = <MFILE>;
    $_ = <MFILE>;

    while ( <MFILE> ) {
        #-- match
        if ( /^\s+(\d+)\s+(\d+)\s+(\d+)\s+\S+\s+\S+$/ ) {
            @align = ($1, $1, $2, $2, 100, $lenR, $lenQ, $idR, $idQ);
            $len = $3 - 1;
            $align[1] += $dR == 1 ? $len : -$len;
            $align[3] += $dQ == 1 ? $len : -$len;
            push @$aref, [ @align ];
            next;
        }

        #-- cluster header
        if ( /^[ \-][1-3] [ \-][1-3]$/ ) {
            $dR = /^-/ ? -1 : 1;
            $dQ = /-[1-3]$/ ? -1 : 1;
            next;
        }

        #-- sequence header
        if ( /^>(\S+) (\S+) (\d+) (\d+)$/ ) {
            ($idR, $idQ, $lenR, $lenQ) = ($1, $2, $3, $4);
            next;
        }

        #-- default
        die "ERROR: Could not parse $OPT_Mfile\n$_";
    }

    close (MFILE)
        or print STDERR "WARNING: Trouble closing $OPT_Mfile, $!\n";
}


#------------------------------------------------------------- ParseMummer ----#
sub ParseMummer ($)
{
    my $aref = shift;

    print STDERR "Reading mummer file $OPT_Mfile (use mummer -c)\n";

    open (MFILE, "<$OPT_Mfile")
        or die "ERROR: Could not open $OPT_Mfile, $!\n";

    my @align;
    my ($dQ, $len);
    my ($lenQ, $idQ);

    while ( <MFILE> ) {
        #-- 3 column match
        if ( /^\s+(\d+)\s+(\d+)\s+(\d+)$/ ) {
            @align = ($1, $1, $2, $2, 100, 0, $lenQ, "REF", $idQ);
            $len = $3 - 1;
            $align[1] += $len;
            $align[3] += $dQ == 1 ? $len : -$len;
            push @$aref, [ @align ];
            next;
        }

        #-- 4 column match
        if ( /^\s+(\S+)\s+(\d+)\s+(\d+)\s+(\d+)$/ ) {
            @align = ($2, $2, $3, $3, 100, 0, $lenQ, $1, $idQ);
            $len = $4 - 1;
            $align[1] += $len;
            $align[3] += $dQ == 1 ? $len : -$len;
            push @$aref, [ @align ];
            next;
        }

        #-- sequence header
        if ( /^> (\S+)/ ) {
            $idQ = $1;
            $dQ = /^> \S+ Reverse/ ? -1 : 1;
            $lenQ = /Len = (\d+)/ ? $1 : 0;
            next;
        }

        #-- default
        die "ERROR: Could not parse $OPT_Mfile\n$_";
    }

    close (MFILE)
        or print STDERR "WARNING: Trouble closing $OPT_Mfile, $!\n";
}


#------------------------------------------------------------- ParseTiling ----#
sub ParseTiling ($)
{
    my $aref = shift;

    print STDERR "Reading tiling file $OPT_Mfile\n";

    open (MFILE, "<$OPT_Mfile")
        or die "ERROR: Could not open $OPT_Mfile, $!\n";

    my @align;
    my ($dR, $dQ, $len);
    my ($lenR, $lenQ, $idR, $idQ);

    while ( <MFILE> ) {
        #-- tile
        if ( /^(\S+)\s+\S+\s+\S+\s+(\d+)\s+\S+\s+(\S+)\s+([+-])\s+(\S+)$/ ) {
            @align = ($1, $1, 1, 1, $3, $lenR, $2, $idR, $5);
            $len = $2 - 1;
            $align[1] += $len;
            $align[($4 eq "-" ? 2 : 3)] += $len;
            push @$aref, [ @align ];
            next;
        }

        #-- sequence header
        if ( /^>(\S+) (\d+) bases$/ ) {
            ($idR, $lenR) = ($1, $2);
            next;
        }

        #-- default
        die "ERROR: Could not parse $OPT_Mfile\n$_";
    }

    close (MFILE)
        or print STDERR "WARNING: Trouble closing $OPT_Mfile, $!\n";

    if ( ! defined $OPT_coverage ) { $OPT_coverage = 1; }
}


#--------------------------------------------------------------- LayoutIDs ----#
# For each reference and query sequence, find the set of alignments that
# produce the heaviest (both in non-redundant coverage and percent
# identity) alignment subset of each sequence using a modified version
# of the longest increasing subset algorithm. Let R be the union of all
# reference LIS subsets, and Q be the union of all query LIS
# subsets. Let S be the intersection of R and Q. Using this LIS subset,
# recursively span reference and query sequences by their smaller
# counterparts until all spanning sequences have been placed. The goal
# is to cluster all the "major" alignment information along the main
# diagonal for easy viewing and interpretation.
sub LayoutIDs ($$)
{
    my $rref = shift;
    my $qref = shift;

    my %rc;          # chains of qry seqs needed to span each ref
    my %qc;          # chains of ref seqs needed to span each qry
    #  {idR} -> [ placed, len, {idQ} -> [ \slope, \loR, \hiR, \loQ, \hiQ ] ]
    #  {idQ} -> [ placed, len, {idR} -> [ \slope, \loQ, \hiQ, \loR, \hiR ] ]

    my @rl;          # oo of ref seqs
    my @ql;          # oo of qry seqs
    #  [ [idR, slope] ]
    #  [ [idQ, slope] ]

    #-- get the filtered alignments
    open (BTAB, "$BIN_DIR/show-coords -B $OPT_Dfile |")
        or die "ERROR: Could not open show-coords pipe, $!\n";

    my @align;
    my ($sR, $eR, $sQ, $eQ, $lenR, $lenQ, $idR, $idQ);
    my ($loR, $hiR, $loQ, $hiQ);
    my ($dR, $dQ, $slope);
    while ( <BTAB> ) {
        chomp;
        @align = split "\t";
        if ( scalar @align != 21 ) {
            die "ERROR: Could not read show-coords pipe, invalid btab format\n";
        }

        $sR   = $align[8];   $eR   = $align[9];
        $sQ   = $align[6];   $eQ   = $align[7];
        $lenR = $align[18];  $lenQ = $align[2];
        $idR  = $align[5];   $idQ  = $align[0];

        #-- skip it if not on include list
        if ( !exists $rref->{$idR} || !exists $qref->{$idQ} ) { next; }

        #-- get orientation of both alignments and alignment slope
        $dR = $sR < $eR ? 1 : -1;
        $dQ = $sQ < $eQ ? 1 : -1;
        $slope = $dR == $dQ ? 1 : -1;

        #-- get lo's and hi's
        $loR = $dR == 1 ? $sR : $eR;
        $hiR = $dR == 1 ? $eR : $sR;

        $loQ = $dQ == 1 ? $sQ : $eQ;
        $hiQ = $dQ == 1 ? $eQ : $sQ;

        if ($OPT_ONLY_USE_FATTEST)
        {
          #-- Check to see if there is another better alignment
          if (exists $qc{$idQ})
          {
            my ($oldR) = keys %{$qc{$idQ}[2]};
            my $val = $qc{$idQ}[2]{$oldR};

            if (${$val->[4]} - ${$val->[3]} > $hiR - $loR)
            {
              #-- Old alignment is better, skip this one
              next;
            }
            else
            {
              #-- This alignment is better, prune old alignment
              delete $rc{$oldR}[2]{$idQ};
              delete $qc{$idQ};
            }
          }
        }

        #-- initialize
        if ( !exists $rc{$idR} ) { $rc{$idR} = [ 0, $lenR, { } ]; }
        if ( !exists $qc{$idQ} ) { $qc{$idQ} = [ 0, $lenQ, { } ]; }

        #-- if no alignments for these two exist OR
        #-- this alignment is bigger than the current
        if ( !exists $rc{$idR}[2]{$idQ} || !exists $qc{$idQ}[2]{$idR} ||
             $hiR - $loR >
             ${$rc{$idR}[2]{$idQ}[2]} - ${$rc{$idR}[2]{$idQ}[1]} ) {

            #-- rc and qc reference these anonymous values
            my $aref = [ $slope, $loR, $hiR, $loQ, $hiQ ];

            #-- rc is ordered [ slope, loR, hiR, loQ, hiQ ]
            #-- qc is ordered [ slope, loQ, hiQ, loR, hiR ]
            $rc{$idR}[2]{$idQ}[0] = $qc{$idQ}[2]{$idR}[0] = \$aref->[0];
            $rc{$idR}[2]{$idQ}[1] = $qc{$idQ}[2]{$idR}[3] = \$aref->[1];
            $rc{$idR}[2]{$idQ}[2] = $qc{$idQ}[2]{$idR}[4] = \$aref->[2];
            $rc{$idR}[2]{$idQ}[3] = $qc{$idQ}[2]{$idR}[1] = \$aref->[3];
            $rc{$idR}[2]{$idQ}[4] = $qc{$idQ}[2]{$idR}[2] = \$aref->[4];
        }
    }

    close (BTAB)
        or print STDERR "WARNING: Trouble closing show-coords pipe, $!\n";

    #-- recursively span sequences to generate the layout
    foreach $idR ( sort { $rc{$b}[1] <=> $rc{$a}[1] } keys %rc ) {
        SpanXwY ($idR, \%rc, \@rl, \%qc, \@ql);
    }

    #-- undefine the current offsets
    foreach $idR ( keys %{$rref} ) { undef $rref->{$idR}[0]; }
    foreach $idQ ( keys %{$qref} ) { undef $qref->{$idQ}[0]; }

    #-- redefine the offsets according to the new layout
    my $roff = 0;
    foreach my $r ( @rl ) {
        $idR = $r->[0];
        $rref->{$idR}[0] = $roff;
        $rref->{$idR}[2] = $r->[1];
        $roff += $rref->{$idR}[1] - 1;
    }
    #-- append the guys left out of the layout
    foreach $idR ( keys %{$rref} ) {
        if ( !defined $rref->{$idR}[0] ) {
            $rref->{$idR}[0] = $roff;
            $roff += $rref->{$idR}[1] - 1;
        }
    }

    #-- redefine the offsets according to the new layout
    my $qoff = 0;
    foreach my $q ( @ql ) {
        $idQ = $q->[0];
        $qref->{$idQ}[0] = $qoff;
        $qref->{$idQ}[2] = $q->[1];
        $qoff += $qref->{$idQ}[1] - 1;
    }
    #-- append the guys left out of the layout
    foreach $idQ ( keys %{$qref} ) {
        if ( !defined $qref->{$idQ}[0] ) {
            $qref->{$idQ}[0] = $qoff;
            $qoff += $qref->{$idQ}[1] - 1;
        }
    }
}


#----------------------------------------------------------------- SpanXwY ----#
sub SpanXwY ($$$$$) {
    my $x   = shift;   # idX
    my $xcr = shift;   # xc ref
    my $xlr = shift;   # xl ref
    my $ycr = shift;   # yc ref
    my $ylr = shift;   # yl ref

    my @post;
    foreach my $y ( sort { ${$xcr->{$x}[2]{$a}[1]} <=> ${$xcr->{$x}[2]{$b}[1]} }
                    keys %{$xcr->{$x}[2]} ) {

        #-- skip if already placed (RECURSION BASE)
        if ( $ycr->{$y}[0] ) { next; }
        else { $ycr->{$y}[0] = 1; }

        #-- get len and slope info for y
        my $len = $ycr->{$y}[1];
        my $slope = ${$xcr->{$x}[2]{$y}[0]};

        #-- if we need to flip, reverse complement all y records
        if ( $slope == -1 ) {
            foreach my $xx ( keys %{$ycr->{$y}[2]} ) {
                ${$ycr->{$y}[2]{$xx}[0]} *= -1;

                my $loy = ${$ycr->{$y}[2]{$xx}[1]};
                my $hiy = ${$ycr->{$y}[2]{$xx}[2]};
                ${$ycr->{$y}[2]{$xx}[1]} = $len - $hiy + 1;
                ${$ycr->{$y}[2]{$xx}[2]} = $len - $loy + 1;
            }
        }

        #-- place y
        push @{$ylr}, [ $y, $slope ];

        #-- RECURSE if y > x, else save for later
        if ( $len > $xcr->{$x}[1] ) { SpanXwY ($y, $ycr, $ylr, $xcr, $xlr); }
        else { push @post, $y; }
    }

    #-- RECURSE for all y < x
    foreach my $y ( @post ) { SpanXwY ($y, $ycr, $ylr, $xcr, $xlr); }
}


#---------------------------------------------------------------- PlotData ----#
sub PlotData ($$$)
{
    my $aref = shift;
    my $rref = shift;
    my $qref = shift;

    print STDERR "Writing plot files $OPT_Ffile, $OPT_Rfile",
    (defined $OPT_Hfile ? ", $OPT_Hfile\n" : "\n");

    open (FFILE, ">$OPT_Ffile")
        or die "ERROR: Could not open $OPT_Ffile, $!\n";
    print FFILE "#-- forward hits sorted by %sim\n0 0 0\n0 0 0\n\n\n";

    open (RFILE, ">$OPT_Rfile")
        or die "ERROR: Could not open $OPT_Rfile, $!\n";
    print RFILE "#-- reverse hits sorted by %sim\n0 0 0\n0 0 0\n\n\n";

    if ( defined $OPT_Hfile ) {
        open (HFILE, ">$OPT_Hfile")
            or die "ERROR: Could not open $OPT_Hfile, $!\n";
        print HFILE "#-- highlighted hits sorted by %sim\n0 0 0\n0 0 0\n\n\n";
    }

    my $fh;
    my $align;
    my $isplotted;
    my $ismultiref;
    my $ismultiqry;
    my ($plenR, $plenQ, $pidR, $pidQ);

    #-- for each alignment sorted by ascending identity
    foreach $align ( sort { $a->[4] <=> $b->[4] } @$aref ) {

        my ($sR, $eR, $sQ, $eQ, $sim, $lenR, $lenQ, $idR, $idQ) = @$align;

        if ( ! defined $pidR ) {
            ($plenR, $plenQ, $pidR, $pidQ) = ($lenR, $lenQ, $idR, $idQ);
        }

        #-- set the sequence offset, length, direction, etc...
        my ($refoff, $reflen, $refdir);
        my ($qryoff, $qrylen, $qrydir);

        if ( %$rref ) {
            #-- skip reference sequence or set atts from hash
            if ( !exists ($rref->{$idR}) ) { next; }
            else { ($refoff, $reflen, $refdir) = @{$rref->{$idR}}; }
        }
        else {
            #-- no reference hash, so default atts
            ($refoff, $reflen, $refdir) = (0, $lenR, 1);
        }

        if ( %$qref ) {
            #-- skip query sequence or set atts from hash
            if ( !exists ($qref->{$idQ}) ) { next; }
            else { ($qryoff, $qrylen, $qrydir) = @{$qref->{$idQ}}; }
        }
        else {
            #-- no query hash, so default atts
            ($qryoff, $qrylen, $qrydir) = (0, $lenQ, 1);
        }

        #-- get the orientation right
        if ( $refdir == -1 ) {
            $sR = $reflen - $sR + 1;
            $eR = $reflen - $eR + 1;
        }
        if ( $qrydir == -1 ) {
            $sQ = $qrylen - $sQ + 1;
            $eQ = $qrylen - $eQ + 1;
        }

        #-- forward file, reverse file, highlight file
        my @fha;

        if ( defined $OPT_breaklen &&
             ( ($sR - 1 > $OPT_breaklen &&
                $sQ - 1 > $OPT_breaklen &&
                $reflen - $sR > $OPT_breaklen &&
                $qrylen - $sQ > $OPT_breaklen)
               ||
               ($eR - 1 > $OPT_breaklen &&
                $eQ - 1 > $OPT_breaklen &&
                $reflen - $eR > $OPT_breaklen &&
                $qrylen - $eQ > $OPT_breaklen) ) ) {
            push @fha, \*HFILE;
        }

        push @fha, (($sR < $eR) == ($sQ < $eQ) ? \*FFILE : \*RFILE);

        #-- plot it
        $sR += $refoff; $eR += $refoff;
        $sQ += $qryoff; $eQ += $qryoff;

        if ( $OPT_coverage ) {
            foreach $fh ( @fha ) {
                print $fh
                    "$sR 10 $sim\n", "$eR 10 $sim\n\n\n",
                    "$sR $sim 0\n", "$eR $sim 0\n\n\n";
            }
        }
        else {
            foreach $fh ( @fha ) {
                print $fh "$sR $sQ $sim\n", "$eR $eQ $sim\n\n\n";
            }
        }            

        #-- set some flags
        if ( !$ismultiref && $idR ne $pidR ) { $ismultiref = 1; }
        if ( !$ismultiqry && $idQ ne $pidQ ) { $ismultiqry = 1; }
        if ( !$isplotted ) { $isplotted = 1; }
    }


    #-- highlight the SNPs
    if ( defined $OPT_SNP ) {

        print STDERR "Determining SNPs from sequence and alignment data\n";

        open (SNPS, "$BIN_DIR/show-snps -H -T -l $OPT_Mfile |")
            or die "ERROR: Could not open show-snps pipe, $!\n";

        my @snps;
        my ($pR, $pQ, $lenR, $lenQ, $idR, $idQ);
        while ( <SNPS> ) {
            chomp;
            @snps = split "\t";
            if ( scalar @snps != 14 ) {
                die "ERROR: Could not read show-snps pipe, invalid format\n";
            }

            $pR   = $snps[0];   $pQ   = $snps[3];
            $lenR = $snps[8];   $lenQ = $snps[9];
            $idR  = $snps[12];  $idQ  = $snps[13];

            #-- set the sequence offset, length, direction, etc...
            my ($refoff, $reflen, $refdir);
            my ($qryoff, $qrylen, $qrydir);
            
            if ( %$rref ) {
                #-- skip reference sequence or set atts from hash
                if ( !exists ($rref->{$idR}) ) { next; }
                else { ($refoff, $reflen, $refdir) = @{$rref->{$idR}}; }
            }
            else {
                #-- no reference hash, so default atts
                ($refoff, $reflen, $refdir) = (0, $lenR, 1);
            }
            
            if ( %$qref ) {
                #-- skip query sequence or set atts from hash
                if ( !exists ($qref->{$idQ}) ) { next; }
                else { ($qryoff, $qrylen, $qrydir) = @{$qref->{$idQ}}; }
            }
            else {
                #-- no query hash, so default atts
                ($qryoff, $qrylen, $qrydir) = (0, $lenQ, 1);
            }

            #-- get the orientation right
            if ( $refdir == -1 ) { $pR = $reflen - $pR + 1; }
            if ( $qrydir == -1 ) { $pQ = $qrylen - $pQ + 1; }

            #-- plot it
            $pR += $refoff;
            $pQ += $qryoff;

            if ( $OPT_coverage ) {
                print HFILE "$pR 10 0\n", "$pR 10 0\n\n\n",
            }
            else {
                print HFILE "$pR $pQ 0\n", "$pR $pQ 0\n\n\n";
            }            
        }

        close (SNPS)
            or print STDERR "WARNING: Trouble closing show-snps pipe, $!\n";
    }


    close (FFILE)
        or print STDERR "WARNING: Trouble closing $OPT_Ffile, $!\n";

    close (RFILE)
        or print STDERR "WARNING: Trouble closing $OPT_Rfile, $!\n";

    if ( defined $OPT_Hfile ) {
        close (HFILE)
            or print STDERR "WARNING: Trouble closing $OPT_Hfile, $!\n";
    }


    if ( !%$rref ) {
        if ( $ismultiref ) {
            print STDERR
                "WARNING: Multiple ref sequences overlaid, try -R or -r\n";
        }
        elsif ( defined $pidR ) {
            $rref->{$pidR} = [ 0, $plenR, 1 ];
        }
    }

    if ( !%$qref ) {
        if ( $ismultiqry && !$OPT_coverage ) {
            print STDERR
                "WARNING: Multiple qry sequences overlaid, try -Q, -q or -c\n";
        }
        elsif ( defined $pidQ ) {
            $qref->{$pidQ} = [ 0, $plenQ, 1 ];
        }
    }

    if ( !$isplotted ) {
        die "ERROR: No alignment data to plot\n";
    }
}


#----------------------------------------------------------------- WriteGP ----#
sub WriteGP ($$)
{
    my $rref = shift;
    my $qref = shift;

    print STDERR "Writing gnuplot script $OPT_Gfile\n";

    open (GFILE, ">$OPT_Gfile")
        or die "ERROR: Could not open $OPT_Gfile, $!\n";

    my ($FWD, $REV, $HLT) = (1, 2, 3);
    my $SIZE = $TERMSIZE{$OPT_terminal}{$OPT_size};

    #-- terminal specific stuff
    my ($P_TERM, $P_SIZE, %P_PS, %P_LW);
    foreach ( $OPT_terminal ) {
        /^$X11/    and do {
            $P_TERM = $OPT_gpstatus == 0 ?
                "$X11 font \"$FFACE,$FSIZE\"" : "$X11";

            %P_PS = ( $FWD => 1.0, $REV => 1.0, $HLT => 1.0 );

            %P_LW = $OPT_coverage || $OPT_color ?
                ( $FWD => 3.0, $REV => 3.0, $HLT => 3.0 ) :
                ( $FWD => 2.0, $REV => 2.0, $HLT => 2.0 );

            $P_SIZE = $OPT_coverage ?
                "set size 1,1" :
                "set size 1,1";

            last;
        };

        /^$PS/     and do {
            $P_TERM = defined $OPT_color && $OPT_color == 0 ?
                "$PS monochrome" : "$PS color";
            $P_TERM .= $OPT_gpstatus == 0 ?
                " solid \"$FFACE\" $FSIZE" : " solid \"$FFACE\" $FSIZE";

            %P_PS = ( $FWD => 0.5, $REV => 0.5, $HLT => 0.5 );

            %P_LW = $OPT_coverage || $OPT_color ?
                ( $FWD => 4.0, $REV => 4.0, $HLT => 4.0 ) :
                ( $FWD => 2.0, $REV => 2.0, $HLT => 2.0 );

            $P_SIZE = $OPT_coverage ?
                "set size ".(1.0 * $SIZE).",".(0.5 * $SIZE) :
                "set size ".(1.0 * $SIZE).",".(1.0 * $SIZE);

            last;
        };

        /^$PNG/    and do {
            $P_TERM = $OPT_gpstatus == 0 ?
                "$PNG tiny size $SIZE,$SIZE" : "$PNG small";
            if ( defined $OPT_color && $OPT_color == 0 ) {
                $P_TERM .= " xffffff x000000 x000000";
                $P_TERM .= " x000000 x000000 x000000";
                $P_TERM .= " x000000 x000000 x000000";
            }
            
            %P_PS = ( $FWD => 1.0, $REV => 1.0, $HLT => 1.0 );

            %P_LW = $OPT_coverage || $OPT_color ?
                ( $FWD => 3.0, $REV => 3.0, $HLT => 3.0 ) :
                ( $FWD => 3.0, $REV => 3.0, $HLT => 3.0 );

            $P_SIZE = $OPT_coverage ?
                "set size 1,.375" :
                "set size 1,1";

            last;
        };

        die "ERROR: Don't know how to initialize terminal, $OPT_terminal\n";
    }

    #-- plot commands
    my ($P_WITH, $P_FORMAT, $P_LS, $P_KEY, %P_PT, %P_LT);

    %P_PT = ( $FWD => 6, $REV => 6, $HLT => 6 );
    %P_LT = defined $OPT_Hfile ?
        ( $FWD => 2, $REV => 2, $HLT => 1 ) :
        ( $FWD => 1, $REV => 3, $HLT => 2 );

    $P_WITH = $OPT_coverage || $OPT_color ? "w l" : "w lp";

    $P_FORMAT = "set format \"$TFORMAT\"";
    if ( $OPT_gpstatus == 0 ) {
        $P_LS = "set style line";
        $P_KEY = "unset key";
        $P_FORMAT .= "\nset mouse format \"$TFORMAT\"";
        $P_FORMAT .= "\nset mouse mouseformat \"$MFORMAT\"";
        #$P_FORMAT .= "\nif(GPVAL_VERSION < 5) { set mouse clipboardformat \"$MFORMAT\" }";
    }
    else {
        $P_LS = "set linestyle";
        $P_KEY = "set nokey";
    }


    my @refk = keys (%$rref);
    my @qryk = keys (%$qref);
    my ($xrange, $yrange);
    my ($xlabel, $ylabel);
    my ($tic, $dir);
    my $border = 0;

    #-- terminal header and output
    print GFILE "set terminal $P_TERM\n";

    if ( defined $OPT_Pfile ) {
        print GFILE "set output \"$OPT_Pfile\"\n";
    }

    if ( defined $OPT_title ) {
        print GFILE "set title \"$OPT_title\"\n";
    }

    #-- set tics, determine labels, ranges (ref)
    if ( scalar (@refk) == 1 ) {
        $xlabel = $refk[0];
        $xrange = $rref->{$xlabel}[1];
    }
    else {
        $xrange = 0;
        print GFILE "set xtics rotate \( \\\n";
        foreach $xlabel ( sort { $rref->{$a}[0] <=> $rref->{$b}[0] } @refk ) {
            $xrange += $rref->{$xlabel}[1];
            $tic = $rref->{$xlabel}[0] + 1;
            $dir = ($rref->{$xlabel}[2] == 1) ? "" : "*";
            print GFILE " \"$dir$xlabel\" $tic.0, \\\n";
        }
        print GFILE " \"\" $xrange \\\n\)\n";
        $xlabel = "REF";
    }
    if ( $xrange == 0 ) { $xrange = "*"; }

    #-- set tics, determine labels, ranges (qry)
    if ( $OPT_coverage ) {
        $ylabel = "%SIM";
        $yrange = 110;
    }
    elsif ( scalar (@qryk) == 1 ) {
        $ylabel = $qryk[0];
        $yrange = $qref->{$ylabel}[1];
    }
    else {
        $yrange = 0;
        print GFILE "set ytics \( \\\n";
        foreach $ylabel ( sort { $qref->{$a}[0] <=> $qref->{$b}[0] } @qryk ) {
            $yrange += $qref->{$ylabel}[1];
            $tic = $qref->{$ylabel}[0] + 1;
            $dir = ($qref->{$ylabel}[2] == 1) ? "" : "*";
            print GFILE " \"$dir$ylabel\" $tic.0, \\\n";
        }
        print GFILE " \"\" $yrange \\\n\)\n";
        $ylabel = "QRY";
    }
    if ( $yrange == 0 ) { $yrange = "*"; }

    #-- determine borders
    if ( $xrange ne "*" && scalar (@refk) == 1 ) { $border |= 10; }
    if ( $yrange ne "*" && scalar (@qryk) == 1 ) { $border |= 5; }
    if ( $OPT_coverage ) { $border |= 5; }

    #-- grid, labels, border
    print GFILE
        "$P_SIZE\n",
        "set grid\n",
        "$P_KEY\n",
        "set border $border\n",
        "set tics scale 0\n",
        "set xlabel \"$xlabel\"\n",
        "set ylabel \"$ylabel\"\n",
        "$P_FORMAT\n";

    #-- ranges
    if ( defined $OPT_xrange ) { print GFILE "set xrange $OPT_xrange\n"; }
    else                       { print GFILE "set xrange [1:$xrange]\n"; }

    if ( defined $OPT_yrange ) { print GFILE "set yrange $OPT_yrange\n"; }
    else                       { print GFILE "set yrange [1:$yrange]\n"; }

    #-- if %sim plot
    if ( $OPT_color ) {
        print GFILE
            "set zrange [0:100]\n",
            "set colorbox default\n",
            "set cblabel \"%similarity\"\n",
            "set cbrange [0:100]\n",
            "set cbtics 20\n",
            "set pm3d map\n",
            "set palette model RGB defined ( \\\n",
            "  0 \"#000000\", \\\n",
            "  4 \"#DD00DD\", \\\n",
            "  6 \"#0000DD\", \\\n",
            "  7 \"#00DDDD\", \\\n",
            "  8 \"#00DD00\", \\\n",
            "  9 \"#DDDD00\", \\\n",
            " 10 \"#DD0000\"  \\\n)\n";
    }

    foreach my $s ( ($FWD, $REV, $HLT) ) {
        my $ss = "$P_LS $s ";
        $ss .= $OPT_color ? " palette" : " lt $P_LT{$s}";
        $ss .= " lw $P_LW{$s}";
        if ( ! $OPT_coverage || $s == $HLT ) {
            $ss .= " pt $P_PT{$s} ps $P_PS{$s}";
        }
        print GFILE "$ss\n";
    }

    #-- plot it
    print GFILE
        ($OPT_color ? "splot \\\n" : "plot \\\n");
    print GFILE
        " \"$OPT_Ffile\" title \"FWD\" $P_WITH ls $FWD, \\\n",
        " \"$OPT_Rfile\" title \"REV\" $P_WITH ls $REV",
        (! defined $OPT_Hfile ? "\n" :
         ", \\\n \"$OPT_Hfile\" title \"HLT\" w lp ls $HLT");
    
    #-- interactive mode
    if ( $OPT_terminal eq $X11 ) {
        print GFILE "\n",
        "print \"-- INTERACTIVE MODE --\"\n",
        "print \"consult gnuplot docs for command list\"\n",
        "print \"mouse 1: coords to clipboard\"\n",
        "print \"mouse 2: mark on plot\"\n",
        "print \"mouse 3: zoom box\"\n",
        "print \"'h' for help in plot window\"\n",
        "print \"enter to exit\"\n",
        "pause -1\n";
    }

    close (GFILE)
        or print STDERR "WARNING: Trouble closing $OPT_Gfile, $!\n";
}


#------------------------------------------------------------------- RunGP ----#
sub RunGP ( )
{
    if ( defined $OPT_Pfile ) {
        print STDERR "Rendering plot $OPT_Pfile\n";
    }
    else {
        print STDERR "Rendering plot to screen\n";
    }

    my $cmd = $GNUPLOT_EXE;

    #-- x11 specifics
    if ( $OPT_terminal eq $X11 ) {
        my $size = $TERMSIZE{$OPT_terminal}{$OPT_size};
        $cmd .= " -geometry ${size}x";
        if ( $OPT_coverage ) { $size = sprintf ("%.0f", $size * .375); }
        $cmd .= "${size}+0+0 -title mummerplot";

        if ( defined $OPT_color && $OPT_color == 0 ) {
            $cmd .= " -mono";
            $cmd .= " -xrm 'gnuplot*line1Dashes: 0'";
            $cmd .= " -xrm 'gnuplot*line2Dashes: 0'";
            $cmd .= " -xrm 'gnuplot*line3Dashes: 0'";
        }

        if ( $OPT_rv ) {
            $cmd .= " -rv";
            $cmd .= " -xrm 'gnuplot*background: black'";
            $cmd .= " -xrm 'gnuplot*textColor: white'";
            $cmd .= " -xrm 'gnuplot*borderColor: white'";
            $cmd .= " -xrm 'gnuplot*axisColor: white'";
        }
    }

    $cmd .= " $OPT_Gfile";

    system ($cmd)
        and print STDERR "WARNING: Unable to run '$cmd', $!\n";
}


#---------------------------------------------------------------- ListenGP ----#
sub ListenGP($$)
{
    my $rref = shift;
    my $qref = shift;

    my ($refc, $qryc);
    my ($refid, $qryid);
    my ($rsock, $qsock);
    my $oldclip = "";

    #-- get IDs sorted by offset
    my @refo = sort { $rref->{$a}[0] <=> $rref->{$b}[0] } keys %$rref;
    my @qryo = sort { $qref->{$a}[0] <=> $qref->{$b}[0] } keys %$qref;

    #-- attempt to connect sockets
    if ( $OPT_rport ) {
        $rsock = IO::Socket::INET->new("localhost:$OPT_rport")
            or print STDERR "WARNING: Could not connect to rport $OPT_rport\n";
    }

    if ( $OPT_qport ) {
        $qsock = IO::Socket::INET->new("localhost:$OPT_qport")
            or print STDERR "WARNING: Could not connect to qport $OPT_qport\n";
    }

    #-- while parent still exists
    while ( getppid != 1 ) {

        #-- query the clipboard
        $_ = `xclip -o -silent -selection primary`;
        if ( $? >> 8 ) {
            die "WARNING: Unable to query clipboard with xclip\n";
        }

        #-- if cliboard has changed and contains a coordinate
        if ( $_ ne $oldclip && (($refc, $qryc) = /^\[(\d+), (\d+)\]/) ) {

            $oldclip = $_;

            #-- translate the reference position
            $refid = "NULL";
            for ( my $i = 0; $i < (scalar @refo); ++ $i ) {
                my $aref = $rref->{$refo[$i]};
                if ( $i == $#refo || $aref->[0] + $aref->[1] > $refc ) {
                    $refid = $refo[$i];
                    $refc -= $aref->[0];
                    if ( $aref->[2] == -1 ) {
                        $refc = $aref->[1] - $refc + 1;
                    }
                    last;
                }
            }

            #-- translate the query position
            $qryid = "NULL";
            for ( my $i = 0; $i < (scalar @qryo); ++ $i ) {
                my $aref = $qref->{$qryo[$i]};
                if ( $i == $#qryo || $aref->[0] + $aref->[1] > $qryc ) {
                    $qryid = $qryo[$i];
                    $qryc -= $aref->[0];
                    if ( $aref->[2] == -1 ) {
                        $qryc = $aref->[1] - $qryc + 1;
                    }
                    last;
                }
            }

            #-- print the info to stdout and socket
            print "$refid\t$qryid\t$refc\t$qryc\n";

            if ( $rsock ) {
                print $rsock "contig I$refid $refc\n";
                print "sent \"contig I$refid $refc\" to $OPT_rport\n";
            }
            if ( $qsock ) {
                print $qsock "contig I$qryid $qryc\n";
                print "sent \"contig I$qryid $qryc\" to $OPT_qport\n";
            }
        }

        #-- sleep for half second
        select undef, undef, undef, .5;
    }

    exit (0);
}


#------------------------------------------------------------ ParseOptions ----#
sub ParseOptions ( )
{
    my ($opt_small, $opt_medium, $opt_large);
    my ($opt_ps, $opt_x11, $opt_png);
    my $cnt;

    #-- Get options
    my $err = $tigr -> TIGR_GetOptions
        (
         "b|breaklen:i" => \$OPT_breaklen,
         "color!"       => \$OPT_color,
         "c|coverage!"  => \$OPT_coverage,
         "f|filter!"    => \$OPT_filter,
         "l|layout!"    => \$OPT_layout,
         "p|prefix=s"   => \$OPT_prefix,
         "rv"           => \$OPT_rv,
         "r|IdR=s"      => \$OPT_IdR,
         "q|IdQ=s"      => \$OPT_IdQ,
         "R|Rfile=s"    => \$OPT_IDRfile,
         "Q|Qfile=s"    => \$OPT_IDQfile,
         "rport=i"      => \$OPT_rport,
         "qport=i"      => \$OPT_qport,
         "s|size=s"     => \$OPT_size,
         "S|SNP"        => \$OPT_SNP,
         "t|terminal=s" => \$OPT_terminal,
         "title=s"      => \$OPT_title,
         "x|xrange=s"   => \$OPT_xrange,
         "y|yrange=s"   => \$OPT_yrange,
         "x11"          => \$opt_x11,
         "postscript"   => \$opt_ps,
         "png"          => \$opt_png,
         "small"        => \$opt_small,
         "medium"       => \$opt_medium,
         "large"        => \$opt_large,
         "fat"          => \$OPT_ONLY_USE_FATTEST,
         );

    if ( !$err  ||  scalar (@ARGV) != 1 ) {
        $tigr -> printUsageInfo( );
        die "Try '$0 -h' for more information.\n";
    }

    $cnt = 0;
    if ( $opt_png ) { $OPT_terminal = $PNG; $cnt ++; }
    if ( $opt_ps  ) { $OPT_terminal = $PS;  $cnt ++; }
    if ( $opt_x11 ) { $OPT_terminal = $X11; $cnt ++; }
    if ( $cnt > 1 ) {
        print STDERR
            "WARNING: Multiple terminals not allowed, using '$OPT_terminal'\n";
    }

    $cnt = 0;
    if ( $opt_large  ) { $OPT_size = $LARGE;  $cnt ++; }
    if ( $opt_medium ) { $OPT_size = $MEDIUM; $cnt ++; }
    if ( $opt_small  ) { $OPT_size = $SMALL;  $cnt ++; }
    if ( $cnt > 1 ) {
        print STDERR
            "WARNING: Multiple sizes now allowed, using '$OPT_size'\n";
    }

    #-- Check that status of gnuplot
    $OPT_gpstatus = system ("gnuplot --version");

    if ( $OPT_gpstatus == -1 ) {
        print STDERR
            "WARNING: Could not find gnuplot, plot will not be rendered\n";
    }
    elsif ( $OPT_gpstatus ) {
        print STDERR
            "WARNING: Using outdated gnuplot, use v4 or later for best results\n";

        if ( $OPT_color ) {
            print STDERR
                "WARNING: Turning off --color option for compatibility\n";
            undef $OPT_color;
        }

        if ( $OPT_terminal eq $PNG  &&  $OPT_size ne $SMALL ) { 
            print STDERR
                "WARNING: Turning off --size option for compatibility\n";
            $OPT_size = $SMALL;
        }
    }

    #-- Check options
    if ( !exists $TERMSIZE{$OPT_terminal} ) {
        die "ERROR: Invalid terminal type, $OPT_terminal\n";
    }

    if ( !exists $TERMSIZE{$OPT_terminal}{$OPT_size} ) {
        die "ERROR: Invalid terminal size, $OPT_size\n";
    }

    if ( $OPT_xrange ) {
        $OPT_xrange =~ tr/,/:/;
        $OPT_xrange =~ /^\[\d+:\d+\]$/
            or die "ERROR: Invalid xrange format, $OPT_xrange\n";
    }

    if ( $OPT_yrange ) {
        $OPT_yrange =~ tr/,/:/;
        $OPT_yrange =~ /^\[\d+:\d+\]$/
            or die "ERROR: Invalid yrange format, $OPT_yrange\n";
    }

    #-- Set file names
    $OPT_Mfile = $ARGV[0];
    $tigr->isReadableFile ($OPT_Mfile)
        or die "ERROR: Could not read $OPT_Mfile, $!\n";

    $OPT_Ffile = $OPT_prefix . $SUFFIX{$FWDPLOT};
    $tigr->isWritableFile ($OPT_Ffile) or $tigr->isCreatableFile ($OPT_Ffile)
        or die "ERROR: Could not write $OPT_Ffile, $!\n";

    $OPT_Rfile = $OPT_prefix . $SUFFIX{$REVPLOT};
    $tigr->isWritableFile ($OPT_Rfile) or $tigr->isCreatableFile ($OPT_Rfile)
        or die "ERROR: Could not write $OPT_Rfile, $!\n";

    if ( defined $OPT_breaklen || defined $OPT_SNP ) {
        $OPT_Hfile = $OPT_prefix . $SUFFIX{$HLTPLOT};
        $tigr->isWritableFile($OPT_Hfile) or $tigr->isCreatableFile($OPT_Hfile)
            or die "ERROR: Could not write $OPT_Hfile, $!\n";
    }

    if ($OPT_ONLY_USE_FATTEST)
    {
      $OPT_layout = 1;
    }

    if ( $OPT_filter || $OPT_layout ) {
        $OPT_Dfile = $OPT_prefix . $SUFFIX{$FILTER};
        $tigr->isWritableFile($OPT_Dfile) or $tigr->isCreatableFile($OPT_Dfile)
            or die "ERROR: Could not write $OPT_Dfile, $!\n";
    }

    $OPT_Gfile = $OPT_prefix . $SUFFIX{$GNUPLOT};
    $tigr->isWritableFile ($OPT_Gfile) or $tigr->isCreatableFile ($OPT_Gfile)
        or die "ERROR: Could not write $OPT_Gfile, $!\n";

    if ( exists $SUFFIX{$OPT_terminal} ) {
        $OPT_Pfile = $OPT_prefix . $SUFFIX{$OPT_terminal};
        $tigr->isWritableFile($OPT_Pfile) or $tigr->isCreatableFile($OPT_Pfile)
            or die "ERROR: Could not write $OPT_Pfile, $!\n";
    }

    if ( defined $OPT_IDRfile ) {
        $tigr->isReadableFile ($OPT_IDRfile)
            or die "ERROR: Could not read $OPT_IDRfile, $!\n";
    }

    if ( defined $OPT_IDQfile ) {
        $tigr->isReadableFile ($OPT_IDQfile)
            or die "ERROR: Could not read $OPT_IDQfile, $!\n";
    }

    if ( (defined $OPT_rport || defined $OPT_qport) &&
         ($OPT_terminal ne $X11 || $OPT_gpstatus ) ) {
        print STDERR
            "WARNING: Port options available only for v4.0 X11 plots\n";
        undef $OPT_rport;
        undef $OPT_qport;
    }


    if ( defined $OPT_color && defined $OPT_Hfile ) {
        print STDERR
            "WARNING: Turning off --color option so highlighting is visible\n";
        undef $OPT_color;
    }
}
