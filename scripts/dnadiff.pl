#!@PERL@

#-------------------------------------------------------------------------------
#   Programmer: Adam M Phillippy, University of Maryland
#         File: dnadiff
#         Date: 11 / 29 / 06
#
#   Try 'dnadiff -h' for more information.
#
#-------------------------------------------------------------------------------

use lib "@LIB_DIR@";
use Foundation;
use File::Spec::Functions;
use warnings;
use strict;

my $BIN_DIR = "@BIN_DIR@";
my $SCRIPT_DIR = "@LIB_DIR@";

my $VERSION_INFO = q~
DNAdiff version 1.3
    ~;

my $HELP_INFO = q~
  USAGE: dnadiff  [options]  <reference>  <query>
    or   dnadiff  [options]  -d <delta file>

  DESCRIPTION:
    Run comparative analysis of two sequence sets using nucmer and its
    associated utilities with recommended parameters. See MUMmer
    documentation for a more detailed description of the
    output. Produces the following output files:

    .report  - Summary of alignments, differences and SNPs
    .delta   - Standard nucmer alignment output
    .1delta  - 1-to-1 alignment from delta-filter -1
    .mdelta  - M-to-M alignment from delta-filter -m
    .1coords - 1-to-1 coordinates from show-coords -THrcl .1delta
    .mcoords - M-to-M coordinates from show-coords -THrcl .mdelta
    .snps    - SNPs from show-snps -rlTHC .1delta
    .rdiff   - Classified ref breakpoints from show-diff -rH .mdelta
    .qdiff   - Classified qry breakpoints from show-diff -qH .mdelta
    .unref   - Unaligned reference IDs and lengths (if applicable)
    .unqry   - Unaligned query IDs and lengths (if applicable)

  MANDATORY:
    reference       Set the input reference multi-FASTA filename
    query           Set the input query multi-FASTA filename
      or
    delta file      Unfiltered .delta alignment file from nucmer

  OPTIONS:
    -d|delta        Provide precomputed delta file for analysis
    -h
    --help          Display help information and exit
    -p|prefix       Set the prefix of the output files (default "out")
    -V
    --version       Display the version information and exit
    ~;


my $USAGE_INFO = q~
  USAGE: dnadiff  [options]  <reference>  <query>
    or   dnadiff  [options]  -d <delta file>
    ~;


my @DEPEND_INFO =
    (
     "$BIN_DIR/delta-filter",
     "$BIN_DIR/show-diff",
     "$BIN_DIR/show-snps",
     "$BIN_DIR/show-coords",
     "$BIN_DIR/nucmer",
     "$SCRIPT_DIR/Foundation.pm"
     );

my $DELTA_FILTER = "$BIN_DIR/delta-filter";
my $SHOW_DIFF = "$BIN_DIR/show-diff";
my $SHOW_SNPS = "$BIN_DIR/show-snps";
my $SHOW_COORDS = "$BIN_DIR/show-coords";
my $NUCMER = "$BIN_DIR/nucmer";

my $SNPBuff         = 20;            # required buffer around "good" snps
my $OPT_Prefix      = "out";         # prefix for all output files
my $OPT_RefFile;                     # reference file
my $OPT_QryFile;                     # query file
my $OPT_DeltaFile;                   # unfiltered alignment file
my $OPT_ReportFile  = ".report";     # report file
my $OPT_DeltaFile1  = ".1delta";     # 1-to-1 delta alignment
my $OPT_DeltaFileM  = ".mdelta";     # M-to-M delta alignment
my $OPT_CoordsFile1 = ".1coords";    # 1-to-1 alignment coords
my $OPT_CoordsFileM = ".mcoords";    # M-to-M alignment coords
my $OPT_SnpsFile    = ".snps";       # snps output file
my $OPT_DiffRFile   = ".rdiff";      # diffile for R
my $OPT_DiffQFile   = ".qdiff";      # diffile for Q
my $OPT_UnRefFile    = ".unref";     # unaligned ref IDs and lengths
my $OPT_UnQryFile    = ".unqry";     # unaligned qry IDs and lengths

my $TIGR;  # TIGR Foundation object


sub RunAlignment();
sub RunFilter();
sub RunCoords();
sub RunSNPs();
sub RunDiff();
sub MakeReport();

sub FastaSizes($$);

sub FileOpen($$);
sub FileClose($$);

sub GetOpt();


#--------------------------------------------------------------------- main ----
 main:
{
    GetOpt();

    RunAlignment() unless defined($OPT_DeltaFile);
    RunFilter();
    RunCoords();
    RunSNPs();
    RunDiff();
    MakeReport();

    exit(0);
}


#------------------------------------------------------------- RunAlignment ----
# Run nucmer
sub RunAlignment()
{
    print STDERR "Building alignments\n";
    my $cmd = "$NUCMER --maxmatch -p $OPT_Prefix $OPT_RefFile $OPT_QryFile";
    my $err = "ERROR: Failed to run nucmer, aborting.\n";

    system($cmd) == 0 or die $err;
    $OPT_DeltaFile = $OPT_Prefix . ".delta";
}


#---------------------------------------------------------------- RunFilter ----
# Run delta-filter
sub RunFilter()
{
    print STDERR "Filtering alignments\n";
    my $cmd1 = "$DELTA_FILTER -1 $OPT_DeltaFile > $OPT_DeltaFile1";
    my $cmd2 = "$DELTA_FILTER -m $OPT_DeltaFile > $OPT_DeltaFileM";
    my $err = "ERROR: Failed to run delta-filter, aborting.\n";

    system($cmd1) == 0 or die $err;
    system($cmd2) == 0 or die $err;
}


#------------------------------------------------------------------ RunSNPs ----
# Run show-snps
sub RunSNPs()
{
    print STDERR "Analyzing SNPs\n";
    my $cmd = "$SHOW_SNPS -rlTHC $OPT_DeltaFile1 > $OPT_SnpsFile";
    my $err = "ERROR: Failed to run show-snps, aborting.\n";

    system($cmd) == 0 or die $err;
}


#---------------------------------------------------------------- RunCoords ----
# Run show-coords
sub RunCoords()
{
    print STDERR "Extracting alignment coordinates\n";
    my $cmd1 = "$SHOW_COORDS -rclTH $OPT_DeltaFile1 > $OPT_CoordsFile1";
    my $cmd2 = "$SHOW_COORDS -rclTH $OPT_DeltaFileM > $OPT_CoordsFileM";
    my $err = "ERROR: Failed to run show-coords, aborting.\n";

    system($cmd1) == 0 or die $err;
    system($cmd2) == 0 or die $err;
}


#------------------------------------------------------------------ RunDiff ----
# Run show-diff
sub RunDiff()
{
    print STDERR "Extracting alignment breakpoints\n";
    my $cmd1 = "$SHOW_DIFF -rH $OPT_DeltaFileM > $OPT_DiffRFile";
    my $cmd2 = "$SHOW_DIFF -qH $OPT_DeltaFileM > $OPT_DiffQFile";
    my $err = "ERROR: Failed to run show-diff, aborting.\n";

    system($cmd1) == 0 or die $err;
    system($cmd2) == 0 or die $err;
}


#--------------------------------------------------------------- MakeReport ----
# Output alignment report
sub MakeReport()
{
    print STDERR "Generating report file\n";

    my ($fhi, $fho);    # filehandle-in and filehandle-out
    my (%refs, %qrys) = ((),());            # R and Q ID->length
    my ($rqnAligns1, $rqnAlignsM) = (0,0);  # alignment counter
    my ($rSumLen1, $qSumLen1) = (0,0);      # alignment length sum
    my ($rSumLenM, $qSumLenM) = (0,0);      # alignment length sum
    my ($rqSumLen1, $rqSumLenM) = (0,0);    # combined alignment length sum
    my ($rqSumIdy1, $rqSumIdyM) = (0,0);    # weighted alignment identity sum
    my ($qnIns, $rnIns) = (0,0);            # insertion count
    my ($qSumIns, $rSumIns) = (0,0);        # insertion length sum
    my ($qnTIns, $rnTIns) = (0,0);          # tandem insertion count
    my ($qSumTIns, $rSumTIns) = (0,0);      # tandem insertion length sum
    my ($qnInv, $rnInv) = (0,0);            # inversion count
    my ($qnRel, $rnRel) = (0,0);            # relocation count
    my ($qnTrn, $rnTrn) = (0,0);            # translocation count
    my ($rnSeqs, $qnSeqs) = (0,0);          # sequence count
    my ($rnASeqs, $qnASeqs) = (0,0);        # aligned sequence count
    my ($rnBases, $qnBases) = (0,0);        # bases count
    my ($rnABases, $qnABases) = (0,0);      # aligned bases count
    my ($rnBrk, $qnBrk) = (0,0);            # breakpoint count
    my ($rqnSNPs, $rqnIndels) = (0,0);      # snp and indel counts
    my ($rqnGSNPs, $rqnGIndels) = (0,0);    # good snp and indel counts
    my %rqSNPs =                            # SNP hash
      ( "."=>{"A"=>0,"C"=>0,"G"=>0,"T"=>0},
        "A"=>{"."=>0,"C"=>0,"G"=>0,"T"=>0},
        "C"=>{"."=>0,"A"=>0,"G"=>0,"T"=>0},
        "G"=>{"."=>0,"A"=>0,"C"=>0,"T"=>0},
        "T"=>{"."=>0,"A"=>0,"C"=>0,"G"=>0} );
    my %rqGSNPs =                           # good SNP hash
      ( "."=>{"A"=>0,"C"=>0,"G"=>0,"T"=>0},
        "A"=>{"."=>0,"C"=>0,"G"=>0,"T"=>0},
        "C"=>{"."=>0,"A"=>0,"G"=>0,"T"=>0},
        "G"=>{"."=>0,"A"=>0,"C"=>0,"T"=>0},
        "T"=>{"."=>0,"A"=>0,"C"=>0,"G"=>0} );

    my $header;                             # delta header

    #-- Get delta header
    $fhi = FileOpen("<", $OPT_DeltaFile);
    $header .= <$fhi>;
    $header .= <$fhi>;
    $header .= "\n";
    FileClose($fhi, $OPT_DeltaFile);

    #-- Collect all reference and query IDs and lengths
    FastaSizes($OPT_RefFile, \%refs);
    FastaSizes($OPT_QryFile, \%qrys);

    #-- Count ref and qry seqs and lengths
    foreach ( values(%refs) ) {
        $rnSeqs++;
        $rnBases += $_;
    }
    foreach ( values(%qrys) ) {
        $qnSeqs++;
        $qnBases += $_;
    }

    #-- Count aligned seqs, aligned bases, and breakpoints for each R and Q
    $fhi = FileOpen("<", $OPT_CoordsFileM);
    while (<$fhi>) {
        chomp;
        my @A = split "\t";
        scalar(@A) == 13
            or die "ERROR: Unrecognized format $OPT_CoordsFileM, aborting.\n";

        #-- Add to M-to-M alignment counts
        $rqnAlignsM++;
        $rSumLenM += $A[4];
        $qSumLenM += $A[5];
        $rqSumIdyM += ($A[6] / 100.0) * ($A[4] + $A[5]);
        $rqSumLenM += ($A[4] + $A[5]);

        #-- If new ID, add to sequence and base count
        if ( $refs{$A[11]} > 0 ) {
            $rnASeqs++;
            $rnABases += $refs{$A[11]};
            $refs{$A[11]} *= -1; # If ref has alignment, length will be -neg
        }
        if ( $qrys{$A[12]} > 0 ) {
            $qnASeqs++;
            $qnABases += $qrys{$A[12]};
            $qrys{$A[12]} *= -1; # If qry has alignment, length will be -neg
        }

        #-- Add to breakpoint counts
        my ($lo, $hi);
        if ( $A[0] < $A[1] ) { $lo = $A[0]; $hi = $A[1]; }
        else                 { $lo = $A[1]; $hi = $A[0]; }
        $rnBrk++ if ( $lo != 1 );
        $rnBrk++ if ( $hi != $A[7] );

        if ( $A[2] < $A[3] ) { $lo = $A[2]; $hi = $A[3]; }
        else                 { $lo = $A[3]; $hi = $A[2]; }
        $qnBrk++ if ( $lo != 1 );
        $qnBrk++ if ( $hi != $A[8] );
    }
    FileClose($fhi, $OPT_CoordsFileM);

    #-- Calculate average %idy, length, etc.
    $fhi = FileOpen("<", $OPT_CoordsFile1);
    while (<$fhi>) {
        chomp;
        my @A = split "\t";
        scalar(@A) == 13
            or die "ERROR: Unrecognized format $OPT_CoordsFile1, aborting.\n";

        #-- Add to 1-to-1 alignment counts
        $rqnAligns1++;
        $rSumLen1 += $A[4];
        $qSumLen1 += $A[5];
        $rqSumIdy1 += ($A[6] / 100.0) * ($A[4] + $A[5]);
        $rqSumLen1 += ($A[4] + $A[5]);
    }
    FileClose($fhi, $OPT_CoordsFile1);

    #-- If you are reading this, you need to get out more...

    #-- Count reference diff features and indels
    $fhi = FileOpen("<", $OPT_DiffRFile);
    while (<$fhi>) {
        chomp;
        my @A = split "\t";
        defined($A[4])
            or die "ERROR: Unrecognized format $OPT_DiffRFile, aborting.\n";
        my $gap = $A[4];
        my $ins = $gap;

        #-- Add to tandem insertion counts
        if ( $A[1] eq "GAP" ) {
            scalar(@A) == 7
                or die "ERROR: Unrecognized format $OPT_DiffRFile, aborting.\n";
            $ins = $A[6] if ( $A[6] > $gap );
            if ( $A[4] <= 0 && $A[5] <= 0 && $A[6] > 0 ) {
                $rnTIns++;
                $rSumTIns += $A[6];
             }
        }

        #-- Remove unaligned sequence from count
        if ( $A[1] ne "DUP" ) {
          $rnABases -= $gap if ( $gap > 0 );
        }

        #-- Add to insertion count
        if ( $ins > 0 ) {
            $rnIns++;
            $rSumIns += $ins;
        }

        #-- Add to rearrangement counts
        $rnInv++ if ( $A[1] eq "INV" );
        $rnRel++ if ( $A[1] eq "JMP" );
        $rnTrn++ if ( $A[1] eq "SEQ" );
    }
    FileClose($fhi, $OPT_DiffRFile);

    #-- Count query diff features and indels
    $fhi = FileOpen("<", $OPT_DiffQFile);
    while (<$fhi>) {
        chomp;
        my @A = split "\t";
        defined($A[4])
            or die "ERROR: Unrecognized format $OPT_DiffRFile, aborting.\n";
        my $gap = $A[4];
        my $ins = $gap;

        #-- Add to tandem insertion counts
        if ( $A[1] eq "GAP" ) {
            scalar(@A) == 7
                or die "ERROR: Unrecognized format $OPT_DiffRFile, aborting.\n";
            $ins = $A[6] if ( $A[6] > $gap );
            if ( $A[4] <= 0 && $A[5] <= 0 && $A[6] > 0 ) {
                $qnTIns++;
                $qSumTIns += $A[6];
            }
        }

        #-- Remove unaligned sequence from count
        if ( $A[1] ne "DUP" ) {
          $qnABases -= $gap if ( $gap > 0 );
        }

        #-- Add to insertion count
        if ( $ins > 0 ) {
            $qnIns++;
            $qSumIns += $ins;
        }

        #-- Add to rearrangement counts
        $qnInv++ if ( $A[1] eq "INV" );
        $qnRel++ if ( $A[1] eq "JMP" );
        $qnTrn++ if ( $A[1] eq "SEQ" );
    }
    FileClose($fhi, $OPT_DiffQFile);

    #-- Count SNPs
    $fhi = FileOpen("<", $OPT_SnpsFile);
    while(<$fhi>) {
        chomp;
        my @A = split "\t";
        scalar(@A) == 12
            or die "ERROR: Unrecognized format $OPT_SnpsFile, aborting\n";

        my $r = uc($A[1]);
        my $q = uc($A[2]);

        #-- Plain SNPs
        $rqSNPs{$r}{$q}++;
        if ( !exists($rqSNPs{$q}{$r}) ) { $rqSNPs{$q}{$r} = 0; }
        if ( $r eq '.' || $q eq '.' ) { $rqnIndels++; }
        else                          { $rqnSNPs++; }

        #-- Good SNPs with sufficient match buffer
        if ( $A[4] >= $SNPBuff ) {
            $rqGSNPs{$r}{$q}++;
            if ( !exists($rqGSNPs{$q}{$r}) ) { $rqGSNPs{$q}{$r} = 0; }
            if ( $r eq '.' || $q eq '.' ) { $rqnGIndels++; }
            else                          { $rqnGSNPs++; }
        }
    }
    FileClose($fhi, $OPT_SnpsFile);


    #-- Output report
    $fho = FileOpen(">", $OPT_ReportFile);

    print  $fho $header;
    printf $fho "%-15s %20s %20s\n", "", "[REF]", "[QRY]";

    print  $fho "[Sequences]\n";

    printf $fho "%-15s %20d %20d\n",
    "TotalSeqs", $rnSeqs, $qnSeqs;
    printf $fho "%-15s %20s %20s\n",
    "AlignedSeqs",
    ( sprintf "%10d(%.2f%%)",
      $rnASeqs, ($rnSeqs ? $rnASeqs / $rnSeqs * 100.0 : 0) ),
    ( sprintf "%10d(%.2f%%)",
      $qnASeqs, ($rnSeqs ? $qnASeqs / $qnSeqs * 100.0 : 0) );
    printf $fho "%-15s %20s %20s\n",
    "UnalignedSeqs",
     ( sprintf "%10d(%.2f%%)",
       $rnSeqs - $rnASeqs,
       ($rnSeqs ? ($rnSeqs - $rnASeqs) / $rnSeqs * 100.0 : 0) ),
     ( sprintf "%10d(%.2f%%)",
       $qnSeqs - $qnASeqs,
       ($qnSeqs ? ($qnSeqs - $qnASeqs) / $qnSeqs * 100.0 : 0) );

    print  $fho "\n[Bases]\n";

    printf $fho "%-15s %20d %20d\n",
    "TotalBases", $rnBases, $qnBases;
    printf $fho "%-15s %20s %20s\n",
    "AlignedBases",
    ( sprintf "%10d(%.2f%%)",
      $rnABases, ($rnBases ? $rnABases / $rnBases * 100.0 : 0) ),
    ( sprintf "%10d(%.2f%%)",
      $qnABases, ($qnBases ? $qnABases / $qnBases * 100.0 : 0) );
    printf $fho "%-15s %20s %20s\n",
    "UnalignedBases",
    ( sprintf "%10d(%.2f%%)",
      $rnBases - $rnABases,
      ($rnBases ? ($rnBases - $rnABases) / $rnBases * 100.0 : 0) ),
    ( sprintf "%10d(%.2f%%)",
      $qnBases - $qnABases,
      ($qnBases ? ($qnBases - $qnABases) / $qnBases * 100.0 : 0) );

    print  $fho "\n[Alignments]\n";

    printf $fho "%-15s %20d %20d\n",
    "1-to-1", $rqnAligns1, $rqnAligns1;
    printf $fho "%-15s %20d %20d\n",
    "TotalLength", $rSumLen1, $qSumLen1;
    printf $fho "%-15s %20.2f %20.2f\n",
    "AvgLength",
    ($rqnAligns1 ? $rSumLen1 / $rqnAligns1 : 0),
    ($rqnAligns1 ? $qSumLen1 / $rqnAligns1 : 0);
    printf $fho "%-15s %20.2f %20.2f\n",
    "AvgIdentity",
    ($rqSumLen1 ? $rqSumIdy1 / $rqSumLen1 * 100.0 : 0),
    ($rqSumLen1 ? $rqSumIdy1 / $rqSumLen1 * 100.0 : 0);

    print  $fho "\n";

    printf $fho "%-15s %20d %20d\n",
    "M-to-M", $rqnAlignsM, $rqnAlignsM;
    printf $fho "%-15s %20d %20d\n",
    "TotalLength", $rSumLenM, $qSumLenM;
    printf $fho "%-15s %20.2f %20.2f\n",
    "AvgLength",
    ($rqnAlignsM ? $rSumLenM / $rqnAlignsM : 0),
    ($rqnAlignsM ? $qSumLenM / $rqnAlignsM : 0);
    printf $fho "%-15s %20.2f %20.2f\n",
    "AvgIdentity",
    ($rqSumLenM ? $rqSumIdyM / $rqSumLenM * 100.0 : 0),
    ($rqSumLenM ? $rqSumIdyM / $rqSumLenM * 100.0 : 0);

    print  $fho "\n[Feature Estimates]\n";

    printf $fho "%-15s %20d %20d\n",
    "Breakpoints", $rnBrk, $qnBrk;
    printf $fho "%-15s %20d %20d\n",
    "Relocations", $rnRel, $qnRel;
    printf $fho "%-15s %20d %20d\n",
    "Translocations", $rnTrn, $qnTrn;
    printf $fho "%-15s %20d %20d\n",
    "Inversions", $rnInv, $qnInv;

    print  $fho "\n";

    printf $fho "%-15s %20d %20d\n",
    "Insertions", $rnIns, $qnIns;
    printf $fho "%-15s %20d %20d\n",
    "InsertionSum", $rSumIns, $qSumIns;
    printf $fho "%-15s %20.2f %20.2f\n",
    "InsertionAvg",
    ($rnIns ? $rSumIns / $rnIns : 0),
    ($qnIns ? $qSumIns / $qnIns : 0);

    print  $fho "\n";

    printf $fho "%-15s %20d %20d\n",
    "TandemIns", $rnTIns, $qnTIns;
    printf $fho "%-15s %20d %20d\n",
    "TandemInsSum", $rSumTIns, $qSumTIns;
    printf $fho "%-15s %20.2f %20.2f\n",
    "TandemInsAvg",
    ($rnTIns ? $rSumTIns / $rnTIns : 0),
    ($qnTIns ? $qSumTIns / $qnTIns : 0);

    print  $fho "\n[SNPs]\n";

    printf $fho "%-15s %20d %20d\n",
    "TotalSNPs", $rqnSNPs, $rqnSNPs;
    foreach my $r (keys %rqSNPs) {
      foreach my $q (keys %{$rqSNPs{$r}}) {
        if ( $r ne "." && $q ne "." ) {
          printf $fho "%-15s %20s %20s\n",
            "$r$q",
              ( sprintf "%10d(%.2f%%)",
                $rqSNPs{$r}{$q},
                ($rqnSNPs ? $rqSNPs{$r}{$q} / $rqnSNPs * 100.0 : 0) ),
                  ( sprintf "%10d(%.2f%%)",
                    $rqSNPs{$q}{$r},
                    ($rqnSNPs ? $rqSNPs{$q}{$r} / $rqnSNPs * 100.0 : 0) );
        }
      }
    }

    print  $fho "\n";

    printf $fho "%-15s %20d %20d\n",
    "TotalGSNPs", $rqnGSNPs, $rqnGSNPs;
    foreach my $r (keys %rqGSNPs) {
      foreach my $q (keys %{$rqGSNPs{$r}}) {
        if ( $r ne "." && $q ne "." ) {
          printf $fho "%-15s %20s %20s\n",
            "$r$q",
              ( sprintf "%10d(%.2f%%)",
                $rqGSNPs{$r}{$q},
                ($rqnGSNPs ? $rqGSNPs{$r}{$q} / $rqnGSNPs * 100.0 : 0) ),
                  ( sprintf "%10d(%.2f%%)",
                    $rqGSNPs{$q}{$r},
                    ($rqnGSNPs ? $rqGSNPs{$q}{$r} / $rqnGSNPs * 100.0 : 0) );
        }
      }
    }

    print  $fho "\n";

    printf $fho "%-15s %20d %20d\n",
    "TotalIndels", $rqnIndels, $rqnIndels;
    foreach my $r (keys %rqSNPs) {
      foreach my $q (keys %{$rqSNPs{$r}}) {
        if ( $q eq "." ) {
          printf $fho "%-15s %20s %20s\n",
            "$r$q",
              ( sprintf "%10d(%.2f%%)",
                $rqSNPs{$r}{$q},
                ($rqnIndels ? $rqSNPs{$r}{$q} / $rqnIndels * 100.0 : 0) ),
                  ( sprintf "%10d(%.2f%%)",
                    $rqSNPs{$q}{$r},
                    ($rqnIndels ? $rqSNPs{$q}{$r} / $rqnIndels * 100.0 : 0) );
        }
      }
    }
    foreach my $r (keys %rqSNPs) {
      foreach my $q (keys %{$rqSNPs{$r}}) {
        if ( $r eq "." ) {
          printf $fho "%-15s %20s %20s\n",
            "$r$q",
              ( sprintf "%10d(%.2f%%)",
                $rqSNPs{$r}{$q},
                ($rqnIndels ? $rqSNPs{$r}{$q} / $rqnIndels * 100.0 : 0) ),
                  ( sprintf "%10d(%.2f%%)",
                    $rqSNPs{$q}{$r},
                    ($rqnIndels ? $rqSNPs{$q}{$r} / $rqnIndels * 100.0 : 0) );
        }
      }
    }

    print  $fho "\n";

    printf $fho "%-15s %20d %20d\n",
    "TotalGIndels", $rqnGIndels, $rqnGIndels;
    foreach my $r (keys %rqGSNPs) {
      foreach my $q (keys %{$rqGSNPs{$r}}) {
        if ( $q eq "." ) {
          printf $fho "%-15s %20s %20s\n",
            "$r$q",
              ( sprintf "%10d(%.2f%%)",
                $rqGSNPs{$r}{$q},
                ($rqnGIndels ? $rqGSNPs{$r}{$q} / $rqnGIndels * 100.0 : 0) ),
                  ( sprintf "%10d(%.2f%%)",
                    $rqGSNPs{$q}{$r},
                    ($rqnGIndels ? $rqGSNPs{$q}{$r} / $rqnGIndels * 100.0 : 0) );
        }
      }
    }
    foreach my $r (keys %rqGSNPs) {
      foreach my $q (keys %{$rqGSNPs{$r}}) {
        if ( $r eq "." ) {
          printf $fho "%-15s %20s %20s\n",
            "$r$q",
              ( sprintf "%10d(%.2f%%)",
                $rqGSNPs{$r}{$q},
                ($rqnGIndels ? $rqGSNPs{$r}{$q} / $rqnGIndels * 100.0 : 0) ),
                  ( sprintf "%10d(%.2f%%)",
                    $rqGSNPs{$q}{$r},
                    ($rqnGIndels ? $rqGSNPs{$q}{$r} / $rqnGIndels * 100.0 : 0) );
        }
      }
    }

    FileClose($fho, $OPT_ReportFile);


    #-- Output unaligned reference and query IDs, if applicable
    if ( $rnSeqs != $rnASeqs ) {
        $fho = FileOpen(">", $OPT_UnRefFile);
        while ( my ($key, $val) = each(%refs) ) {
            print $fho "$key\tUNI\t1\t$val\t$val\n" unless $val < 0;
        }
        FileClose($fho, $OPT_UnRefFile);
    }
    if ( $qnSeqs != $qnASeqs ) {
        $fho = FileOpen(">", $OPT_UnQryFile);
        while ( my ($key, $val) = each(%qrys) ) {
            print $fho "$key\tUNI\t1\t$val\t$val\n" unless $val < 0;
        }
        FileClose($fho, $OPT_UnQryFile);
    }
}


#--------------------------------------------------------------- FastaSizes ----
# Compute lengths for a multi-fasta file and store in hash reference
sub FastaSizes($$)
{

    my $file = shift;
    my $href = shift;
    my ($tag, $len);

    my $fhi = FileOpen("<", $file);
    while (<$fhi>) {
        chomp;

        if ( /^>/ ) {
            $href->{$tag} = $len if defined($tag);
            ($tag) = /^>(\S+)/;
            $len = 0;
        } else {
            if ( /\s/ ) {
                die "ERROR: Whitespace found in FastA $file, aborting.\n";
            }
            $len += length;
        }
    }
    $href->{$tag} = $len if defined($tag);
    FileClose($fhi, $file);
}


#----------------------------------------------------------------- FileOpen ----
# Open file, return filehandle, or die
sub FileOpen($$)
{
    my ($mode, $name) = @_;
    my $fhi;
    open($fhi, $mode, $name)
        or die "ERROR: Could not open $name, aborting. $!\n";
    return $fhi;
}


#---------------------------------------------------------------- FileClose ----
# Close file, or die
sub FileClose($$)
{
    my ($fho, $name) = @_;
    close($fho) or die "ERROR: Could not close $name, aborting. $!\n"
}


#------------------------------------------------------------------- GetOpt ----
# Get command options and check file permissions
sub GetOpt()
{
    #-- Initialize TIGR::Foundation
    $TIGR = new TIGR::Foundation;
    if ( !defined($TIGR) ) {
        print STDERR "ERROR: TIGR::Foundation could not be initialized";
        exit(1);
    }

    #-- Set help and usage information
    $TIGR->setHelpInfo($HELP_INFO);
    $TIGR->setUsageInfo($USAGE_INFO);
    $TIGR->setVersionInfo($VERSION_INFO);
    $TIGR->addDependInfo(@DEPEND_INFO);

    #-- Get options
    my $err = !$TIGR->TIGR_GetOptions
        (
         "d|delta=s"  => \$OPT_DeltaFile,
         "p|prefix=s" => \$OPT_Prefix,
         );

    #-- Check if the parsing was successful
    if ( $err
         || (defined($OPT_DeltaFile) && scalar(@ARGV) != 0)
         || (!defined($OPT_DeltaFile) && scalar(@ARGV) != 2) ) {
        $TIGR->printUsageInfo();
        print STDERR "Try '$0 -h' for more information.\n";
        exit(1);
    }

    my @errs;

    $TIGR->isExecutableFile($DELTA_FILTER)
        or push(@errs, $DELTA_FILTER);

    $TIGR->isExecutableFile($SHOW_DIFF)
        or push(@errs, $SHOW_DIFF);

    $TIGR->isExecutableFile($SHOW_SNPS)
        or push(@errs, $SHOW_SNPS);

    $TIGR->isExecutableFile($SHOW_COORDS)
        or push(@errs, $SHOW_COORDS);

    $TIGR->isExecutableFile($NUCMER)
        or push(@errs, $NUCMER);

    if ( defined($OPT_DeltaFile) ) {
        $TIGR->isReadableFile($OPT_DeltaFile)
            or push(@errs, $OPT_DeltaFile);

        my $fhi = FileOpen("<", $OPT_DeltaFile);
        $_ = <$fhi>;
        FileClose($fhi, $OPT_DeltaFile);

        ($OPT_RefFile, $OPT_QryFile) = /^(.+) (.+)$/;
    }
    else {
        $OPT_RefFile = File::Spec->rel2abs($ARGV[0]);
        $OPT_QryFile = File::Spec->rel2abs($ARGV[1]);
    }

    $TIGR->isReadableFile($OPT_RefFile)
        or push(@errs, $OPT_RefFile);

    $TIGR->isReadableFile($OPT_QryFile)
        or push(@errs, $OPT_QryFile);

    $OPT_ReportFile = $OPT_Prefix . $OPT_ReportFile;
    $TIGR->isCreatableFile("$OPT_ReportFile")
        or $TIGR->isWritableFile("$OPT_ReportFile")
        or push(@errs, "$OPT_ReportFile");

    $OPT_DeltaFile1 = $OPT_Prefix . $OPT_DeltaFile1;
    $TIGR->isCreatableFile("$OPT_DeltaFile1")
        or $TIGR->isWritableFile("$OPT_DeltaFile1")
        or push(@errs, "$OPT_DeltaFile1");

    $OPT_DeltaFileM = $OPT_Prefix . $OPT_DeltaFileM;
    $TIGR->isCreatableFile("$OPT_DeltaFileM")
        or $TIGR->isWritableFile("$OPT_DeltaFileM")
        or push(@errs, "$OPT_DeltaFileM");

    $OPT_CoordsFile1 = $OPT_Prefix . $OPT_CoordsFile1;
    $TIGR->isCreatableFile("$OPT_CoordsFile1")
        or $TIGR->isWritableFile("$OPT_CoordsFile1")
        or push(@errs, "$OPT_CoordsFile1");

    $OPT_CoordsFileM = $OPT_Prefix . $OPT_CoordsFileM;
    $TIGR->isCreatableFile("$OPT_CoordsFileM")
        or $TIGR->isWritableFile("$OPT_CoordsFileM")
        or push(@errs, "$OPT_CoordsFileM");

    $OPT_SnpsFile = $OPT_Prefix . $OPT_SnpsFile;
    $TIGR->isCreatableFile("$OPT_SnpsFile")
        or $TIGR->isWritableFile("$OPT_SnpsFile")
        or push(@errs, "$OPT_SnpsFile");

    $OPT_DiffRFile = $OPT_Prefix . $OPT_DiffRFile;
    $TIGR->isCreatableFile("$OPT_DiffRFile")
        or $TIGR->isWritableFile("$OPT_DiffRFile")
        or push(@errs, "$OPT_DiffRFile");

    $OPT_DiffQFile = $OPT_Prefix . $OPT_DiffQFile;
    $TIGR->isCreatableFile("$OPT_DiffQFile")
        or $TIGR->isWritableFile("$OPT_DiffQFile")
        or push(@errs, "$OPT_DiffQFile");

    $OPT_UnRefFile = $OPT_Prefix . $OPT_UnRefFile;
        $TIGR->isCreatableFile("$OPT_UnRefFile")
        or $TIGR->isWritableFile("$OPT_UnRefFile")
        or push(@errs, "$OPT_UnRefFile");

    $OPT_UnQryFile = $OPT_Prefix . $OPT_UnQryFile;
        $TIGR->isCreatableFile("$OPT_UnQryFile")
        or $TIGR->isWritableFile("$OPT_UnQryFile")
        or push(@errs, "$OPT_UnQryFile");

    if ( scalar(@errs) ) {
        print STDERR "ERROR: The following critical files could not be used\n";
        while ( scalar(@errs) ) { print(STDERR pop(@errs),"\n"); }
        print STDERR "Check your paths and file permissions and try again\n";
        exit(1);
    }
}
