# MUMmer4.x README 

**NOTE**

Further documentation, but potentially out of date, is in the [docs](../../tree/master/docs) directory. Please refer to the [INSTALL.md](INSTALL.md) file for installation
instructions.

This file contains brief descriptions of all
executables in the base directory and general information about the
MUMmer package.

## DESCRIPTION
MUMmer is a system for rapidly aligning entire genomes.  The current
version (release 4.x) can find all 20 base pair maximal exact matches between
two bacterial genomes of ~5 million base pairs each in 20 seconds, using 90 MB
of memory, on a typical 1.8 GHz Linux desktop computer.  MUMmer can also align
incomplete genomes; it handles the 100s or 1000s of contigs from a shotgun
sequencing project with ease, and will align them to another set of contigs or
a genome, using the nucmer utility included with the system.  The promer
utility takes this a step further by generating alignments based upon the
six-frame translations of both input sequences.  promer permits the alignment
of genomes for which the proteins are similar but the DNA sequence is too
divergent to detect similarity.  See the nucmer and promer readme files in the
"docs/" subdirectory for more details.  MUMmer is open source, so all we ask
is that you cite our most recent paper in any publications that use this
system:

**(Version 3.0 described)**
 
 >Versatile and open software for comparing large genomes.</br>
 >S. Kurtz, A. Phillippy, A.L. Delcher, M. Smoot, M. Shumway, C. Antonescu, and S.L. Salzberg.</br>
 >Genome Biology (2004), 5:R12.

**(Version 2.1 described)**
>Fast algorithms for large-scale genome alignment and comparison.</br>
>A.L. Delcher. A. Phillippy, J. Carlton, and S.L. Salzberg.</br>
>Nucleic Acids Research 30:11 (2002), 2478-2483.

**(Version 1.0 described)**
>Alignment of Whole Genomes.</br>
>A.L. Delcher, S. Kasif,R.D. Fleischmann, J. Peterson, O. White, and S.L. Salzberg.</br>
>Nucleic Acids Research, 27:11 (1999), 2369-2376.


## RUNNING MUMmer4.x

MUMmer4.x is comprised of many various utilities
and scripts.  For general purposes, the programs `nucmer`, and
`promer` will be all that is needed.  See their descriptions in the
"RUNNING THE MUMmer PROGRAMS" section, or refer to their individual
documentation in the "docs/" subdirectory.  Refer to the "RUNNING THE
MUMmer UTILITIES" section for a brief description of all of the
utilities in this directory.

### Simple use case

Given a file containing a single reference sequence (ref.seq) in
FASTA format and another file containing multiple sequences in FastA
format (qry.seq) type the following at the command line:

`./nucmer  -p <prefix>  ref.seq  qry.seq`

To produce the following files:
   
    <prefix>.delta

Please read the utility-specific documentation in the "docs/" subdirectory
for descriptions of these files and information on how to change the
alignment parameters for the scripts (minimum match length, etc.), or see
the notes below in the "RUNNING THE MUMmer SCRIPTS" section for a brief
explanation.

To see a simple gnuplot output, if you have gnuplot installed, run
the perl script `mummerplot` on the output files. This script can be run
on mummer output (.out), or nucmer/promer output (.delta). Edit the
<prefix>.gp file that is created to change colors, line thicknesses, etc. or
explore the \<prefix>.[fr]plot file to see the data collection.

`./mummerplot  -p <prefix>  <prefix>.out`

## RUNNING THE MUMmer SCRIPTS
Because of MUMmer's modular design, it may be necessary to use a number
of separate programs to produce the desired output.  The MUMmer scripts
attempt to simplify this process by wrapping various utilities into packages
that can perform standard alignment requests.  Listed below are brief
descriptions and usage definitions for these scripts.  Please refer to the
"docs/" subdirectory for a more detailed description of each script.


### nucmer

**Description:**

`nucmer` is for the all-vs-all comparison of nucleotide sequences
contained in multi-FastA data files.  It is best used for highly
similar sequence that may have large rearrangements.  Common use
cases are: comparing two unfinished shotgun sequencing assemblies,
mapping an unfinished sequencing assembly to a finished genome, and
comparing two fairly similar genomes that may have large
rearrangements and duplications.  Please refer to "docs/nucmer.README"
for more information regarding this script and its output, or type
`nucmer -h` for a list of its options.

        USAGE:
        nucmer  [options]  <reference>  <query>

        [options]    type 'nucmer -h' for a list of options.
        <reference>  specifies the multi-FastA sequence file that contains
                     the reference sequences, to be aligned with the queries.
        <query>      specifies the multi-FastA sequence file that contains
                     the query sequences, to be aligned with the references.

        OUTPUT:
        out.delta    the delta encoded alignments between the reference and
                     query sequences.  This file can be parsed with any of
                     the show-* programs which are described in the "RUNNING
                     THE MUMmer UTILITIES" section.

**Notes:**

All output coordinates reference the forward strand of the involved
sequence, regardless of the match direction. Also, nucmer now uses
only matches that are unique in the reference sequence by default,
use the '--mum' or '--maxmatch' options to change this behavior.


### promer

**Description:**

`promer` is for the protein level, all-vs-all comparison of nucleotide
sequences contained in multi-FastA data files.  The nucleotide input
files are translated in all 6 reading frames and then aligned to one
another via the same methods as nucmer.  It is best used for highly
divergent sequences that may have moderate to high similarity on the
protein level.  Common use cases are: identifying syntenic regions
between highly divergent genomes, comparative genome annotation i.e.
using an already annotated genome to help in the annotation of a
newly sequenced genome, and the general comparison of two fairly
divergent genomes that have large rearrangements and may only be
similar on the protein level. Please refer to "docs/promer.README"
for more information regarding this script and its output, or type
`promer -h` for a list of its options.

        USAGE:
        promer  [options]  <reference>  <query>

        [options]    type 'promer -h' for a list of options.
        <reference>  specifies the multi-FastA sequence file that contains
                     the reference sequences, to be aligned with the queries.
        <query>      specifies the multi-FastA sequence file that contains
                     the query sequences, to be aligned with the references.

        OUTPUT:
        out.delta    the delta encoded alignments between the reference and
                     query sequences.  This file can be parsed with any of
                     the show-* programs which are described in the "RUNNING
                     THE MUMmer UTILITIES" section.

**Notes:**

All output coordinates reference the forward strand of the involved
sequence, regardless of the match direction, and are measured in
nucleotides with the exception of the delta integers which are
measured in amino acids (1 delta int = 3 nucleotides). Also, promer
now uses only matches that are unique in the reference sequence by
default, use the '--mum' or '--maxmatch' options to change this
behavior.


**Notes:**

All output coordinates reference their respective strand.  This means
that for all reverse matches, the coordinates that reference the
query sequence will be relative to the reverse complement of the
query sequence.  Please use nucmer or promer if this coordinate
system is confusing.


### dnadiff

**Description:**

This script is a wrapper around nucmer that builds an
alignment using default parameters, and runs many of nucmer's
helper scripts to process the output and report alignment
statistics, SNPs, breakpoints, etc. It is designed for
evaluating the sequence and structural similarity of two
highly similar sequence sets. E.g. comparing two different
assemblies of the same organism, or comparing two strains of
the same species.  Please refer to "docs/dnadiff.README" for
more information regarding this script and its output, or type
'dnadiff -h' for a list of its options.

        USAGE: dnadiff  [options]  <reference>  <query>
          or   dnadiff  [options]  -d <delta file>

        <reference>       Set the input reference multi-FASTA filename
        <query>           Set the input query multi-FASTA filename
           or
        <delta file>      Unfiltered .delta alignment file from nucmer

        OUTPUT:
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

**Notes:**

The report file generated by this script can be useful for
comparing the differences between two similar genomes or
assemblies. The other outputs generated by this script are in
unlabeled tabular format, so please refer to the utility
specific documentation for interpreting them. A full
description of the report file is given in "docs/dnadiff.README".


## RUNNING THE MUMmer UTILITIES
The MUMmer package consists of various utilities that can interact with
the `mummer` program.  `mummer` performs all maximal and maximal unique
matching, and all other utilities were designed to process the input and
output of this program and its related scripts, in order to extract
additional information from the output.  Listed below are the descriptions
and usage definitions for these utilities.


### annotate

**Description:**
This program reads the output of the `gaps` program and adds alignment
information to it.  Part of the original MUMmer1.0 pipeline and can
only be used on the output of the `gaps` program.

        USAGE:
        annotate  <gapsfile>  <seq2>

        <gapsfile>  the output of the 'gaps' program.
        <seq2>      the file containing the second sequence in the comparison.

        OUTPUT:
        stdout           the 'gaps' output interspersed with the alignments of
                         the gaps between adjacent MUMs.  An alignment of a
                         gap comes after the second MUM defining the gap, and
                         alignment errors are marked with a '^' character.
        witherrors.gaps  the 'gaps' output with an appended column that lists
                         the number of alignment errors for each gap.

**Notes:**

This program will eventually be dropped in favor of the combineMUMs
or nucmer match extenders, but persists for the time being.


### combineMUMs

**Description:**

This program reads the output of the `mgaps` program and adds alignment
information to it.  Part of the MUMmer4.x pipeline and can only be
used on the output of the `mgaps` program. This -D option alters this
behavior and only outputs the positions of difference, e.g. SNPs.

        USAGE:
        combineMUMs  [options]  <reference>  <query>  <mgapsfile>

        [options]    type 'combineMUMs -h' for a list of options.
        <reference>  the FastA reference file used in the comparison.
        <query>      the multi-FastA reference file used in the comparison.
        <mgapsfile>  the output of the 'mgaps' program run on the match
                     list produced by 'mummer' for the reference and query
                     files.

        OUTPUT:
        stdout           the 'mgaps' output interspersed with the alignments
                         of the gaps between adjacent MUMs.  An alignment of a
                         gap comes after the second MUM defining the gap, and
                         alignment errors are marked with a '^' character.  At
                         the end of each cluster is a summary line (keyword
                         "Region") noting the bounds of the cluster in the
                         reference and query sequences, the total number of
                         errors for the region, the length of the region and
                         the percent error of the region.
        witherrors.gaps  the 'mgaps' output with an appended column that lists
                         the number of alignment errors for each gap.


### delta-filter

**Description:**

This program filters a delta alignment file produced by either
nucmer or promer, leaving only the desired alignments which
are output to stdout in the same delta format as the
input. Its primary function is the LIS algorithm which
calculates the longest increasing subset of alignments. This
allows for the calculation of a global set of alignments
(i.e. 1-to-1 and mutually consistent order) with the -g option
or locally consistent with -1 or -m. Reference sequences can
be mapped to query sequences with -r, or queries to references
with -q. This allows the user to exclude chance and repeat
induced alignments, leaving only the "best" alignments between
the two data sets. Filtering can also be performed on length,
identity, and uniquenes.

        USAGE:
        delta-filter  [options]  <deltafile>

        [options]    type 'delta-filter -h' for a list of options.
        <deltafile>  the .delta output file from either nucmer or promer.

        OUTPUT:
        stdout  The same delta alignment format as output by nucmer and promer.

**Notes:**

For most cases the -m option is recommended, however -1 is
useful for applications that require a 1-to-1 mapping, such as
SNP finding. Use the -q option for mapping query contigs to
their best reference location.


### exact-tandems

**Description:**

This script finds exact tandem repeats in a specified FastA sequence
file.  It is a post-processor for `repeat-match` and provides a simple
interface and output for tandem repeat detection.

        USAGE:
        exact-tandems  <file>  <min match>

        <file>       the single sequence in FastA format to search for repeats.
        <min match>  the minimum match length for the tandems.

        OUTPUT:
        stdout  4 columns, the start of the tandem repeat, the total extent
                of the repeat region, the length of each repetitive unit, and
                to total copies of the repetitive unit involved.

### mgaps

**Description:**

This program reads a list of matches between a single-FastA reference
and a multi-FastA query file and outputs clusters of matches that lie
on similar diagonals and within a reasonable distance.  Part of the
MUMmer4.x pipeline and the output of `mummer` need not be processed
before passing it to this program, so long as `mummer` was run on a
1-vs-many or 1-vs-1 dataset.

        USAGE:
        mgaps  [options]  <  <matchlist>

        [options]    type 'mgaps -h' for a list of options.
        <matchlist>  A list of matches separated by their sequence FastA tags.
                     The columns of the match list should be start in
                     reference, start in query, and length of the match.

        OUTPUT:
        stdout  An ordered set of the input matches, separated by headers.
                Individual clusters are separated by a '#' character and
                sets of clusters from different sequences are separated by
                the FastA header tag for the query sequence.

**Notes:**

It is often very helpful to adjust the clustering parameters.  Check
`mgaps -h` for the list of parameters and check the source for a
better idea of how each parameter affects the result.  Often, it is
helpful to run this program a number of times with different
parameters until the desired result is achieved.


### mummer

**Description:**

This is the core program of the MUMmer package.  It is the suffix-tree
based match finding routine, and the main part of every MUMmer script.
For a detailed manual describing how to use this program, please refer
to "docs/maxmat3man.pdf" or in LaTeX format "docs/maxmat3man.tex". By
default, `mummer` now finds maximal matches regardless of their
uniqueness. Limiting the output to only unique matches can be specified
as a command line switch.

        USAGE:
        mummer  [options]  <reference>  <query> ...

        [options]    type 'mummer -help' for a list of options.
        <reference>  specifies the single or multi-FastA sequence file that
                     contains the reference sequence(s), to be aligned with
                     the queries.
        <query>      specifies the multi-FastA sequence file that contains
                     the query sequences, to be aligned with the references.
                     Multiple query files are allowed, up to 32.

        OUTPUT:
        stdout  a list of exact matches. Varies depending on input, refer to
                the manual specified in the description above.

**Notes:**

Many thanks to Stefan Kurtz for the latest mummer version. `mummer`
now behaves like the old `mummer2` program by default. The -mum switch
forces it to behave like `mummer1`, the -mumreference switch forces it
to behave like `mummer2` while the -maxmatch switch forces it to behave
like the old `max-match` program.


### mummerplot

**Description:**

`mummerplot` is a perl script that generates gnuplot scripts and data
collections for plotting with the gnuplot utility.  It can generate
2-d dotplots and 1-d coverage plots for the output of mummer, nucmer,
promer or show-tiling. It can also color dotplots with an identity
color gradient.

        USAGE:
        mummerplot  [options]  <matchfile>

        [options]    type 'mummerplot -h' for a list of options.
        <matchfile>  the output of 'mummer', 'nucmer', 'promer', or
                     'show-tiling'. 'mummerplot' will automatically determine
                     the format of the data it was given and produce the plot
                     accordingly.

        OUTPUT:
        out.gp     The gnuplot script, type 'gnuplot out.gp' to evaluate the
                   the gnuplot script.
        out.fplot
        out.rplot
        out.hplot  The forward, reverse and highlighted match information for
                   plotting with gnuplot.

        out.ps
        out.png    The plotted image file, postscript or png depending on the
                   selected terminal type.

**Notes:**

For alignments with multiple reference or query sequences, be sure to
use the -r -q or -R -Q options to avoid overlaying multiple plots in
the same space. For better looking color gradient plots, try the
postscript terminal and avoid the png terminal.


### nucmer2xfig

**Description:**

Script for plotting nucmer hits against a reference sequence. See top
of script for more information, or see if `mummerplot` or `mapview`
has the functionality required as they are properly maintained.


### repeat-match

**Description:**

Finds exact repeats within a single sequence.

        USAGE:
        repeat-match  [options]  <seq>

        [options]  type 'repeat-match -h' for a list of options.
        <seq>      the single sequence in FastA format to search for repeats.

        OUTPUT:
        stdout  3 columns, the start of the first copy of the repeat, the
                start of the second copy of the repeat, and the length of the
                repeat respectively.

**Notes:**

REPuter (freely available for universities) may be better suited for
most repeat matching, but `repeat-match` is open-source and has some
functionality that REPuter does not so we include it along with the
MUMmer package.


### show-aligns

**Description:**

This program parses the delta alignment output of nucmer and promer
and displays all of the pairwise alignments from the two sequences
specified on the command line.

        USAGE:
        show-aligns  [options]  <deltafile>  <IdR>  <IdQ>

        [options]    type 'show-aligns -h' for a list of options.
        <deltafile>  the .delta output file from either nucmer or promer.
        <IdR>        the FastA header tag of the desired reference sequence.
        <IdQ>        the FastA header tag of the desired query sequence.

        OUTPUT:
        stdout  each alignment header and footer describes the frame of the
                alignment in each sequence, and the start and finish
                (inclusive) of the alignment in each sequence.  At the
                beginning of each line of aligned sequence are two numbers, the
                top is the coordinate of the first reference base on that line
                and the bottom is the coordinate of the first query base on
                that line.  ALL coordinates reference the forward strand of the
                DNA sequence, even if it is a protein alignment.  A gap caused
                by an insertion or deletion is filled with a '.' character.
                Errors in a DNA alignment are marked with a '^' below the
                error.  Errors in an amino acid alignment are marked with a
                whitespace in the middle consensus line, while matches are
                marked with the consensus base and similarities are marked with
                a '+' in the consensus line.


### show-coords

**Description:**

This program parses the delta alignment output of nucmer and promer
and displays the coordinates, and other useful information about the
alignments.

        USAGE:
        show-coords  [options]  <deltafile>

        [options]    type 'show-coords -h' for a list of options.
        <deltafile>  the .delta output file from either nucmer or promer.

        OUTPUT:
        stdout  run 'show-coords' without the -H option to see the column
                header tags.  Here is a description of each tag.  Note that
                some of the below tags do not apply to nucmer data, and that
                all coordinates are inclusive and relative to the forward DNA
                strand.

        [S1]    Start of the alignment region in the reference sequence.

        [E1]    End of the alignment region in the reference sequence.

        [S2]    Start of the alignment region in the query sequence.

        [E2]    End of the alignment region in the query sequence.

        [LEN 1] Length of the alignment region in the reference sequence,
        measured in nucleotides.

        [LEN 2] Length of the alignment region in the query sequence, measured
        in nucleotides.

        [% IDY] Percent identity of the alignment, calculated as the
        (number of exact matches) / ([LEN 1] + insertions in the query).

        [% SIM] Percent similarity of the alignment, calculated like the above
        value, but counting positive BLOSUM matrix scores instead of exact
        matches.

        [% STP] Percent of stop codons of the alignment, calculated as
        (number of stop codons) / (([LEN 1] + insertions in the query) * 2).

        [LEN R] Length of the reference sequence.

        [LEN Q] Length of the query sequence.

        [COV R] Percent coverage of the alignment on the reference sequence,
        calculated as [LEN 1] / [LEN R].

        [COV Q] Percent coverage of the alignment on the query sequence,
        calculated as [LEN 2] / [LEN Q].

        [FRM]   Reading frame for the reference sequence and the reading frame
        for the query sequence respectively.  This is one of the columns
        absent from the nucmer data, however, match direction can easily be
        determined by the start and end coordinates.

        [TAGS]  The reference FastA ID and the query FastA ID.

                There is also an optional final column (turned on with the -w
        or -o option) that will contain some 'annotations'. The -o option will
        annotate alignments that represent overlaps between two sequences,
        while the -w option is antiquated and should no longer be used.
        Sometimes, nucmer or promer will extend adjacent clusters past one
        another, thus causing a somewhat redundant output, this option will
        notify users of such rare occurrences.

**Notes:**

The -c and -l options are useful when comparing two sets of assembly
contigs, in that these options help determine if an alignment spans an
entire contig, or is just a partial hit to a different read.  The -b
option is useful when the user wishes to identify sytenic regions
between two genomes, but is not particularly interested in the actual
alignment similarity or appearance.  This option also disregards match
orientation, so should not be used if this information is needed.


### show-diff

**Description:**

This program classifies alignment breakpoints for the
quantification of macroscopic differences between two
genomes. It takes a standard, unfiltered delta file as input,
determines the best mapping between the two sequence sets, and
reports on the breaks in that mapping.

        USAGE:
        show-diff  [options]  <deltafile>

        [options]    type 'show-diff -h' for a list of options.
        <deltafile>  the .delta output file from nucmer

        OUTPUT:
        stdout  Classified breakpoints are output one per line with
                the following types and column definitions. The first
                five columns of every row are seq ID, feature type,
                feature start, feature end, and feature length.

        Feature Columns

        IDR GAP gap-start gap-end gap-length-R gap-length-Q gap-diff
        IDR DUP dup-start dup-end dup-length
        IDR BRK gap-start gap-end gap-length
        IDR JMP gap-start gap-end gap-length
        IDR INV gap-start gap-end gap-length
        IDR SEQ gap-start gap-end gap-length prev-sequence next-sequence

        Feature Types

        [GAP] A gap between two mutually consistent ordered and
        oriented alignments. gap-length-R is the length of the
        alignment gap in the reference, gap-length-Q is the length of
        the alignment gap in the query, and gap-diff is the difference
        between the two gap lengths. If gap-diff is positive, sequence
        has been inserted in the reference. If gap-diff is negative,
        sequence has been deleted from the reference. If both
        gap-length-R and gap-length-Q are negative, the indel is
        tandem duplication copy difference.

        [DUP] A duplicated sequence in the reference that occurs more
        times in the reference than in the query. The coordinate
        columns specify the bounds and length of the
        duplication. These features are often bookended by BRK
        features if there is unique sequence bounding the duplication.

        [BRK] An insertion in the reference of unknown origin, that
        indicates no query sequence aligns to the sequence bounded by
        gap-start and gap-end. Often found around DUP elements or at
        the beginning or end of sequences.

        [JMP] A relocation event, where the consistent ordering of
        alignments is disrupted. The coordinate columns specify the
        breakpoints of the relocation in the reference, and the
        gap-length between them. A negative gap-length indicates the
        relocation occurred around a repetitive sequence, and a
        positive length indicates unique sequence between the
        alignments.

        [INV] The same as a relocation event, however both the
        ordering and orientation of the alignments is disrupted. Note
        that for JMP and INV, generally two features will be output,
        one for the beginning of the inverted region, and another for
        the end of the inverted region.

        [SEQ] A translocation event that requires jumping to a new
        query sequence in order to continue aligning to the
        reference. If each input sequence is a chromosome, these
        features correspond to inter-chromosomal translocations.

**Notes:**

The estimated number of features, take inversions for example,
represents the number of breakpoints classified as bordering
an inversion. Therefore, since there will be a breakpoint at
both the beginning and the end of an inversion, the feature
counts are roughly double the number of inversion events. In
addition, all counts are estimates and do not represent the
exact number of each evolutionary event.

Summing the fifth column (ignoring negative values) yeilds an
estimate of the total inserted sequence in the
reference. Summing the fifth column after removing DUP
features yields an estimate of the total amount of unique
(unaligned) sequence in the reference. Note that unaligned
sequences are not counted, and could represent additional
"unique" sequences. Use the `dnadiff` script if you must
recover this information. Finally, the -q option switches
references for queries, and uses the query coordinates for the
analysis.


### show-snps

**Description:**

This program reports polymorphism contained in a delta encoded
alignment file output by either nucmer or promer. It catalogs
all of the single nucleotide polymorphisms (SNPs) and
insertions/deletions within the delta file
alignments. Polymorphisms are reported one per line, in a
delimited fashion similar to show-coords. Pairing this program
with the appropriate MUMmer tools can create an easy to use
SNP pipeline for the rapid identification of putative SNPs
between any two sequence sets.

        USAGE:
        show-snps  [options]  <deltafile>

        [options]    type 'show-snps -h' for a list of options.
        <deltafile>  the .delta output file from either nucmer or promer.

        OUTPUT:
        stdout  Standard output has column headers with the following
                meanings. Not all columns will be output by default,
                see 'show-snps -h' for switch to control the output.

        [P1]    SNP position in the reference.

        [SUB]   Character in the reference.

        [SUB]   Character in the query.

        [P2]    SNP position in the query.

        [BUFF]  Distance from this SNP to the nearest mismatch (end of
        alignment, indel, SNP, etc) in the same alignment.

        [DIST]  Distance from this SNP to the nearest sequence end.

        [R]     Number of repeat alignments which cover this reference
        position, >0 means repetitive sequence.

        [Q]     Number of repeat alignments which cover this query
        position, >0 means repetitive sequence.

        [LEN R] Length of the reference sequence.

        [LEN Q] Length of the query sequence.

        [CTX R] Surrounding context sequence in the reference.

        [CTX Q] Surrounding context sequence in the query.

        [FRM]   Reading frame for the reference sequence and the
        reading frame for the query sequence respectively. Simply
        'forward' 1, or 'reverse' -1 for nucmer data.

        [TAGS]  The reference FastA ID and the query FastA ID.

**Notes:**

It is often helpful to run this with the -C option to assure
reported SNPs are only reported from uniquely aligned regions.


### show-tiling

**Description:**

This program attempts to construct a tiling path out of the query
contigs as mapped to the reference sequences.  Given the delta
alignment information of a few long reference sequences and many small
query contigs, `show-tiling` will determine the best location on a
reference for each contig.  Note that each contig may only be tiled
once, so repetitive regions may cause this program some difficulty.
This program is useful for aiding in the scaffolding and closure of an
unfinished set of contigs, if a suitable, high similarity, reference
genome is available.  Or, if using promer, `show-tiling` will help
in the identification of syntenic regions and their contig's mapping
the the references.

        USAGE:
        show-tiling  [options]  <deltafile>

        [options]    type 'show-tiling -h' for a list of options.
        <deltafile>  the .delta output file from either nucmer or promer.

        OUTPUT:
        stdout  Standard output has 8 columns: start in reference, end in
                reference, gap between this contig and the next, length of this
                contig, alignment coverage of this contig, average percent
                identity of the alignments for this contig, orientation of this
                contig, contig ID. All matches to a reference are headed by the
                FASTA tag of that reference.  Output with the -a option is the
                same as 'show-coords -cl' when run on nucmer data.

**Notes:**

When run with the -x option, `show-tiling` will produce an XML output
format that can be accepted by TIGR's open source scaffolding software
'Bambus' as contig linking information.

## Obsolete programs

The programs `mapview`, `run-mummer1`, `run-mummer3` and `nucmer2xfig`
are now obsolete. The original documentation is still available in
[OBSOLETE.md](OBSOLETE.md).

## CONTACT INFORMATION

Please address questions and bug reports via the [github issue tracker](https://github.com/gmarcais/mummer/issues).
