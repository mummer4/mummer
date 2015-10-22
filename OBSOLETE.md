##Obsolete software

The following programs were present in earlier version of MUMmer, but
have been superseeded by other programs. They are not very useful
anymore, deemed obsolete, and not maintained anymore. They are not
built and would have to be built by hand.

###run-mummer1

**Description:**

This script is taken directly from MUMmer1.0 and is best used to
align two sequences in which there is high similarity and no re-
arrangements.  Common use cases are: aligning two finished bacterial
chromosomes.  Please refer to "docs/run-mummer1.README" for the
original documentation for this script and its output.

        USAGE:
        run-mummer1  <seq1>  <seq2>  <tag>  [-r]

        <seq1>  specifies the file with the first sequence in FastA format.
                No more than one sequence is allowed.
        <seq2>  specifies the file with the second sequence in FastA format.
                No more than one sequence is allowed.
        <tag>   specifies the prefix to be used for the output files.
        [-r]    is an optional parameter that will reverse complement the
                second sequence.

        OUTPUT:
        out.align       the out.gaps file interspersed with the alignments
                        of the gaps.
        out.errorsgaps  the out.gaps file with an extra column stating the
                        number of errors contained in each gap.
        out.gaps        an ordered (clustered) list of matches with position
                        information, and gap distances between each match.
        out.out         a list of all maximal unique matches between the two
                        input sequences ordered by their start position in the
                        second sequence.

**Notes:**

All output coordinates reference their respective strand.  This means
that if the -r switch is active, coordinates that reference the
second sequence will be relative to the reverse complement of the
second sequence.  Please use nucmer or promer if this coordinate
system is confusing.

Eventually, this script's components will be rewritten to work
with the new MUMmer format standards and phased out in favor of the
new components and wrapping script.


###run-mummer3

**Description:**

This script is the improved version of the MUMmer1.0 run-mummer1
script.  It uses a new clustering algorithm that appropriately
handles multiple sequence rearrangements and inversions.  Because
of this, it can handle more divergent sequences better than
run-mummer1.  In addition, it allows a multi-FastA query file for
1-vs-many sequence comparisons.  Please refer to
"docs/run-mummer3.README" for more detailed documentation of this
script and its output.

        USAGE:
        run-mummer3  <reference>  <query>  <prefix>

        <reference>  specifies the file with the reference sequence in FastA
                     format.  No more than one sequence is allowed.
        <query>      specifies the multi-FastA sequence file that contains
                     the query sequences.
        <prefix>     specifies the file prefix for the output files.

        OUTPUT:
        out.align       the out.gaps file interspersed with the alignments
                        of the gaps.
        out.errorsgaps  the out.gaps file with an extra column stating the
                        number of errors contained in each gap.
        out.gaps        an ordered (clustered) list of matches with position
                        information, and gap distances between each match.
        out.out         a list of all maximal unique matches between the two
                        input sequences ordered by their start position in the
                        second sequence.

###gaps

**Description:**

This program reads a list of unique matches between two strings and
outputs the longest consistent set of matches, followed by all the
other matches.  Part of the MUMmer1.0 pipeline and the output of the
`mummer` program needs to be processed (to strip all non-match lines)
before it can be passed to this program.

        USAGE:
        gaps  <seq1>  [-r]  <  <matchlist>

        <seq1>       The first sequence file that the match list represents.
        <matchlist>  A simple list of matches and NO header lines or other
                     mumbo jumbo.  The columns of the match list should be
                     start in the reference, start in the query, and length
                     of the match.
        [-r]         Simply puts the string "reverse" on the header of the
                     output so 'annotate' knows to reverse the second
                     sequence.

        OUTPUT:
        stdout  an ordered set of the input matches, separated by headers.
                The first set is the longest consistent set of matches and
                the second set is all other matches.

**Notes:**

This program will eventually be rewritten to be interchangeable with
`mgaps`, so that it may be plugged into the nucmer or promer
pipelines.

###mapview

**Description:**

`mapview` is a utility program for displaying sequence alignments as
provided by MUMmer, nucmer or promer. This program takes the output
from these alignment routines and converts it to a FIG, PDF or PS
file for visual analysis. It can also break the output into multiple
files for easier viewing and printing. Please refer to
"docs/mapview.README" for a more detailed description and explination.

        USAGE:
        mapview  [options]  <coords file>  [UTR coords]  [CDS coords]

        [options]       type 'mapview -h' for a list of options.
        <coords file>   show-coords output file
        [UTR coords]    UTR coordinate file in GFF format
        [CDS coords]    CDS coordinate file in GFF format

        OUTPUT:
        Default output format is an xfig file, however this can be changed to
        a postscript of PDF file with the -f option. See 'mapview -h' for a
        list of available formatting options.

**Notes:**

The produce the coords file input, `show-coords` must be run with the
-r -l options. To reduce redundant matches in promer output, run
show-coords with the -k option. To generate output formats other than
xfig, the fig2dev utility must be available from the system path. For
very large reference genomes, FIG format may be the only option that
will allow the entire display to be stored in one file, as fig2dev has
problems if the output is too large.
