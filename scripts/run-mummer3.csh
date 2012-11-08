#!__CSH_PATH -f
#
# for running the basic mummer 3 suite, should use nucmer instead when possible
# to avoid the confusing reverse coordinate system of the raw programs.
#
# NOTE:  be warned that all reverse matches will then
# be relative to the reverse complement of the query sequence.
#
# Edit this script as necessary to alter the matching and clustering values
#

set ref = $1
set qry = $2
set pfx = $3

set bindir = __BIN_DIR

if($ref == '' || $qry == '' || $pfx == '') then
    echo "USAGE: $0 <fasta reference> <multi-fasta query> <prefix>"
    exit(-1)
endif

echo "Find MUMs"
$bindir/mummer -mumreference -b -l 20 $ref $qry > $pfx.out
echo "Determine gaps"
$bindir/mgaps -l 100 -f .12 -s 600 < $pfx.out > $pfx.gaps
echo "Align gaps"
$bindir/combineMUMs -x -e .10 -W $pfx.errorsgaps $ref $qry $pfx.gaps > $pfx.align
