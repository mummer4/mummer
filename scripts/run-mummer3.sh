#!/bin/sh
#
# for running the basic mummer 3 suite, should use nucmer instead when possible
# to avoid the confusing reverse coordinate system of the raw programs.
#
# NOTE:  be warned that all reverse matches will then
# be relative to the reverse complement of the query sequence.
#
# Edit this script as necessary to alter the matching and clustering values
#

ref=$1
qry=$2
pfx=$3

bindir="@BIN_DIR@"
libexecdir="@LIBEXEC_DIR@"

if [ -z "$ref" -o -z "$qry" -o -z "$pfx" ]; then
    echo "USAGE: $0 <fasta reference> <multi-fasta query> <prefix>"
    exit 1
fi

echo "Find MUMs"
$bindir/mummer -mumreference -b -l 20 $ref $qry > $pfx.out
echo "Determine gaps"
$libexecdir/mgaps -l 100 -f .12 -s 600 < $pfx.out > $pfx.gaps
echo "Align gaps"
$bindir/combineMUMs -x -e .10 -W $pfx.errorsgaps $ref $qry $pfx.gaps > $pfx.align
