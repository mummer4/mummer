#!/bin/sh
#
# **SEVERELY** antiquated script for running the mummer 1 suite
# -r option reverse complements the query sequence, coordinates of the reverse
# matches will be relative to the reversed sequence
#

ref=$1
qry=$2
pfx=$3
rev=$4

bindir="@BIN_DIR@"

if [ -z "$ref" -o -z "$qry" -o -z "$pfx" ]; then
    echo "USAGE: $0 <fasta reference> <fasta query> <prefix> [-r]"
    exit 1
fi

echo >&2 "*** The run-mummer1 script is deprecated ***"
echo "Find MUMs"
$bindir/mummer -mum -l 20 $rev $ref $qry | tail -n +2 > $pfx.out
echo "Determine gaps"
$bindir/gaps $ref $rev < $pfx.out > $pfx.gaps
echo "Align gaps"
$bindir/annotate $pfx.gaps $qry > $pfx.align
mv witherrors.gaps $pfx.errorsgaps
