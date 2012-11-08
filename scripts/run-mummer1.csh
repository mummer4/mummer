#!__CSH_PATH -f
#
# **SEVERELY** antiquated script for running the mummer 1 suite
# -r option reverse complements the query sequence, coordinates of the reverse
# matches will be relative to the reversed sequence
#

set ref = $1
set qry = $2
set pfx = $3
set rev = $4

set bindir = __BIN_DIR

if($ref == '' || $qry == '' || $pfx == '') then
    echo "USAGE: $0 <fasta reference> <fasta query> <prefix> [-r]"
    exit(-1)
endif

echo "Find MUMs"
$bindir/mummer -mum -l 20 $rev $ref $qry | tail +2 > $pfx.out
echo "Determine gaps"
$bindir/gaps $ref $rev < $pfx.out > $pfx.gaps
echo "Align gaps"
$bindir/annotate $pfx.gaps $qry > $pfx.align
mv witherrors.gaps $pfx.errorsgaps
