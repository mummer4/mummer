#!__CSH_PATH -f
#
# Find exact tandem repeats in specified file involving an
# exact duplicate of at least the specified length

set filename = $1
set matchlen = $2

set bindir = __BIN_DIR
set scriptdir = __SCRIPT_DIR

if  ($filename == '' || $matchlen == '')  then
    echo "USAGE:  $0 <file> <min-match-len>"
    exit -1
endif

echo "Finding matches"
$bindir/repeat-match -t -n $matchlen $filename | tail +3 > $$.tmp.matches
if  ($status != 0)  exit -1

echo "Tandem repeats"
sort -k1n -k2n $$.tmp.matches | awk -f $scriptdir/tandem-repeat.awk
rm -f $$.tmp.matches
