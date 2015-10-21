#!/bin/sh
#
# Find exact tandem repeats in specified file involving an
# exact duplicate of at least the specified length

filename=$1
matchlen=$2

bindir="@BIN_DIR@"
libdir="@LIB_DIR@"

if [ -z "$filename" -o -z "$matchlen" ]; then
   echo "USAGE:  $0 <file> <min-match-len>"
   exit 1
fi

#echo "Finding matches and tandem repeats"
# Trick: pipe exit status of first command to exit at end of the pipeline
exec 4>&1; { \
    { $bindir/repeat-match -t -n $matchlen $filename; echo $? >&3; } | \
        tail -n +3 | sort -k1n -k2n | awk -f $libdir/tandem-repeat.awk >&4; \
} 3>&1 | exit `cat`

