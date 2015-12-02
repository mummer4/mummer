#! /bin/sh

set -e

NCPUS=$(grep -c '^processor' /proc/cpuinfo)

autoreconf -i
./configure YAGGO=~/bin/yaggo
make -j $NCPUS
make -j $NCPUS check || { cat test-suite.log; false; }
