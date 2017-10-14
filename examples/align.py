#! /usr/bin/env python

import mummer
import sys

fd = open(sys.argv[1])
ref = fd.read()
fd.close()

fd = open(sys.argv[2])
qry = fd.read()
fd.close()

o = mummer.Options()
o.minmatch(10)
o.mincluster(10)
aligns = mummer.align_sequences(ref, qry, o)

for a in aligns:
    print("%d %d %d %d %d %d %d" % (a.sA, a.eA, a.sB, a.eB, a.Errors, a.SimErrors, a.NonAlphas))
    for d in a.delta:
        print(d)
    print("0")
