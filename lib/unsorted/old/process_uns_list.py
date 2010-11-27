#!/usr/bin/python
# From decoystats output first grep "DS " <file> to get just the
# UNS entries. This is the first input file.
# Second file should be the total_residue for each pdb in the 
# decoystats run with one number per line
# output is a verbose list of the number of burried unsatisfied
# hbonds for every residue (most are zero)

import sys

if len(sys.argv) != 3:
    print "usage: process_uns_hbond.py <UNS list> <pdb lengths>"
    sys.exit()

pdbnresfile = open(sys.argv[2])
pdbnres = [int(x) for x in pdbnresfile]

unsfile = open(sys.argv[1])

unsid = []
unsres = []
unsaa = []
npdbs = 0
for line in unsfile.xreadlines():
    x = line.split()
#    print x
    if x[1] == "START_FILE":
        npdbs += 1
        unsid.append(x[2])
        unsres.append([])
        unsaa.append([])
    else:
        unsres[npdbs-1].append(int(x[4]))
        unsaa[npdbs-1].append(x[3])

#print unsid
#print unsres
#print unsaa

for i in range(npdbs):
    k = 0
    for j in range(pdbnres[i]):
        count = 0
        aa = '_'
        while k < len(unsres[i]) and unsres[i][k] == j+1:
            count += 1
            aa = unsaa[i][k]
            k += 1
            #if k >= len(unsres[i]):
            #    break
        print "%i\t%s\t%i\t%i\t%s"%(i+1,unsid[i],j+1,count,aa)
