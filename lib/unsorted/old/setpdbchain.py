#!/usr/bin/python
import sys,os

if len(sys.argv) < 4:
    print "usage: setpdbchain.py <pdb files> <chainfrom> <chainto>"

chainfrom   = sys.argv[1]
chainto     = sys.argv[2]

pdbcount = 0
for pdbfile in sys.argv[3:]:
    pdb = open(pdbfile)
    lines = pdb.readlines()
    pdb.close()
    pdb = open(pdbfile,'w')
    count = 0
    for line in lines:
        if line[:5] == "ATOM " and \
           (line[20:23] == " %s "%chainfrom or chainfrom=='*'):
            pdb.write(line[:21] + chainto + line[22:])
            count += 1
        else:
            pdb.write(line)
    pdb.close()
    print "replaced %i chains in %s"%(count,pdbfile)
    pdbcount += 1
print "replaced chain in %i pdbs"%pdbcount

