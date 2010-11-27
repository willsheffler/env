#!/usr/bin/python

import string
from sys import stderr,argv
from os import popen

def Help():
    print 'uggh'
    assert 0==1


if len(argv) not in [2,3]:
    Help()
elif len(argv) == 2:
    N = 20
else:
    N = int(argv[2])

prefix = argv[1]

lines = map(string.split,popen('grep "CLUSTER_SCORE" %s.info'%prefix).readlines())

for line in lines[:N]:
    print 'load %s.%s.pdb\necho "%s"\n'%(prefix,line[-1],line[-1])
    print 'backbone\ncolor group\nselect */1\nbackbone 200\npause'
    print 'select cys\nspacefill 500\nrestrict */1\npause'
    print 'select hydrophobic\ncolor blue\nselect polar\ncolor green\nselect charged\ncolor red\nselect gly\ncolor gold\nselect cys and */1\ncolor purple\npause'
    print 'select all\ncolor temperature\necho "colored by conservation"\npause'
    print 'zap'
    
