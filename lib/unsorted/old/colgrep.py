#!/usr/bin/python
import os,sys


if(len(sys.argv)!=3):
    print "usage: colgrep pattern file"
    sys.exit()
pattern  = sys.argv[1]
fname    = sys.argv[2]
restrict = pattern       # take out lines with this char

f = open(fname)
header = f.readline().split()

f.close()
cols = [str(i) for i in range(len(header)) if header[i].find(pattern)!=-1 ]

feilds = ' '.join(cols)

os.system( "echo rowname " + ' '.join([header[int(i)] for i in cols])  )
os.system( "grep -v %(restrict)s %(fname)s | cut -d ' ' -f '1 %(feilds)s'"%vars() )

