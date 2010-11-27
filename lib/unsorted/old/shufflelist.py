import sys,os,random

l = open(sys.argv[1]).readlines()
random.shuffle(l)
for x in l:
    print x.strip()
