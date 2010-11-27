#!/usr/bin/python
import sys,os

if len(sys.argv) != 4:
    print "usage: copypdbs.py <file list> <source dir> <dest dir>"
    sys.exit()

list   = open(sys.argv[1])
source = sys.argv[2]
dest   = sys.argv[3]

if source[-1] is not '/': source = source+'/'
if dest[-1]   is not '/': dest = dest+'/'

for fname in list.readlines():
    if fname[:8] == 'filename':
        continue
    fname = fname.split()[0]
    fname = fname.split('.')[0]+'.pdb'
    #print "cp "+source+'/'+fname+" "+dest+'/'
    os.system("cp "+source+fname+" "+dest)

list.close()
