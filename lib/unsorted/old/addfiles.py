import sys,os
from Numeric import array

if len(sys.argv) == 1:
    print "at least one arg!"

try:
    skip  = int(sys.argv[1])
    fnames = sys.argv[2:]
except:
    skip  = 0
    fnames = sys.argv[1:]

if len(fnames) < 2:
    print "dont' call this unless you have at least two files!"


files = [open(f) for f in fnames]
for ii in range(skip):
    for f in files:
        f.readline()
        
while True:
    try:
        n = array([int(x) for x in files[0].readline().split()])
        if len(n) < 100:
            break
        for file in files[1:]:
            n += array([int(x) for x in file.readline().split()])
        print " ".join(str(n)[1:-1].split())
    except:
        break

        
