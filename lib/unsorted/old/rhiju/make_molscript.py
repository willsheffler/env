#!/usr/bin/python
## make molscript images for list of proteins
## images are put in <pdb-file>.molscript.ps
##

## From Phil BRadley -- modified by Rhiju, Sep 2005.


from phil import *
from rasmol import Rclean
from sys import stderr

def Help():
    print 'Usage:  pdb name  (or *.pdb)'
    print ''
    print 'Runs rasmol and then makes a postscript of protein'
    print 'in desired view. From Phil Bradley.'
    print ''
    assert 0==1


## def Get_scop_class( id ):
##     lines = popen('grep %s /data/pbradley/scop/dir.cla.scop.txt_1.65'%id).readlines()
##     if len(lines) != 1:
##         return '<NOSCOP>'
##     else:
##         return string.split(lines[0])[3]


if len(argv) not in [2,3]:
    Help()

rout,rin = popen2('/users/pbradley/rasmol_32BIT')

files = argv[1:]


for file in files:

    system('rm junk.mol junk2.mol')

    ## will appear in the middle of the picture
    label = string.split(file,'/')[-1][:-4]
    ps_file = file+'.molscript.ps'

    ## get good view by hand:
    rin.write('zap\nload %s\ncolor group\nset background white\nwireframe off\ncartoons\n'%file)
    rin.flush()
    Rclean(rin,rout)
    raw_input()

    rin.write('write molscript "junk.mol"\n')
    rin.flush()
    Rclean(rin,rout)

    lines = open('junk.mol','r').readlines()
    for line in lines:
        print id,line[:-1]
    out = open('junk2.mol','w')
    done = 0
    for line in lines:
        if line[:4] in ['  co', '  he','  st','  tu'] and not done:
            out.write('  set residuecolour amino-acids rainbow, colourparts on;\n')
            out.write('  set labeloffset 0 -20 0;\n')
            out.write('  set labelsize 50;\n')
            out.write('  label 0 0 0 "%s";\n'%label)
            done = 1
        out.write(line)
    out.close()

    system('/users/rohl/bin/molscript -ps < junk2.mol > %s'%ps_file)


rin.write('quit\n')
rin.flush()


