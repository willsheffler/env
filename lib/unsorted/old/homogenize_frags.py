#!/usr/bin/python
from sys import argv,exit
from  os import system as run,listdir


if len(argv) > 1:
    sourcedir = argv[1]
else:
    print "usage homogenize_frags.py <frag dir> [<target dir>]"
    exit()
    
if len(argv) > 2:
    destdir   = argv[2]
else:
    destdir = sourcedir

F = listdir(sourcedir)

P = [f[:4] for f in F if f.endswith('.pdb') and len(f) is 8]

for f in F:
    if f[:4] in P and f[5] in ('.','-','_') and f[4] != '_':
        newname = f[:4]+'_'+f[5:]
        if newname not in F:
            cmd = "rsync -z %s/%s %s/%s"%(sourcedir,f,destdir,newname)
        print cmd
        run(cmd)
    if f[2:6] in P and f[7:9] in ('03','09') and f[6] != '_':
        newname = "aa%s_"+f[7:9]+"_99.999_will"
        if newname not in F:
            cmd = 'rsync -z %s/%s %s/%s'%(sourcedir,f,destdir,newname)
        print cmd
        run(cmd)
    if f in [p+'.pdb' for p in P]:
        cmd = "rsync -z %s/%s %s/%s"%(sourcedir,f,destdir,f)
        exe = "/users/sheffler/scripts/setpdbchain.py"
        cmd = "%s '*' ' ' %s/%s"%(exe,destdir,f)
        print cmd
        run(cmd)


