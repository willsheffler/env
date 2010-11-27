import sys,os

list = open(sys.argv[1])

for line in list.readlines():
    pdb = line[:-1]
    dssp = line[:-5]+'.dssp'
    err  = dssp+'.err'
    print pdb
    cmd = "~pbradley/dssp "+pdb+" > "+dssp+" 2> "+err
    os.system(cmd)
