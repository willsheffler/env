import sys,os

nstruct = int(sys.argv[1])
pdbs = sys.argv[2:]

for ii in range(1,nstruct+1):
    if ii >= len(pdbs):
        break
    os.system("ln -s "+pdbs[ii]+" decoy_"+`ii`+".pdb")


