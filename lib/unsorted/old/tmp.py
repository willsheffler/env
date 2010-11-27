import sys,os

actual = os.listdir("/work/sheffler/data/pdb/rosetta_pdb_chains")
newlist = [x.strip() for x in open("/work/sheffler/pdb_energies/pdb_small200.list").readlines()]


for n in newlist:
    if n+".gz" in actual:
        pass #print n
    else:
        print n
