import os,sys,time

if len(sys.argv) < 2:
    feature = "atomdist"
else:
    feature = sys.argv[1]

date = time.strftime('%x').replace('/','-')

pdbcmd = "ssh -n whip%(whip)02i 'cd pdb_energies; nice ~/rosetta/rosetta_features/rosetta.gcc -score -decoyfeatures %(feature)s -fa_input -try_both_his_tautomers -l pdb.list%(ii)i > /users/sheffler/data/raw/pdb_%(feature)s_%(date)s.log%(ii)i 2> /users/sheffler/data/raw/pdb_%(feature)s_%(date)s.err%(ii)i' &"
ds2cmd = "ssh -n whip%(whip)02i 'cd decoy_energy_polar; nice ~/rosetta/rosetta_features/rosetta.gcc -score -decoyfeatures %(feature)s -fa_input -try_both_his_tautomers -l ds2_cull.list%(ii)i > /users/sheffler/data/raw/ds2_%(feature)s_%(date)s.log%(ii)i 2> /users/sheffler/data/raw/ds2_%(feature)s_%(date)s.err%(ii)i' &"

for ii in range(10):
    whip = ii+1
    #print vars()
    #print pdbcmd%vars()
    #print ds2cmd%vars()
    os.system(pdbcmd%vars())
    os.system(ds2cmd%vars())



