import sys,os

template = """Executable 	= /users/sheffler/rosetta/rosetta.gcc
Universe 	= vanilla
Log		= %(prot)s_relax.log
Output		= %(prot)s_relax.out.$(Process)
Error		= %(prot)s_relax.err.$(Process)

Arguments	=  rn %(prot)s _ -relax -farlx -ex1 -ex2 -use_trie -s rosetta_data/%(prot)s_0001.pdb -nstruct 200
Queue 5
"""

bqianprots = ('t196','t199','t205','t206','t223','t224','t232','t234','t249','t262','t279')
decoydir = '/data/sheffler/decoys/bqian/'

for prot in bqianprots:
    os.chdir(decoydir+prot)
    condorfile = open(prot+'_relax.condor','w')
    condorfile.write(template%{'prot':prot})
    condorfile.close()
