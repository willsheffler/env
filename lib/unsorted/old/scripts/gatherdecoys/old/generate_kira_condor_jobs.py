import sys,os

template = """Executable 	= /users/sheffler/rosetta/rosetta.gcc
Universe 	= vanilla
Log		= %(prot)s_relax.log
Output		= %(prot)s_relax.out.$(Process)
Error		= %(prot)s_relax.err.$(Process)

Arguments	=  rn %(prot)s _ -relax -farlx -ex1 -ex2 -use_trie -s rosetta_data/%(prot)s_0001.pdb -nstruct 200
Queue 5
"""

kiraprots = ('1aa3','1acp','1ail','1b72','1bf4','1erv','1ig5',
             '1pgx','1r69','1tif','1tig','1utg','1vii','256b')
decoydir = '/users/sheffler/data/decoys/kira/'

for prot in kiraprots:
    os.chdir(decoydir+prot)
    condorfile = open(prot+'_relax.condor','w')
    condorfile.write(template%{'prot':prot})
    condorfile.close()
