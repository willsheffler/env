import sys,os,PathsFile

def setupJump(prot,chain,native,pfile,nstruct=1,queue=1,paths='paths.txt'):

    template = """Executable      = /users/sheffler/rosetta/rosetta.gcc
Universe        = vanilla
Log             = jump.%(series)s%(prot)s%(chain)s.log
Output          = jump.%(series)s%(prot)s%(chain)s.out.$(Process)
Error           = jump.%(series)s%(prot)s%(chain)s.err.$(Process)

Arguments       =  %(series)s %(prot)s %(chain)s -paths %(paths)s -jumping -n %(native)s -pairing_file %(pairfile)s -sheet1 2 -nstruct %(nstruct)i -seed_offset $(Process)
Queue %(queue)i
"""

    vanilla = """Executable      = /users/sheffler/rosetta/rosetta.gcc
Universe        = vanilla
Log             = jump.%(series)s%(prot)s%(chain)s.log
Output          = jump.%(series)s%(prot)s%(chain)s.out.$(Process)
Error           = jump.%(series)s%(prot)s%(chain)s.err.$(Process)

Arguments       =  %(series)s %(prot)s %(chain)s -paths %(paths)s -nstruct %(nstruct)i -seed_offset $(Process)
Queue %(queue)i
"""
    
#     prot    = sys.argv[1]
#     chain   = sys.argv[2]
#     native  = sys.argv[3]
#     pfile   = sys.argv[4]
#     nstruct = int(sys.argv[5])
    
    pairs = [x.split() for x in open(pfile).readlines() if len(x.split())==4]
    
    for ii in range(len(pairs)):                        
        if ii > 9:
            print "can only handle up to 10 pairings!"
            break
        pairfile = "%s%s_pairs_%i.pair"%(prot,chain,ii)
        f = open(pairfile,'w')
        f.write("%s %s %s %s\n"%(pairs[ii][0],pairs[ii][1],pairs[ii][2],pairs[ii][3]))            
        f.close()
        condorfile = "%s%s_pairs_%i.condor"%(prot,chain,ii)
        f = open(condorfile,'w')
        series = "p"+`ii`
        f.write(template%vars())
        f.close()        
        print "~/rosetta/rosetta_051029/rosetta.gcc %(series)s %(prot)s %(chain)s -trajectory -jumping -n %(native)s -pairing_file %(pairfile)s -sheet1 2 -nstruct 1 > test.%(series)s%(prot)s%(chain)s.log 2> test.%(series)s%(prot)s%(chain)s.err ; "%vars()

    condorfile = "%s%s_allpairs.condor"%(prot,chain)
    f = open(condorfile,'w')
    series = "ap"
    f.write(template%vars())
    f.close()        

    condorfile = "%s%s_vanilla.condor"%(prot,chain)
    f = open(condorfile,'w')
    series = "vn"
    f.write(vanilla%vars())
    f.close()        


prots = "1acf_ 1fkb_ 1fna_ 1gvp_ 1kpeA 1npsA 1tul_ 1urnA 1who_ 2acy_".split()
for p in prots:
    prot  = p[:4]
    chain = p[4]
    nat   = p+'/'+p[:4]+'.pdb'
    pfile = p[:4]+'.pairings'
    nstruct = 1000
    queue = 1
    d = './'+p
    os.system("mkdir -p "+'./'+p+"_traj")
    pargs = {'pdb1':d,'fragments':d,'structure':d,'sequence':d,'movie':'./'+p+"_traj"}
    paths = PathsFile.makePathsFile(pargs,filename='paths.%s.txt'%p)
    setupJump(prot,chain,nat,pfile,nstruct,queue,paths)
