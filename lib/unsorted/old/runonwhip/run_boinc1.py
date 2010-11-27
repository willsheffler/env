import os,sys,time

if len(sys.argv) < 2:
    feature = "atomdist"
else:
    feature = sys.argv[1]

date = time.strftime('%x').replace('/','-')

outfiles = os.popen("/bin/ls /users/boinc/results/BARCODE_/BARCODE_*299_0.out").read().strip().split()
outnames = [ o.replace('/','-')[1:] for o in outfiles ]

print outfiles
print outnames

workdir = "/work/sheffler/project/rosetta_tmp"
binary  = "/users/sheffler/svn/rosetta_fix_decoyfeatures/rosetta.gcc"
options = "-pdbstats distances -fa_input -try_both_his_tautomers -silent_input -l %(outfile)s > %(outname)s-pdbstats-distances.log"
whips   = [8,9,10,11]

cmd = "ssh -n whip%(whip)02i 'cd %(workdir)s ; nice %(binary)s "+options+"' &"

outfile = outfiles[0]
outname = outnames[0]
whip = whips[0]
print cmd%vars()

#for ii in range(10):
#    whip = ii+1
    #print vars()
    #print pdbcmd%vars()
    #print ds2cmd%vars()
#    os.system(pdbcmd%vars())
#    os.system(ds2cmd%vars())



