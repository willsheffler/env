import sys,os,re

native   = sys.argv[1]
chain    = sys.argv[2]
outfiles = sys.argv[3:]

os.system("mkdir -p results")

for outfile in outfiles:
    if out[-4:] != ".out":
        print "out files should end with .out!"
        sys.exit()
    name = outfile[:-4]
    cmd = []
    cmd.append("cp %(name)s.out results/"%vars())
    cmd.append("grep SCORE %(name)s.out > results/%(name)s.score"%vars())
    cmd.append("grep FOLD_TREE %(name)s.out > results/%(name)s.foldtree"%vars())
    cmd.append("egrep -v 'FOLD_TREE|JUMP' %(name)s.out > results/%(name)s.cleanout"%vars())
    cmd.append("~pbradley/python/make_coords_file.py %(native)s %(chain)s %(name).cleanout > %(name)s.coords"%vars())
    cmd.append("mkdir -p %(name)s_cluster"%vars())
    cmd.append("~pbradley/C/cluster_info_silent.out %(name)s.cleanout %(name)s.coords %(name)s_cluster/%(name)s 5,25,75 3,4"%vars())
    cmd.append("~pbradley/python/make_color_trees.py %(name)s_cluster/%(name)s 1 25 -cartoons 12"%vars())
    cmd.append("~pbradley/tools/make_new_plot.py %(name)s_cluster/%(name)s.cluster.contacts"%vars())
               
    print cmd
    os.system(cmd)
