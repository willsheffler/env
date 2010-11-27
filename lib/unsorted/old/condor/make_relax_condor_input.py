#!/usr/bin/python
import sys,os

if len(sys.argv) < 4 or len(sys.argv) > 7:
    print "usage make_relax_condor_input.py <series> <prot> <chain> [ <idealized start> <nstruct> <nproc> ]"
    sys.exit()

args = {}
args['series']  = sys.argv[1]
args['prot']    = sys.argv[2]
args['chain']   = sys.argv[3]
args['start']   = args['prot'] + "_0001.pdb"
args['nstruct'] = 200
args['queue']   = 5
if len(sys.argv) >= 5: args['start']   = sys.argv[4]
if len(sys.argv) >= 6: args['nstruct'] = int(sys.argv[5])
if len(sys.argv) == 7: args['queue']   = int(sys.argv[6])

template = """Executable      = /users/sheffler/rosetta/rosetta.gcc
Universe        = vanilla
Log             = relax_native.log
Output          = relax_native.out.$(Process)
Error           = relax_native.err.$(Process)

Arguments       =  %(series)s %(prot)s %(chain)s -relax -farlx -ex1 -ex2 -use_trie -s %(start)s -nstruct %(nstruct)i
Queue %(queue)i
"""

print template%args
