import pymol
from pymol import cmd
import sys,os,random


def useOccRadii(sel="all"):
    for a in cmd.get_model(sel).atom:
        cmd.alter("%s and resi %s and name %s"%(sel,a.resi,a.name),"vdw=%f"%(a.q))
        
cmd.extend("useOccRadii",useOccRadii)


for a in cmd.get_model('tmp').atom:
    print a.resi_number
