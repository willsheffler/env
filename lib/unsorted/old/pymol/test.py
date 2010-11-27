import pymol
from pymol import cmd
import sys,os,random

#PyMOL>print dir(cmd.get_model('r1').atom[0])
#['__cmp__', '__doc__', '__getattr__', '__module__', 'b', 'bohr', 'chain', 'coord', 'defaults', 'flags', 'formal_charge', 'get_implicit_valence', 'get_mass', 'get_number', 'get_signature', 'has', 'hetatm', 'id', 'in_same_residue', 'index', 'name', 'new_in_residue', 'numeric_type', 'partial_charge', 'q', 'resi', 'resi_number', 'resn', 'segi', 'ss', 'stereo', 'symbol', 'vdw']

def useRosettaRadii():
	cmd.alter("element C", "vdw=2.00")
	cmd.alter("element N", "vdw=1.75")
	cmd.alter("element O", "vdw=1.55")
	cmd.alter("element H", "vdw=1.00")
	cmd.alter("element P", "vdw=1.90")
	cmd.set("sphere_scale", 1.0)
	
cmd.extend('useRosettaRadii', useRosettaRadii)

def useOccRadii(sel="all"):
	for a in cmd.get_model(sel).atom:
		q = a.q
		if q >= 2:
			print "shrik radius"
			q <- 0.1
		cmd.alter("%s and resi %s and name %s"%(sel,a.resi,a.name),"vdw=%f"%(q))

cmd.extend("useOccRadii",useOccRadii)

def which(v,x):
	r = []
	for ii in range(len(v)):
		if v[ii]==x:
			r.append(ii)
	return r

## cmd.reinitialize()
## cmd.do("")

## wd = "/users/sheffler/project/modelpicking/pdbs/"
## pdbs = os.listdir(wd)
## pdbs = [wd+p for p in pdbs if p.startswith('cf') and p.endswith('.pdb')]

## cmd.create('n','none')
## cmd.load(wd+'1hz6_0001.pdb','n')
## #cmd.hide('lines')
## #cmd.show("cartoon")

## for ii in range(50):#len(pdbs)):
## 	print ii
## 	p = pdbs[ii]
## 	cmd.create('d'+`ii`,'none')
## 	cmd.load(p,'d'+`ii`)
## 	cmd.hide('all')
## 	#cmd.show("cartoon")

## cmd.do("dss")
## cmd.show('lines','name ca or name n or name c')
