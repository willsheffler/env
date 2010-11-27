import sys,os,re,math
from VMD import *
from molecule import *
from Molecule import *
from AtomSel import AtomSel

#print sys.version

mols = []
nativeSS = []

WIDTH  = 4.0#1.9
HEIGHT = 3.0
if display.get("stereo") != "Off":
    WIDTH /= 2
    
def stereon():
    if dislpay.get("stereo") == "Off":
        WIDTH = WIDTH/2
        display.set('stereo',"SideBySide")
        
def stereoff():
    if display.get("stereo") != "Off":
        WIDTH = WIDTH*2
        display.set('stereo',"Off")
        
    

def mean(x):
    t = 0
    for v in x:
        t += v
    return t/len(x)

def test():
    print 'test!'


def reset():
    for m in mols:
        trans.resetview(m)
    for m in mols:
        trans.set_center(m,[0,0,0])


def show(first,*args):
    if type(first) == type([]) or type(first) == type((1,)):
          args = list(first)+list(args)
    else: args = [first]+list(args)
    for ii in range(len(mols)):
	trans.show(int(mols[ii]),0)
    for ii in args:
	trans.show(mols[ii],1)
    if len(args) > 1:
        for ii in range(len(args)):
            mols[args[ii]].reps()[0].changeColor("colorid "+`ii`)
    #for ii in args:
    #    print ii, trans.get_center(mols[ii])
    
def showpairs(val):
    for m in mols:
        n = molrep.num(int(m))
        for r in range(1,n):
            molrep.set_visible(int(m),r,val)

def changerep(val):
    for m in mols:
        molrep.modrep(int(m),0,style=val)

def changecolor(val):
    for m in mols:
        molrep.modrep(int(m),0,color=val)


def spread(*args):
    reset()
    if len(args) == 0:
        args = range(len(mols))
    n = len(args)
    if HEIGHT > WIDTH:
        h = int(math.ceil(math.sqrt(n)))+1
        w = int(math.ceil(float(n)/(h-1)))+1
    else:
        w = int(math.ceil(math.sqrt(n)))+1
        h = int(math.ceil(float(n)/(w-1)))+1
    print n , h , w
    X = [WIDTH *(-0.5+x/float(w-1)) for x in range(w)]
    Y = [HEIGHT*(0.5-y/float(h-1)) for y in range(h)]
    print X
    print Y
    for ii in range(len(args)):
        ix = (ii)%(w-1) 
        x  = (X[ix]+X[ix+1])/2
        iy = int(math.floor((ii)/(w-1)))
        y  = (Y[iy]+Y[iy+1])/2
        print x,y 
	m = int(mols[args[ii]])
	#print ii , m
	for jj in args:
	    trans.fix(int(mols[jj]),1)
	trans.fix(int(mols[args[ii]]),0)
        #s = trans.get_scale(int(mols[ii]))
        #trans.set_scale(int(mols[ii]),  1.0/w)
        #if WIDTH > HEIGHT:            x,y = y,x
	trans.translate(x,y,0)

    for ii in args:
	trans.fix(int(mols[ii]),0)
        mols[ii].reps()[0].changeColor('index')
    for ii in range(len(mols)):
        if ii in args:
            trans.show(int(mols[ii]),1)
        else:
            trans.show(int(mols[ii]),0)   

    trans.scale(1.40/w)


def loadSS(dsspfile):
    global nativeSS
    lines = open(dsspfile).readlines()
    for ii in range(len(lines)):
        del lines[0]
        if lines[0].split()[:4] == ['#','RESIDUE','AA','STRUCTURE']:
            break    
    res = [line[6:10] for line in lines]
    nativeSS  = [line[16]   for line in lines]
    for ii in range(len(nativeSS)):
        if not nativeSS[ii] in ("H","E"):
            nativeSS[ii] = "L"
    print "LOAD SS",nativeSS

def setBeta(beta):
    print len(beta)
    for m in mols:
        for ii in range(len(beta)):
            sel = AtomSel('resid '+`ii`,int(m))
            sel.set('beta',beta[ii])
            del sel

def colorss():
    beta = [0]*len(nativeSS)
    for ii in range(len(nativeSS)):
        if nativeSS[ii] == 'H': beta[ii] = 0.0
        if nativeSS[ii] == 'E': beta[ii] = 0.4
        if nativeSS[ii] == 'L': beta[ii] = 1.0
    setBeta(beta)
    changecolor('beta')

def loadJumpPDB(pdbfile,pairfile,nativepdb=None,recenter=1,ref=None):
    res = []

#(backbone and resid 27 and name O C) or (backbone and resid 100 and name O C)
#(backbone and resid 27 and name N H) or (backbone and resid 100 and name N H)

    pdb = open(pdbfile).readlines()
    cutpoints = []
    for line in pdb:
        if line[:17] == "REMARK cutpoints:":
            for cp in line[17:].split():
                cutpoints.append(int(cp))

	#print cutpoints
    
    def jpdonsel(pair):
        pos1   = int(pair[0])
        pos2   = int(pair[1])
        orient = pair[2]
        pleat  = int(pair[3])
        if orient in ('1','A'):
            #print "anti-parallel",pair[3]
            #if pair[3] == 'P' or pair[3] == 2:
            return "(resid %i and name N H) or (resid %i and name N H)"%(pos1,pos2)
        else:
            #print "parallel",pair[3]
            if pleat == 1:
                #print 'pleat 1'
                return "(resid %i and name N H) or (resid %i and name N H)"%(pos1+1,pos2)
            else:
                #print 'pleat 2'
                return "(resid %i and name N H) or (resid %i and name N H)"%(pos1,pos2+1)
    
    def jpaccsel(pair):
        pos1   = int(pair[0])
        pos2   = int(pair[1])
        orient = pair[2]
        pleat  = int(pair[3])
        if orient in ('1','A'):
            #print "anti-parallel",pair[3]
            #if pair[3] == 'P' or pair[3] == 2:
            if pleat == 1:
                return "(resid %i and name O C) or (resid %i and name O C)"%(pos1,pos2)
            else:
                return "(resid %i and name O C) or (resid %i and name O C)"%(pos1,pos2)
        else:
            #print "parallel",pleat
            if pleat == 1:
                #print 'pleat 1'
                return "(resid %i and name O C) or (resid %i and name O C)"%(pos1-1,pos2)
            else:
                #print 'pleat 2'
                return "(resid %i and name O C) or (resid %i and name O C)"%(pos1,pos2-1)

    def jpsel(pair):
	a = int(pair[1])
	b = int(pair[0])
	c,d = a,a
	if pair[2] == 'P' or pair[2] == 2:
	    c,d = a+1,a-1
	sels = []
	sels.append("resid %i and name N"%c)
	sels.append("resid %i and name N"%b)
	sels.append("resid %i and name C"%d)
	sels.append("resid %i and name C"%b)
	return sels
    
    def jpsssel(ss,lhe):
        h = []
        for ii in range(len(ss)):
            if ss[ii] == lhe:
                h.append(ii)
        return "backbone and resid "+" ".join([str(i) for i in h])
    
    # print jpsssel(ss,'H')
    
    pairings = [x.split() for x in open(pairfile).read().split('\n') if len(x.split())==4]
    # print pairings
    
    pdb = open(pdbfile).readlines()
    
    m = Molecule()
    m.load(pdbfile)
    m.clearReps()

    if ref is not None:
        print ref._AtomSel__text
        sel = AtomSel(ref._AtomSel__text,m)
        sel.align(ref)
    elif recenter:
	sel   = AtomSel("all",m)
	x,y,z = sel.get('x','y','z')
	mx,my,mz = mean(x), mean(y), mean(z)
	x = [t - mx for t in x]
	y = [t - my for t in y]
	z = [t - mz for t in z]
	sel.set('x',x)
	sel.set('y',y)
	sel.set('z',z)
	print AtomSel("backbone",m).center()
	print trans.set_center(m,[0,0,0])

    repbb  = MoleculeRep('lines 2',      'index',   'backbone'   )
    m.addRep(repbb)
    
    repdonjp = []
    repaccjp = []
    represjp = []
    for ii in range(len(pairings)):
	c = `(3,4,6,7,8,9,10,12,13,14)[ii]`
        print pairings[ii]
        repdonjp.append(MoleculeRep('vdw 0.5',     'colorid '+c,
                                 jpdonsel(pairings[ii])))
        repaccjp.append(MoleculeRep('licorice 0.4','colorid '+c,
                                 jpaccsel(pairings[ii])))
        represjp.append(MoleculeRep('licorice 0.2','colorid '+c,
                                    "backbone and resid "+pairings[ii][0]+" "+pairings[ii][1]))
        m.addRep(repdonjp[ii])
        m.addRep(repaccjp[ii])
        m.addRep(represjp[ii])

    for ii in range(len(pairings)):
	c = `(3,4,6,7,8,9,10,12,13,14)[ii]`
	xyz = [tuple([x[0] for x in AtomSel(s,m).get('x','y','z')]) for s in jpsel(pairings[ii])]
	graphics.color(m,c)
        if pairings[ii][2] == 'P' or pairings[ii][2] == 2:
            print "TRI parallel"
            graphics.triangle(m,xyz[0],xyz[1],xyz[2])
            graphics.triangle(m,xyz[0],xyz[1],xyz[3])
        else:
            graphics.triangle(m,xyz[0],xyz[1],xyz[2])
            graphics.triangle(m,xyz[1],xyz[2],xyz[3])
	#print xyz

    for cp in cutpoints:
	v1 = AtomSel('resid %i and name C'%(cp),m).get('x','y','z')
	v1 = ( v1[0][0], v1[1][0], v1[2][0] )
	v2 = AtomSel('resid %i and name N'%(cp+1),m).get('x','y','z')
	v2 = ( v2[0][0], v2[1][0], v2[2][0] )
	graphics.color(m,'white')
	graphics.cylinder(m,v1,v2,radius=0.4,filled=1)

    return m




def visualize(which):
    global mols
    if which == 'rhiju':

        pairfile = 'rhiju1iib/threepair_rev.pdat'
        mols.append(loadJumpPDB('rhiju1iib/1iib.pdb',pairfile))
        pairings = [x.split() for x in open(pairfile).read().split('\n') if len(x.split())==4]
        res = []
        ext = 2
        for p in pairings:
            for x in range(-ext,ext+1):
                res.append(str(int(p[0])+x))
                res.append(str(int(p[1])+x))
        ref = AtomSel("backbone and resid "+" ".join(res),mols[0])        
        for ii in range(1,3):
            mols.append(loadJumpPDB('rhiju1iib/xx1iib%04i.pdb'%ii,
                                    pairfile,
                                    ref=ref))#
        for ii in range(9,12):
            mols.append(loadJumpPDB('rhiju1iib/S_%i.pdb'%ii,
                                    pairfile,
                                    ref=ref))#
        #showpairs(0)
    
    else:
        pairfile = '/users/sheffler/project/jumping/1fkb/1fkb.pairings'
        natfile  = '/users/sheffler/project/jumping/1fkb/1fkb_/1fkb.pdb'
        natdssp  = '/users/sheffler/project/jumping/1fkb/1fkb_/1fkb.dssp'
        loadSS(natdssp)
        mols.append(loadJumpPDB(natfile,pairfile))
        #mols.append(loadJumpPDB(natfile,pairfile))
        #mols.append(loadJumpPDB(natfile,pairfile))
        #mols.append(loadJumpPDB(natfile,pairfile))
        #ref = None
        ref = AtomSel('backbone',mols[0])
        mols.append(loadJumpPDB("results/ap1fkb_good/ap1fkb0304.pdb",pairfile,ref=ref))
        mols.append(loadJumpPDB("results/ap1fkb_good/ap1fkb1931.pdb",pairfile,ref=ref))
        #mols.append(loadJumpPDB("results/ap1fkb_good/ap1fkb2694.pdb",pairfile,ref=ref))
        #mols.append(loadJumpPDB("results/ap1fkb_bad/ap1fkb0157.pdb",pairfile,ref=ref))
        #mols.append(loadJumpPDB("results/ap1fkb_bad/ap1fkb1268.pdb",pairfile,ref=ref))
        mols.append(loadJumpPDB("results/ap1fkb_bad/ap1fkb2287.pdb",pairfile,ref=ref))
        mols.append(loadJumpPDB("results/ap1fkb_rand/ap1fkb0157.pdb",pairfile))
        #mols.append(loadJumpPDB("results/ap1fkb_rand/ap1fkb1268.pdb",pairfile))
        #mols.append(loadJumpPDB("results/ap1fkb_rand/ap1fkb2287.pdb",pairfile))
        
    
    for m in mols:
        trans.set_center(int(m),[0,0,0])
    

