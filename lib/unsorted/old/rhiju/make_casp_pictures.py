#!/usr/bin/python

## used to be ~/python/2002/scratch11_29.py

## run maxsub, get superpositions

import string
from glob import glob
from os import system, popen
from amino_acids import extra_longer_names
from popen2 import popen2
from blast import NBAlign
from sys import argv

def Rclean(rin,rout):
    rin.write('select none\n')
    rin.flush()
    line = rout.readline()
    while line and (len(line)<10 or line[:10] != 'No atoms s'):
        #print line[:-1]
        line = rout.readline()
    return

def Mammoth_align(data): ## this version returns empty if z<MIN_Z_SCORE

    BLOCK = 10
    NBLOCK = 5
    MIN_Z_SCORE = -20

    line = data.readline()
    while line and line[:4] not in ['Z-sc',' Z-s']:
        line = data.readline()

    if len(string.split(line))<2:
        data.close()
        log('Bad mammoth data\n')
        return 0.0,'','','','',{}

    z = float(string.split(line)[1])

    if z<MIN_Z_SCORE:
        data.close()
        return z,'','','','',{}

    while line and line[0:4] != 'Pred':line = data.readline()
    al1 = ''
    al2 = ''
    bars = ''
    stars = ''
    ss1 = ''
    ss2 = ''
    align = []
    while line and line[0:4] == 'Pred':
        #print line[:-1]
        line = line[11:-1] + ' '*(NBLOCK*BLOCK+NBLOCK+10 - len(line) +1)
        #print line
        for i in range(NBLOCK):
            for j in range(BLOCK):
                #print j+i*BLOCK+i,len(line)
                al1 = al1 + line[j+i*BLOCK+i]
        line = data.readline()
        line = line[11:-1] + ' '*(NBLOCK*BLOCK+NBLOCK+10 - len(line) +1)
        for i in range(NBLOCK):
            for j in range(BLOCK):
                #print j+i*BLOCK+i,len(line)
                ss1 = ss1 + line[j+i*BLOCK+i]
        line = data.readline()
        line = line[11:-1] + ' '*(NBLOCK*BLOCK+NBLOCK+10 - len(line) +1)
        for i in range(NBLOCK):
            for j in range(BLOCK):
                bars = bars + line[j+i*BLOCK+i]
        line = data.readline()
        line = line[11:-1] + ' '*(NBLOCK*BLOCK+NBLOCK+10 - len(line) +1)
        for i in range(NBLOCK):
            for j in range(BLOCK):
                #print j+i*BLOCK+i,len(line)
                ss2 = ss2 + line[j+i*BLOCK+i]
        line = data.readline()
        line = line[11:-1] + ' '*(NBLOCK*BLOCK+NBLOCK+10 - len(line) +1)
        for i in range(NBLOCK):
            for j in range(BLOCK):
                al2 = al2 + line[j+i*BLOCK+i]
        line = data.readline()
        line = line[11:-1] + ' '*(NBLOCK*BLOCK+NBLOCK+10 - len(line) +1)
        for i in range(NBLOCK):
            for j in range(BLOCK):
                stars = stars + line[j+i*BLOCK+i]
        line = data.readline()
        line = data.readline()
        line = data.readline()
    data.readlines()
    data.close()

    assert string.count(al1,'-')==0
    assert len(al1) == len(al2) == len(stars) == len(bars) == len(ss1) == len(ss2)

    align = []
    for i in range(len(al1)):
        if stars[i] == '*':
            ## numbering starts at 0
            pos1 = i - string.count(al1[:i],'.')
            pos2 = i - string.count(al2[:i],'.')
            align.append([pos1,pos2])

    seq1 = string.join(string.split(string.join(string.split(al1,'.'),''),' '),'')
    ss1 = string.join(string.split(string.join(string.split(ss1,'.'),''),' '),'')
    seq2 = string.join(string.split(string.join(string.split(al2,'.'),''),' '),'')
    ss2 = string.join(string.split(string.join(string.split(ss2,'.'),''),' '),'')
    return z,seq1,seq2,ss1,ss2,align
################################################################

if 1:
    p_file = argv[1]
    e_file = argv[2]

    files = [p_file,e_file]


    R = float(argv[3])

    system('/users/pbradley/maxsub/maxsub -e %s -p %s -R %f -o %s'\
           %(e_file,p_file,R,
             p_file+'.mammoth'))

    ## read the transformations
    file = 'maxsub_sup.pdb'
    matrix = map(lambda x:map(float,string.split(x)[1:]),
                 popen('grep -A3 "Transformation Matrix" %s'%file).readlines()[1:])

    translation = {}
    translation[p_file] = map(float,
             string.split(popen('grep -A1 "Translation vector (Pred" %s'%file).readlines()[1])[1:])

    translation[e_file] = map(float,
             string.split(popen('grep -A1 "Translation vector (Exp" %s'%file).readlines()[1])[1:])




    def P_transform(v,matrix,t):
        ans = [0.0]*3
        for i in range(3):
            ans[i] = v[i] - t[i]
        return ans

    def E_transform(v,matrix,t):
        ans = [0.0]*3
        for i in range(3):
            for j in range(3):
                ans[i] = ans[i] + matrix[i][j]*(v[j]-t[j])


        return ans

    Transform = {p_file:P_transform,
                 e_file:E_transform}


    seq = {}
    numbers = {}
    for file in files:
        lines = popen('grep "^ATOM.*CA" '+file).readlines()
        seq[file] = ''
        numbers[file] = []

        for line in lines:
            rsd = line[17:20]
            number = int(line[22:26])
            if rsd in extra_longer_names.keys():
                seq[file] = seq[file] + extra_longer_names[rsd]
            else:
                seq[file] = seq[file] + 'X'
            numbers[file].append(number)

    s1 = seq[e_file]
    s2 = seq[p_file]

    align_E2P = NBAlign(s1,s2)
    print 'len(e_seq) %d len(p_seq) %d aligned: %d'\
          %(len(s1),len(s2),len(align_E2P.keys()))



    good_numbers = {}
    good_numbers[e_file] = {}
    good_numbers[p_file] = {}
    for i in align_E2P.keys():
        good_numbers[e_file][numbers[e_file][i]] = numbers[e_file][i]
        good_numbers[p_file][numbers[p_file][align_E2P[i]]] = numbers[e_file][i]



    chain = 'A'

    for file in files:

        out = open(file+'.sup','w')
        lines = popen('grep "^ATOM" '+file).readlines()
        for line in lines:
            number = int(line[22:26])
            if good_numbers[file].has_key(number):
                pos = map(float,[line[30:38],
                                 line[38:46],
                                 line[46:54]])

                pos = Transform[file](pos,matrix,translation[file])
                new_number = good_numbers[file][number]

                out.write('%s%s%4d%s%8.3f%8.3f%8.3f%s'\
                          %(line[:21],chain,new_number,line[26:30],
                            pos[0],pos[1],pos[2],
                            line[54:]))
        out.close()

    ## read the mammoth alignment, map back to seq[e_file]
    z,seq1,seq2,ss1,ss2,align = Mammoth_align(open(p_file+'.mammoth','r'))

    assert string.count(seq[e_file],seq2)
    offset = string.find(seq[e_file],seq2)

    aligned_numbers = []
    for pos in map(lambda x:x[1]+offset,align): ## starts at 0
        if pos in align_E2P.keys():
            aligned_numbers.append( numbers[e_file][pos] )



    ## get the right orientation, write a molscript file for p
    rout2,rin2 = popen2('/users/pbradley/rasmol_16BIT')
    rin2.write('load %s\nstructure\ncartoons\nset background grey\ncolor group\n'%(e_file+'.sup'))
    rin2.write('wireframe off\n')
    rin2.flush()
    Rclean(rin2,rout2)

    rout,rin = popen2('/users/pbradley//rasmol_16BIT')
    rin.write('load %s\nstructure\ncartoons\nset background white\ncolor group\n'%(p_file+'.sup'))
    rin.write('wireframe off\n')
    rin.flush()

    ## color aligned portions
    ual = 0
    for pos in good_numbers[p_file].values():
        if pos not in aligned_numbers:
            if not ual:
                rin.write('select %d\ndefine ual selected\n'%pos)
            else:
                rin.write('select %d or ual\ndefine ual selected\n'%pos)
            ual = ual+1
            if not ual%10:
                Rclean(rin,rout)

    rin.write('select ual\ncartoons off\nbackbone 75\n')


    Rclean(rin,rout)
    raw_input()
    rin.write('select all\nbackbone off\ncartoons\ncolor group\n')
    rin.write('write molscript %s\n'%(p_file+'.mol_temp'))
    rin.flush()

    ## write a molscript file for e
    rin.write('zap\nload %s\nstructure\ncartoons\nset background white\ncolor group\n'%(e_file+'.sup'))
    rin.write('wireframe off\n')
    rin.write('write molscript %s\n'%(e_file+'.mol_temp'))

    Rclean(rin,rout)
    rin.write('quit\n')
    Rclean(rin,rout)

    rin2.write('quit\n')
    Rclean(rin2,rout2)

    ## parse out from the p-molscript file: orientation info, ss info

    head = {}
    foot = {}
    transform = {}

    for file in files:
        lines = open(file+'.mol_temp','r').readlines()


        current = []
        in_transform = 0
        for line in lines:
            if string.count(line,'transform'):
                head[file] = current[:]
                in_transform = 1
                current = []
            elif (in_transform and line[:6] != '    by'):
                transform[file] = current[:]
                current = []
                in_transform = 0

            current.append(line)

        foot[file] = current
        transform[file].append('  set residuecolour amino-acids rainbow, colourparts on;\n')

        for pos in good_numbers[file].values():
            if pos not in aligned_numbers:
                transform[file].append('  set residuecolour residue A%d rgb 0.9 0.9 0.9;\n'\
                                       %pos)


    transform[e_file] = transform[p_file]
    head[e_file][6] = head[p_file][6]

    for file in files:
        out = open(file+'.mol','w')
        out.writelines(head[file]+transform[file]+foot[file])
        out.close()

        out = open(file+'.mol_rainbow','w')
        for line in head[file]+transform[file]+foot[file]:
            if not string.count(line,'residuecolour residue'):
                out.write(line)
        out.close()

        for base in [file+'.mol',file+'.mol_rainbow']:

            system('/users/rohl/bin/molscript -ps < %s > %s'\
                   %(base,
                     base+'.ps'))

            system('/users/rohl/bin/molscript -r < %s > %s'\
                   %(base,
                     base+'.r3d'))

            size_tag = ''
            size_tag = ' -size 1000x1000'
            system('/users/pbradley/raster3d/Raster3D_2.6e/render -jpeg %s < %s > %s'\
                   %(size_tag,
                     base+'.r3d',
                     base+'.r3d.jpeg'))

    for tag in ['.mol','.mol_rainbow']:
        system('cp %s %s'\
               %(e_file+'%s.ps'%tag,
                 p_file+'.native%s.ps'%tag))

        system('cp %s %s'\
               %(e_file+'%s.r3d'%tag,
                 p_file+'.native%s.r3d'%tag))

        system('cp %s %s'\
               %(e_file+'%s.r3d.jpeg'%tag,
                 p_file+'.native%s.r3d.jpeg'%tag))

        system('cp %s %s'\
               %(e_file+tag,
                 p_file+'.native'+tag))

