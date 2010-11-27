#!/usr/bin/python

## An edit of Phil's awesome "make_casp_pictures.py" which
## has examples of pretty much everything you'd need to do
## to make nice molscript pictures.

## used to be ~/python/2002/scratch11_29.py

## run maxsub, get superpositions

import string
from glob import glob
from os import system, popen
from amino_acids import extra_longer_names
from popen2 import popen2
from blast import NBAlign
from sys import argv


def Help():
    print '%'*70
    print 'Usage: make_aligned_pictures filelist rmsdcutoff'
    print ''
    print 'Input a filelist and rmsd cutoff for maxsub, and'
    print 'you''ll get 01.ps, 02.ps, ... molscript cartoons'
    print ''
    print '%'*70



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

#############################################################
if len(argv)<4:
    Help()

if 1:
    filelist = argv[1]
    fid = open( filelist, 'r')
    allfiles = fid.readlines()


    #In Phil's original script, the prediction (p_file) was first, before the
    # experimental file (e_file). I've kind of got it backwards where the
    # native should go first, and all others later. Anyway...

    p_file = allfiles[0][:-1] #Get rid of newline at end of line.

    R = float(argv[2])

    good_numbers = {}
    translation = {}
    seq = {}
    numbers = {}
    matrix={p_file:[]}

#Define a list of all the files that need to be aligned
    e_files = []
    for file in allfiles[1:]:
        e_files.append(  file[:-1])

    Transform = {p_file:P_transform}

    for e_file in e_files:
       files = [p_file, e_file]
       Transform[e_file] = E_transform

       e_file_nopath = string.split(e_file,'/')[-1]
       system('/users/pbradley/maxsub/maxsub -e %s -p %s -R %f -o %s'\
              %(e_file,p_file,R,
                e_file_nopath+'.mammoth'))

    ## read the transformations
       file = 'maxsub_sup.pdb'
       matrix[e_file] = map(lambda x:map(float,string.split(x)[1:]),
                    popen('grep -A3 "Transformation Matrix" %s'%file).readlines()[1:])

       translation[p_file] = map(float,
                                 string.split(popen('grep -A1 "Translation vector (Pred" %s'%file).readlines()[1])[1:])

       translation[e_file] = map(float,
                                 string.split(popen('grep -A1 "Translation vector (Exp" %s'%file).readlines()[1])[1:])

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

       good_numbers[e_file] = {}
       good_numbers[p_file] = {}
       for i in align_E2P.keys():
           good_numbers[e_file][numbers[e_file][i]] = numbers[e_file][i]
           good_numbers[p_file][numbers[p_file][align_E2P[i]]] = numbers[e_file][i]

       chain = 'A'


    files = [p_file]+e_files
    for file in files:
        file_nopath = string.split(file,'/')[-1]
        out = open(file_nopath+'.sup','w')
        lines = popen('grep "^ATOM" '+file).readlines()
        for line in lines:
            number = int(line[22:26])
            if good_numbers[file].has_key(number):
                pos = map(float,[line[30:38],
                                 line[38:46],
                                 line[46:54]])

                pos = Transform[file](pos,matrix[file],translation[file])
                new_number = good_numbers[file][number]

                out.write('%s%s%4d%s%8.3f%8.3f%8.3f%s'\
                          %(line[:21],chain,new_number,line[26:30],
                            pos[0],pos[1],pos[2],
                            line[54:]))
        out.close()


    ## read the mammoth alignment, map back to seq[e_file]
     #  z,seq1,seq2,ss1,ss2,align = Mammoth_align(open(p_file+'.mammoth','r'))

     #  assert string.count(seq[e_file],seq2)
     #  offset = string.find(seq[e_file],seq2)

     #  aligned_numbers = []
     #  for pos in map(lambda x:x[1]+offset,align): ## starts at 0
     #      if pos in align_E2P.keys():
     #          aligned_numbers.append( numbers[e_file][pos] )



    ## get the right orientation, write a molscript file for p
#    rout2,rin2 = popen2('/users/pbradley/rasmol_16BIT')
#    rin2.write('load %s\nstructure\ncartoons\nset background grey\ncolor group\n'%(e_file+'.sup'))
#    rin2.write('wireframe off\n')
#    rin2.flush()
#    Rclean(rin2,rout2)

    rout,rin = popen2('/users/pbradley//rasmol_16BIT')
    file_nopath = string.split(p_file,'/')[-1]
    rin.write('load %s\nstructure\ncartoons\nset background white\ncolor group\n'%(file_nopath+'.sup'))
    rin.write('wireframe off\n')
    rin.flush()

    ## color aligned portions
#    ual = 0
#    for pos in good_numbers[p_file].values():
#        if pos not in aligned_numbers:
#            if not ual:
#                rin.write('select %d\ndefine ual selected\n'%pos)
#            else:
#                 rin.write('select %d or ual\ndefine ual selected\n'%pos)
#             ual = ual+1
#             if not ual%10:
#                 Rclean(rin,rout)

#    rin.write('select ual\ncartoons off\nbackbone 75\n')


    Rclean(rin,rout)
    raw_input()
    rin.write('select all\nbackbone off\ncartoons\ncolor group\n')

    file_nopath = string.split(p_file,'/')[-1]
    rin.write('write molscript %s\n'%(file_nopath+'.mol_temp'))
    rin.flush()

    ## write a molscript file for e
    for e_file in e_files:
       file_nopath = string.split(e_file,'/')[-1]
       print file_nopath
       rin.write('zap\nload %s\nstructure\ncartoons\nset background white\ncolor group\n'%(file_nopath+'.sup'))
       rin.write('wireframe off\n')
       file_nopath = string.split(e_file,'/')[-1]
       rin.write('write molscript %s\n'%(file_nopath+'.mol_temp'))
       Rclean(rin,rout)

    rin.write('quit\n')
    Rclean(rin,rout)

 #   rin2.write('quit\n')
 #   Rclean(rin2,rout2)

    ## parse out from the p-molscript file: orientation info, ss info

    head = {}
    foot = {}
    transform = {}

    files = [p_file]+e_files
    for file in files:
        file_nopath = string.split(file,'/')[-1]
        lines = open(file_nopath+'.mol_temp','r').readlines()


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
            if string.count(line,'window'):
                current.append('  frame off;')

        score  = 0.
        rmsd   = 0.
        maxsub = len(s1)
        env  = 0
        pair = 0
        hs   = 0
        ss   = 0

        temp =  popen('grep rms ' + file + '| grep -v maxsub').readlines()
        if len(temp)>0:
            rmsd = float(string.split(temp[0][:-1])[1]); #First line, ignore newline at end, second word.

        temp =  popen('grep maxsub ' + file).readlines()
        if len(temp)>0:
            maxsub = float(string.split(temp[0][:-1])[2]);

        temp =  popen('grep score ' + file).readlines()
        if len(temp)>0:
            score = float(string.split(temp[0][:-1])[1]);

        temp =  popen('grep env ' + file).readlines()
        if len(temp)>0:
            env = float(string.split(temp[0][:-1])[1]);

        temp =  popen('grep pair ' + file).readlines()
        if len(temp)>0:
            pair = float(string.split(temp[0][:-1])[1]);

        temp =  popen('grep hs ' + file).readlines()
        if len(temp)>0:
            hs = float(string.split(temp[0][:-1])[1]);

        temp =  popen('grep ss ' + file).readlines()
        if len(temp)>0:
            ss = float(string.split(temp[0][:-1])[1]);

        foot[file]=[]
        foot[file].append('  set labelsize 40;\n')
        foot[file].append('  label 0 -20 0 "%s";\n'%file)

        otherusefulcrap = 'Score %6.1f; RMSD %5.2f; MaxSub %d' % (score,rmsd,maxsub)
        foot[file].append('  label 0 25 0 "%s";\n'%otherusefulcrap)

        otherusefulcrap = 'Env %5.1f; Pair %5.1f; HS %5.1f; SS %5.1f' % (env,pair,hs,ss)
        foot[file].append('  label 0 22 0 "%s";\n'%otherusefulcrap)


        foot[file] = foot[file]+current
        transform[file].append('  set residuecolour amino-acids rainbow, colourparts on;\n')


    for e_file in e_files:
        transform[e_file] = transform[p_file]
        head[e_file][6] = head[p_file][6]


    count=0
    for file in files:
        file_nopath = string.split(file,'/')[-1]

        out = open(file_nopath+'.mol_rainbow','w')
        for line in head[file]+transform[file]+foot[file]:
            if not string.count(line,'residuecolour residue'):
                out.write(line)

        out.close()

        base = file_nopath+'.mol_rainbow';

        count = count+1
        system('/users/rohl/bin/molscript -ps < %s > %02d.ps'\
               %(base,
                 count))

#Now get rid of all the crap files that have been generated!
    for file in files:
        file_nopath = string.split(file,'/')[-1]
        system('rm '+file_nopath+'.mammoth')
        system('rm '+file_nopath+'.mol_rainbow')
        system('rm '+file_nopath+'.mol_temp')
        system('rm '+file_nopath+'.sup')
