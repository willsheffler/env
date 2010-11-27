## tools
from os import chdir,remove,popen
import string
from time import time
from os.path import exists
from whrandom import random


def bl2seq(seq1,seq2):
    chdir('/users/pbradley/blast/')
    
    out = open('junk1_bl2.seq','w')
    out.write('>s1\n'+seq1+'\n')
    out.close()
    out = open('junk2_bl2.seq','w')
    out.write('>s2\n'+seq2+'\n')
    out.close()
    lines = popen('./bl2seq -i junk1_bl2.seq -j junk2_bl2.seq -p blastp '+\
                  '| grep "Identities"').readlines()
    if lines and string.count(lines[0],'Ident'):
        id = float(string.split(lines[0])[3][1:-3])
    else:
        #print 0,seq1
        #print 0,seq2
        id = 0
    chdir('/users/pbradley/python/')
    return id

    
def al2seq(seq1,seq2):


    file1 = '/scratch/Phil/al2seq'+str(random())
    file2 = '/scratch/Phil/al2seq'+str(random())
    out = open(file1,'w')
    out.write('>s1\n'+seq1+'\n')
    out.close()
    out = open(file2,'w')
    out.write('>s2\n'+seq2+'\n')
    out.close()
    lines = popen('/users/pbradley/blast/bl2seq -p blastp -F F -i '+file1+\
                  ' -j '+file2).readlines()
    remove(file1)
    remove(file2)
    al  = []
    a1 = ''
    a2 = ''
    start1 = 1000000
    start2 = 1000000
    done = 0
    id = 0
    for line in lines:
        if line[0:10] == ' Identitie':
            if done: break
            id = int(string.split(line)[3][1:-3])
            done = 1
        elif line[0:6] == 'Query:':
            l = string.split(line)
            assert len(l) == 4
            a1 = a1 + l[2]
            start1 = min(start1,int(l[1]))
        elif line[0:6] == 'Sbjct:':
            l = string.split(line)
            assert len(l) == 4
            a2 = a2 + l[2]
            start2 = min(start2,int(l[1]))
            
    s1 = string.join(string.split(a1,'-'),'')
    s2 = string.join(string.split(a2,'-'),'')
    off1 = string.find(seq1,s1)
    off2 = string.find(seq2,s2)
    if 0:
        print a1
        print s1
        print seq1
        print '---'
        print a2
        print s2
        print seq2
        print off1,off2
    assert off1 >= 0 and off1 == start1-1
    assert off2 >= 0 and off2 == start2-1
    for i in range(len(a1)):
        if a1[i] != '-' and a2[i] != '-':
            pos1 = i + off1 - string.count(a1[:i],'-')
            pos2 = i + off2 - string.count(a2[:i],'-')
            al.append([pos1,pos2])
    return al,id

    
def al2seqD(seq1,seq2):
    chdir('/users/pbradley/blast/')
    
    out = open('/scratch/Phil/junk1_al2.seq','w')
    out.write('>s1\n'+seq1+'\n')
    out.close()
    out = open('/scratch/Phil/junk2_al2.seq','w')
    out.write('>s2\n'+seq2+'\n')
    out.close()
    lines = popen('./bl2seq -F F -i /scratch/Phil/junk1_al2.seq '+\
                  '-j /scratch/Phil/junk2_al2.seq -p blastp').readlines()
    al  = {}
    a1 = ''
    a2 = ''
    start1 = 1000000
    start2 = 1000000
    done = 0
    id = 0
    for line in lines:
        if line[0:10] == ' Identitie':
            if done: break
            id = int(string.split(line)[3][1:-3])
            done = 1
        elif line[0:6] == 'Query:':
            l = string.split(line)
            assert len(l) == 4
            a1 = a1 + l[2]
            start1 = min(start1,int(l[1]))
        elif line[0:6] == 'Sbjct:':
            l = string.split(line)
            assert len(l) == 4
            a2 = a2 + l[2]
            start2 = min(start2,int(l[1]))
            
    s1 = string.join(string.split(a1,'-'),'')
    s2 = string.join(string.split(a2,'-'),'')
    off1 = string.find(seq1,s1)
    off2 = string.find(seq2,s2)
    if 0:
        print a1
        print s1
        print seq1
        print '---'
        print a2
        print s2
        print seq2
        print off1,off2
    assert off1 >= 0 and off1 == start1-1
    assert off2 >= 0 and off2 == start2-1
    for i in range(len(a1)):
        if a1[i] != '-' and a2[i] != '-':
            pos1 = i + off1 - string.count(a1[:i],'-')
            pos2 = i + off2 - string.count(a2[:i],'-')
            al[pos1] = pos2
    remove('/scratch/Phil/junk1_al2.seq')
    remove('/scratch/Phil/junk2_al2.seq')
    chdir('/users/pbradley/python/')
    return al

    
def allal2seq(seq1,seq2):
    chdir('/users/pbradley/blast/')

    f1 = '/scratch/Phil/junk1_al2_'+str(time())+'.seq'
    f2 = '/scratch/Phil/junk2_al2_'+str(time())+'.seq'
    count = 1
    while exists(f1) or exists(f2):
        f1 = f1[:-1 * len(str(count-1))] + str(count)
        f2 = f1[:-1 * len(str(count-1))] + str(count)
        count = count + 1
        
    out = open(f1,'w')
    out.write('>s1\n'+seq1+'\n')
    out.close()
    out = open(f2,'w')
    out.write('>s2\n'+seq2+'\n')
    out.close()
    
    lines = popen('./bl2seq -F F -i '+f1+' -j '+f2+' -p blastp').readlines()
    alignments  = []
    a1 = ''
    a2 = ''
    start1 = 1000000
    start2 = 1000000

    for line in lines:
        if line[0:10] == ' Identitie' and a1:
            al = {}
            s1 = string.join(string.split(a1,'-'),'')
            s2 = string.join(string.split(a2,'-'),'')
            off1 = string.find(seq1,s1)
            off2 = string.find(seq2,s2)
            if off1 == start1-1 and off2 == start2-1:
                for i in range(len(a1)):
                    if a1[i] != '-' and a2[i] != '-':
                        pos1 = i + off1 - string.count(a1[:i],'-')
                        pos2 = i + off2 - string.count(a2[:i],'-')
                        al[pos1] = pos2
                alignments.append(al)
            a1 = ''
            a2 = ''
            start1 = 1000000
            start2 = 1000000

        elif line[0:6] == 'Query:':
            l = string.split(line)
            assert len(l) == 4
            a1 = a1 + l[2]
            start1 = min(start1,int(l[1]))
        elif line[0:6] == 'Sbjct:':
            l = string.split(line)
            assert len(l) == 4
            a2 = a2 + l[2]
            start2 = min(start2,int(l[1]))
            
    al = {}
    s1 = string.join(string.split(a1,'-'),'')
    s2 = string.join(string.split(a2,'-'),'')
    off1 = string.find(seq1,s1)
    off2 = string.find(seq2,s2)
    if off1 == start1-1 and off2 == start2-1:
        for i in range(len(a1)):
            if a1[i] != '-' and a2[i] != '-':
                pos1 = i + off1 - string.count(a1[:i],'-')
                pos2 = i + off2 - string.count(a2[:i],'-')
                al[pos1] = pos2
        alignments.append(al)

    remove(f1)
    remove(f2)
    chdir('/users/pbradley/python/')
    return alignments

    
def Align(a,b):
    ## first use subword matching, length 6
    a2b = {}
    b2a = {}
    if string.count(b,a) == 1:
        pos = string.find(b,a)
        for i in range(len(a)):
            a2b[i] = pos+i
        return a2b
    elif string.count(a,b):
        pos = string.find(a,b)
        for i in range(len(b)):
            a2b[i+pos] = i
        return a2b
    
    W = 6
    done = 0
    while not done:
        for i in range(len(a)-W+1):
            word = a[i:i+W]
            if string.count(b,word) == 1:
                pos = string.find(b,word)
                for j in range(W):
                    p1 = i+j
                    p2 = pos+j
                    if p1 not in a2b.keys() and p2 not in b2a.keys():
                        a2b[p1] = p2
                        b2a[p2] = p1

        for i in range(len(b)-W+1):
            word = b[i:i+W]
            if string.count(a,word) == 1:
                pos = string.find(a,word)
                for j in range(W):
                    p1 = pos+j
                    p2 = i+j
                    if p1 not in a2b.keys() and p2 not in b2a.keys():
                        a2b[p1] = p2
                        b2a[p2] = p1

        ks = a2b.keys()
        ks.sort()
        done = 1
        for i in range(len(ks)-1):
            if a2b[ks[i+1]] < a2b[ks[i]]:
                a2b = {}
                b2a = {}
                done = 0
                W = W + 1
                break

    if len(ks)==min(len(a),len(b)):
        return a2b

    prev = -1
    before = {}
    for i in range(len(a)):
        if i in a2b.keys():
            prev = a2b[i]
        else:
            before[i] = prev

    after = {}
    prev = len(b)
    for j in range(len(a)):
        i = len(a)-j-1
        if i in a2b.keys():
            prev = a2b[i]
        else:
            after[i] = prev



    ## now use the bl2seq output to fill in gaps
    alignments = allal2seq(a,b)
    for al in alignments:
        #print 'a2b',a2b
        #print 'al',al
        #print 'b4',before
        #print 'after',after
        
        for i in al.keys():
            if i not in a2b.keys() and before[i] < al[i] < after[i]:
                a2b[i] = al[i]

        #reset the before and after
        prev = -1
        before = {}
        for i in range(len(a)):
            if i in a2b.keys():
                prev = a2b[i]
            else:
                before[i] = prev

        after = {}
        prev = len(b)
        for j in range(len(a)):
            i = len(a)-j-1
            if i in a2b.keys():
                prev = a2b[i]
            else:
                after[i] = prev

        ks = a2b.keys()
        ks.sort()
        for i in range(len(ks)-1):
            assert a2b[ks[i+1]]>a2b[ks[i]]
        
            
    ks = a2b.keys()
    ks.sort()

    for i in range(len(ks)-1):
        assert a2b[ks[i+1]]>a2b[ks[i]]
    return a2b

def NoBlastAlign(a,b):
    ## first use subword matching, length 6
    a2b = {}
    b2a = {}
    if string.count(b,a) == 1:
        pos = string.find(b,a)
        for i in range(len(a)):
            a2b[i] = pos+i
        return a2b
    elif string.count(a,b):
        pos = string.find(a,b)
        for i in range(len(b)):
            a2b[i+pos] = i
        return a2b
    
    W = 6
    done = 0
    while not done:
        for i in range(len(a)-W+1):
            word = a[i:i+W]
            if string.count(b,word) == 1:
                pos = string.find(b,word)
                for j in range(W):
                    p1 = i+j
                    p2 = pos+j
                    if p1 not in a2b.keys() and p2 not in b2a.keys():
                        a2b[p1] = p2
                        b2a[p2] = p1

        ks = a2b.keys()
        ks.sort()
        done = 1
        for i in range(len(ks)-1):
            if a2b[ks[i+1]] < a2b[ks[i]]:
                a2b = {}
                b2a = {}
                done = 0
                W = W + 1
                break

    return a2b

def NBAlign(a,b):
    LB = len(b)
    LA = len(a)
    
    a2b = {}
    if string.count(b,a) == 1:
        pos = string.find(b,a)
        for i in range(LA):
            a2b[i] = pos+i
        return a2b
    elif string.count(a,b):
        pos = string.find(a,b)
        for i in range(LB):
            a2b[i+pos] = i
        return a2b

    W = 6
    counts = {}
    for i in range(LA-W+1):
        word = a[i:i+W]
        start = 0
        while string.count(b[start:],word):
            pos = string.find(b[start:],word)+start
            offset = pos-i
            if not counts.has_key(offset):counts[offset] = []
            for j in range(W):
                p1 = i+j
                counts[offset].append(p1)
            start = pos+1

    l = []
    for k in counts.keys():l.append([len(counts[k]),k])
    l.sort()
    l.reverse()
    

    before = {}
    for i in range(LA):before[i] = -1
    after = {}
    for i in range(LA):after[i] = LB

    for pair in l:
        #print pair
        offset = pair[1]
        for p1 in counts[offset]:
            if not a2b.has_key(p1) and \
               before[p1]<p1+offset and \
               after[p1]>p1+offset:
                a2b[p1] = p1+offset
    
        #reset before,after
        prev = -1
        before = {}
        for i in range(LA):
            if i in a2b.keys():
                prev = a2b[i]
            else:
                before[i] = prev

        after = {}
        prev = LB
        for j in range(LA):
            i = len(a)-j-1
            if i in a2b.keys():
                prev = a2b[i]
            else:
                after[i] = prev
    
    if len(a2b.keys()) == min(LA,LB):
        return a2b

    ## fill in gaps
    prev = -1
    before = {}
    for i in range(LA):
        if i in a2b.keys():
            prev = i ## different from above!!!!!!!!!!!
        else:
            before[i] = prev

    after = {}
    prev = LB
    for j in range(LA):
        i = len(a)-j-1
        if i in a2b.keys():
            prev = i ## different from above!!!!!!!!!!!!!!
        else:
            after[i] = prev

    for i in range(LA):
        if not a2b.has_key(i) and before[i]>=0 and after[i]<LB:
            o1 = a2b[before[i]] - before[i]
            o2 = a2b[after[i]] - after[i]
            if o1==o2:
                a2b[i] = i+o1
    return a2b
