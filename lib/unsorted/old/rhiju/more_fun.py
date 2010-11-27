#!/usr/bin/python
import string
from popen2 import popen2
from sys import argv,exit


rasmol_exe = '~/bin/rasmol_32BIT'

if len(argv)<3:
    print 'usage: %s <rasmol command file> <decoy list file>'%(argv[0])
    exit()

#### some routines for managing the rasmol pipe #########
def clean_pipe(rin,rout):
    rin.flush()
    rin.write('GARBAGE GARBAGE GARBAGE\n')
    rin.flush()
    line = rout.readline()
    while not string.count(line,'GARBAGE'):
        if Good_line(line):
            print 'RasMol says:',line[:-1]
        line = rout.readline()
        
def Good_line(line): ## the lines that get output ---
    return not string.count(line,'GARBAGE') and \
           not string.count(line,'^') and \
           not string.count(line,'Unrecognised') and \
           not string.count(line,'RasMol>') and \
           not string.count(line,'atoms selected')
#########################################################


## read the commands:
print 'Reading the rasmol command file:', argv[1]
command_lines = open(argv[1],'r').readlines()


## now look at the decoys
rout,rin = popen2(rasmol_exe)

print 'Reading the decoy list file:',argv[2]
decoys = map(lambda x:string.split(x)[0],
             open(argv[2],'r').readlines())

print 'found %s decoys'%(len(decoys))

for decoy in decoys:
    print 'loading:',decoy
    rin.write('zap\nload %s\n'%decoy)
    clean_pipe(rin,rout)
    rin.writelines(command_lines)

    rin.flush()

    print 'RasMol> ',
    ans = raw_input()

    while ans != 'q' and ans != 'quit':
        clean_pipe(rin,rout)

        print 'you typed: (%s)'%ans

        rin.write(ans+'\n')
        rin.flush()
    
        ans = raw_input()
        
