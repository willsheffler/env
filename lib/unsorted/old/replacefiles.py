#!/usr/bin/python

import os,sys,re

def replacefiles(files,regexp,rep):
    regexp = '('+regexp+')'
    print regexp
    print rep
    for f in files:
        print f
        fin = open(f)
        s = fin.read()
        fin.close()
        while True:
            #regexp = "Arg"
            print 'matching "'+regexp+'"'
            m = re.compile(regexp).search(s)
            if m is None:
                break
            else:
                tmp = rep
                g = m.groups()
                print 'matches',g
                for ii in range(1,len(g)):
                    print tmp,'\\'+`ii`,g[ii]
                    tmp = tmp.replace('\\'+`ii`,g[ii])
                s = s[:m.start()]+tmp+s[m.end():]
                print s
        #fout = open(f,'w')
        #fout.write(s)
        #fout.close()
        
replacefiles(sys.argv[3:],sys.argv[1],sys.argv[2])
