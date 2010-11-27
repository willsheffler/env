import sys,os

datafiles = os.popen('find ~/ -name *.rdata').read().split()
print datafiles


for d in datafiles:
    n = "/users/sheffler/data/rdata/"+"_".join(d.split('/')[3:])
    #print d
    #print n
    #print
    print 'mv '+d+' '+n
    os.system('mv '+d+' '+n)

    
