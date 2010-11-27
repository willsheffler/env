import os,sys

#print sys.argv[1]

filelist = open(sys.argv[1]).readlines()
#print filelist

for ii in range(len(filelist)):
    f = filelist[ii]
    try:
        d = int(os.popen("du "+f).read().split()[0])
        if d > 4:
            print f,
    except:
        pass
    

