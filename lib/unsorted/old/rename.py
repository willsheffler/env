import sys,os

for i in sys.argv[1:]:
    if i.endswith('.pdb.pdb'):
        print ('mv %s %s'%( i, i[:-4] ) )
        os.system('mv %s %s'%( i, i[:-4] ) )
