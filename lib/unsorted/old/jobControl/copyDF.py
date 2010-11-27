import sys,os,time
from DecoyData import decoyProts
from batchjobs import FeaturesJob



for dset in decoyProts:
    for prot in decoyProts[dset]:
        protdir = '/users/sheffler/data/decoys/%s/%s'%(dset,prot[:4])
        d = os.listdir(protdir)
        d = [r for r in d if r[:5]=='rand_' and not '.' in r]        
        if len(d): rand = d[0]
        else:         rand = "rand"
        for group in ('gsgr','gsbr','relax_native','all',rand):
            #print "checking feature files"
            done = False
            while not done:
                try:
                    l = os.listdir(protdir)
                    done = True
                except OSError:
                    print 'failed to list dir', protdir
            if group in os.listdir(protdir):
                d = os.listdir(protdir+'/'+group)
                pdb = [f for f in d if f[ -4:]=='.pdb']
                df  = [f for f in d if f[-14:]=='.decoyfeatures']
                dc  = [f for f in d if f[-11:]=='.decoycheat']
                #if len(df) != len(pdb):
                if '-fix' in sys.argv:
                    if len(df) < len(pdb) or len(dc) < len(pdb):
                        print "%20s"%dset,prot,"%30s"%group,len(pdb),len(df),len(dc)
                        j = FeaturesJob(dset,prot,group)
                        j.run()
                        while not j.done():
                            print j.status()
                            time.sleep(10)
                else:
                    print "%20s"%dset,prot,"%30s"%group,len(pdb),len(df),len(dc)
                done = False
                if '-copy' in sys.argv:
                    print "copying files"
                    while not done:
                        try:
                            d = "/users/sheffler/data/decoys/data_set_3/"+dset+'/'+prot
                            os.system("mkdir -p "+d)
                            #os.system("cp "+protdir+'/'+group+' '+d)
                            os.system("cd "+protdir+"; tar -c "+group+" > "+d+'/'+group+'.tar')
                            os.system("gzip "+d+'/'+group+'.tar')
                            done = True
                        except OSError:
                            print 'failed, trying again',dset,prot,group
                    
