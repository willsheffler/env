#!/users/sheffler/software/bin/python
import sys,os,random,time
from subprocess import Popen,PIPE
from DecoyData import decoyProts
from PathsFile import makePathsFile

class Job:

    def __init__(me,command,out,err,number,name,append=False,timeout=1000):
        me.command = command
        me.process = None
        me.start = None
        me.end   = None
        me.host = None
        me.out  = out
        me.err  = err
        me.result = None
        me.append = append
        me.number = number
        me.name = name
        me.ID = "%04d "%number+name
        me.prejob()
        me.timeout = timeout
    def prejob(me):
        pass
    def postjob(me):
        pass
    def status(me):
        if me.process is None:
            return "job '"+me.ID+"' not started"%me
        elif me.done():            
            return "job '"+me.ID+"' done after %d seconds"%(me.end-me.start)
        else:
            return "job '"+me.ID+"' running for %d seconds"%(time.time()-me.start)
    def started(me):
        if me.process is None:
            return False
        else:
            return True
    def done(me):
        if me.process is None:
            return False
        if time.time() - me.start > me.timeout:
            me.end = time.time()
            me.result = None
            return True
        if me.process.poll() is not None:
            me.end = time.time()
            me.result = me.process.poll()
            me.postjob()
            return True
        else:
            return False
    def run(me,host="localhost"):
        me.host = host
        if me.append:
            me.out = open(me.out,'a')
            me.err = open(me.err,'a')
        else:
            me.out = open(me.out,'w')
            me.err = open(me.err,'w')
        me.out.write("#################################################################\n"+
                     "# STDOUT from job started by python on "+time.ctime()+"\n# Command: "+me.command+
                     "\n# Host: "+host+"\n"+
                     "#################################################################\n")
        me.err.write("#################################################################\n"+
                     "# STDERR from job started by python on "+time.ctime()+"\n# Command: "+me.command+
                     "\n# Host: "+host+"\n"+
                     "#################################################################\n")
        me.out.flush()
        me.err.flush()
        if host == 'localhost':
            torun = '%s'%(me.command)
        else:
            torun = 'ssh %s "%s"'%(host,me.command)
        me.process = Popen(torun,shell=True,stdout=me.out,stderr=me.err,close_fds=True)
        me.start   = time.time()
        return True
    def halt(me,verbose=True):
        if me.start is None:
            if verbose: print "can't halt: not started"
        elif me.end is None:
            os.system('ssh %(me.host)s kill %i'%me.process.pid)
            if verbose: print "killed pid %i"%me.process.pid
        else:
            if verbose: print "can't halt: already finished"
    def __getitem__(me,item):
        return getattr(me,item)
            

class RegionJob(Job):
    def __init__(me,dataset,prot,group,number=0):
        from JobTemplates import RegionJobTemplate
        me.dataset = dataset
        me.prot    = prot
        me.group   = group
        me.series   = 'rg'
        me.protein = prot[:4]
        me.chain   = prot[4]        
        me.location = '/users/sheffler/data/decoys/%(dataset)s/%(protein)s/'%me
        me.paths    = makePathsFile(location=me.location)
        command = RegionJobTemplate%me
        out = '%(location)s/log/%(group)s_region.log'%me
        err = '%(location)s/log/%(group)s_region.err'%me
        name = prot+'.'+group
        Job.__init__(me,command,out,err,number=number,name=name)
    def prejob(me):
        os.system("""mkdir -p %(location)s/log ;
        mkdir -p %(location)s/data ;
        if [ -e %(location)s/data/%(group)s_region.fasc ]; then
        rm -f %(location)s/log/tmp_score_file.fasc ;
        fi;
        if [ -e %(location)s/data/%(group)s_region.data ]; then
        rm -f %(location)s/data/%(group)s_region.data ;
        fi;
        """%me)
    def postjob(me):
        os.system("nice -n14 grep DF_REGION %(location)s/log/%(group)s_region.log \
        > %(location)s/data/%(group)s_region.data "%me)
    def status(me,npdbs=-1):
        if not me.started() or me.done():
            return Job.status(me)
        me.out.flush()
        me.err.flush()
        out = me.out.name
        err = me.err.name
        ndone = int(os.popen2("grep 'FEATURES output features region for' %s | wc -l"%out)[1].read())
        nerr  = int(os.popen2("egrep -v '^#' %s | wc -l "%err)[1].read())
        if(npdbs>0):
            f = ndone/npdbs
            return "job '"+me.ID+"' time: %8.2f  host: %s  errors: %2d  progress: %1.2f"%(
                time.time()-me.start,me.host,nerr,f)
        else:
            return "job '"+me.ID+"' time: %8.2f  host: %s  errors: %2d  progress: %5d"%(
                time.time()-me.start,me.host,nerr,ndone)



class FeaturesJob(Job):
    def __init__(me,dataset,prot,group,number=0):
        from JobTemplates import FeaturesJobTemplate
        me.dataset = dataset
        me.prot    = prot
        me.group   = group
        me.series   = 'ft'
        me.protein = prot[:4]
        me.chain   = prot[4]        
        me.location = '/users/sheffler/data/decoys/%(dataset)s/%(protein)s/'%me
        pdbout = me.location + '/' + me.group
        me.paths    = makePathsFile({'pdbout':pdbout},location=me.location)
        command = FeaturesJobTemplate%me
        out = '%(location)s/log/%(group)s_features.log'%me
        err = '%(location)s/log/%(group)s_features.err'%me
        name = prot+'.'+group
        Job.__init__(me,command,out,err,number=number,name=name)
    def postjob(me):
        os.system("rm -f %(location)s/%(paths)s"%me)
    def status(me,npdbs=-1):
        if not me.started() or me.done():
            return Job.status(me)
        me.out.flush()
        me.err.flush()
        out = me.out.name
        err = me.err.name
        ndone = int(os.popen2("grep 'FEATURES output features features' %s | wc -l"%out)[1].read())
        nerr  = int(os.popen2("egrep -v '^#' %s | wc -l "%err)[1].read())
        if(npdbs>0):
            f = ndone/npdbs
            return "host: %s  time: %8.2f  errors: %2d  progress: %1.2f"%(
                me.host,time.time()-me.start,nerr,f)+ " job '"+me.ID+"'" 
        else:
            return "host: %s  time: %8.2f  errors: %2d  progress: %4d"%(
                me.host,time.time()-me.start,nerr,ndone)+ " job '"+me.ID+"'" 






#out = open('test.log','w')
#err = open('test.err','w')
#p = Popen(regionjob%vars(),stdout=out,stderr=err,cwd=location,shell=True)

#out,err = os.popen2(job%vars())

#from JobTemplates import regionJobTemplate
#j = Job(regionJobTemplate%vars(),'test.out','test.err')

#j = RegionJob('bqian_caspbench','t196_','gsbr')#
#print j.status()
#j.run('whip01')
#while j.result is None:
#    print j.statu#s()
#    time.sleep(1)#
#
#def test(*args):
#    print args
#test('asdf','rthr')


if 'batchjobs.py' in sys.argv[0]:
    jobs = []
    ii = 1
    for dset in decoyProts:
        for prot in decoyProts[dset]:
            groups = ['relax_native','all','gsbr','gsgr']
            rand = os.listdir('/users/sheffler/data/decoys/%s/%s/'%(dset,prot[:-1]))
            rand = [r for r in rand if r[:5]=="rand_" and '.' not in r]    
            #rms = os.listdir('/users/sheffler/data/decoys/%s/%s/'%(dset,prot[:-1]))
            #rms = [r for r in rms if r[:3]=='rms' and len(r)==4]
            if len(rand):
                groups.append(rand[0])
            #if len(rms):
            #    groups.extend(rms)
            for group in groups:
                print "Batching",dset,prot,group
                jobs.append(FeaturesJob(dset,prot,group,number=ii))
                ii += 1

    #print jobs[19].command
    #jobs = jobs[19:20]
    
    #jobs = [ FeaturesJob("bqian_caspbench","t196_","gsbr",0) ]
    whipjobs = [None]*11
    for ii in range(11):
        if ii >= len(jobs):
            break
        jobs[ii].run("whip%02d"%(ii+1))
        whipjobs[ii] = jobs[ii]
        jobs = jobs[1:]

    print whipjobs
    
    while True:
        for ii in range(11):
            whip = "whip%02d"%(ii+1)
            job = whipjobs[ii]
            if job is None:
                print whip,None
            else:
                print whip,job.status()
            if job is not None:
                if not job.started():
                    job.run(whip)
                elif job.done():
                    whipjobs[ii] = None
            elif len(jobs) > 0:
                whipjobs[ii] = jobs.pop()
        if min([j is None for j in whipjobs]) is True:
            break
        sys.stdout.flush()
        time.sleep(5)
        print '----------------------------------'
    
