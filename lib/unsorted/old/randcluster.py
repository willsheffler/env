#/usr/bin/python
import random

clusters = "gebb maat yah ptah niau gebb maat yah ptah niau bes apep".strip().split()

jobs = """condor/b2_1b72.condor  condor/b2_1hz6.condor  condor/b2_1r69.condor
condor/b2_1dcj.condor  condor/b2_1mky.condor  condor/b2_2reb.condor
condor/b2_1di2.condor  condor/b2_1n0u.condor  condor/b2_2tif.condor
condor/b2_1dtj.condor  condor/b2_1ogw.condor
condor/n2_1b72.condor  condor/n2_1hz6.condor  condor/n2_1r69.condor
condor/n2_1dcj.condor  condor/n2_1mky.condor  condor/n2_2reb.condor
condor/n2_1di2.condor  condor/n2_1n0u.condor  condor/n2_2tif.condor
condor/n2_1dtj.condor  condor/n2_1ogw.condor
""".strip().split()

for ii in range(len(jobs)):
    c = clusters[ii%len(clusters)]
    j = jobs[ii]
    print 'ssh %(c)s -n " cd boinc_relax; condor_submit %(j)s " &'%vars()
