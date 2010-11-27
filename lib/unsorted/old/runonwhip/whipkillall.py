import sys,os

killstring = sys.argv[1]

for ii in (1,2,3,4,5,6,7,8,9,10,11):
    if str(ii) in sys.argv:
        continue
    pscmd = "ssh whip%02i 'ps aux' | grep sheffler | egrep whipkillall | grep %s"
    ps = os.popen(pscmd%(ii,killstring)).read().strip() 
    if ps != "":
        print "ps '"+ps+"'"
        ps = [ x.split()[1] for x in ps.split('\n') ]
        for p in ps:
            os.system("ssh whip%02i 'killall -9 %s'"%(ii,p))
    
