import sys,os,time

where = "/work/sheffler"
if len(sys.argv) > 1:
    where = sys.argv[1]

while True:
    print "checking",where
    d = os.popen("/bin/ls "+where).read()
    #print d.find("No such file or directory")
    if -1 == d.find("No such file or directory"):
        break
    time.sleep(10)
    
