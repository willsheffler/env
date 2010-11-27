import sys,os,time

where = "/work/sheffler"
if len(sys.argv) > 1:
    where = sys.argv[1]

while True:
    print "checking for "+where+"...",
    sys.stdout.flush()
    d = os.popen("/bin/ls "+where+" 2> /dev/null").read()
    #print d.find("No such file or directory")
    if "" != d:
        print "it's there! exiting."
        break
    print "sorry, no luck"
    time.sleep(10)
    
