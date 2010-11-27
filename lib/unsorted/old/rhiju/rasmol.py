def Rclean(rin,rout):
    rin.write('select none\n')
    rin.flush()
    line = rout.readline()
    while line and (len(line)<10 or line[:10] != 'No atoms s'):
        #print line[:-1]
        line = rout.readline()
    return
