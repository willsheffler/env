#!/usr/bin/python

## we assume that the long dimension is closer to filling a page than the short
##  dimension

import string
from os import system,popen,chdir
from os.path import exists
from glob import glob
from whrandom import random
import sys
from math import floor


def Help():
    print
    print '-'*75
    print 'Usage:\n\n',sys.argv[0],'<composite-plot> <N> <input-plot1> {<input-plot2> ... }\n'
    print 'for N by N plotting'
    print '-'*75
    print

    assert 0==1

if len(sys.argv) == 1:
    Help()

args = sys.argv[1:]

by_width = 1

if '-l' in args:
    by_width = 0
    pos = args.index('-l')
    del args[pos]
if '-w' in args:
    by_width = 1
    pos = args.index('-w')
    del args[pos]

out_file = args[0]
N = int(args[1])
plots = args[2:]

if exists(out_file):
    print
    print '!!!!!!!!!!!!!!!!!!!!',out_file,'already exists!!!!!!!!!!!!\n'
#    system('rm '+out_file);
#    if out_file in plots:
#       plots.remove(out_file)
#       print plots
    Help()

print
print '-'*75
print 'Combining',string.join(plots,'\n          '),'\n\n...into',out_file
print '-'*75
print

def MakePage(out,old_plots,L,W,size,by_width):
    plots = old_plots+['']*(L*W-len(old_plots))

    out.write('\\begin{tabular}{'+'c'*W+'}\n')
    for i in range(L):
        for j in range(W):
            if plots[W*i+j]:
                if by_width:
                    out.write('\\epsfxsize='+size+'\n\\epsfbox{'+plots[W*i+j]+'}\n')
                else:
                    out.write('\\epsfysize='+size+'\n\\epsfbox{'+plots[W*i+j]+'}\n')

            if (j+1)%W: out.write('&\n')

        out.write('\\\\\n')
    out.write('\\end{tabular}\n')
    return


base_name = 'junk_combine_plots_'+str(random())
out = open(base_name+'.tex','w')
out.write('\\documentclass{article}\n')

out.write('\\usepackage{epsf}\n')


out.write('\\textwidth=8.5in\n')
out.write('\\textheight=11in\n')
out.write('\\topmargin=-0.5in\n')
out.write('\\oddsidemargin=-1in\n')
out.write('\\evensidemargin=-1in\n')
out.write('\\pagestyle{empty}\n')
out.write('\\begin{document}\n')


PLOTS_PER_PAGE = N*N
W = N
L = N

if not by_width:
    SIZE = str( 0.1*(10*10/L))+'in'
else:
    SIZE = str( 0.1*(floor(70.0/L)))+'in'


for i in range((len(plots)-1)/PLOTS_PER_PAGE +1):

    MakePage(out,plots[i*PLOTS_PER_PAGE:i*PLOTS_PER_PAGE+PLOTS_PER_PAGE],
             L,W,SIZE,by_width)
    out.write('\n\n\n')

out.write('\\end{document}\n')
out.close()

system('latex '+base_name+'.tex 2> /dev/null > /dev/null')
assert exists(base_name+'.dvi')

system('dvips -o '+out_file+' '+base_name+'.dvi 2> /dev/null > /dev/null')
system('rm '+base_name+'.dvi')
system('rm '+base_name+'.tex')
system('rm '+base_name+'.log')
system('rm '+base_name+'.aux')


