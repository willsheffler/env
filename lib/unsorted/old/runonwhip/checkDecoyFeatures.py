import sys,os

for d,dirs,files in os.walk("/users/sheffler/data/decoys"):
    if     'gsbr' in d or \
           'gsgr' in d or \
           'all'  in d or \
           'relax_native' in d
           'rand_' in d:
        print dir
