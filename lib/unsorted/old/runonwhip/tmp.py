import sys,os,time


outfiles = """
/work/baker/boinc_analysis/3-7/BARCODE_301scj_top500.out
/work/baker/boinc_analysis/3-7/BARCODE_301shf_top500.out
/work/baker/boinc_analysis/3-7/BARCODE_301ten_top500.out
/work/baker/boinc_analysis/3-7/BARCODE_301tig_top500.out
/work/baker/boinc_analysis/3-7/BARCODE_301ubi_top500.out
/work/baker/boinc_analysis/3-7/BARCODE_30256b_top500.out
/work/baker/boinc_analysis/3-7/BARCODE_302chf_top500.out
/work/baker/boinc_analysis/3-7/BARCODE_302ci2_top500.out
/work/baker/boinc_analysis/3-7/BARCODE_304ubp_top500.out
/work/baker/boinc_analysis/3-7/BARCODE_305cro_top500.out
/work/baker/boinc_analysis/2-22/HBLR_1.01b72_top500.out
/work/baker/boinc_analysis/2-22/HBLR_1.01dcj_top500.out
/work/baker/boinc_analysis/2-22/HBLR_1.01di2_top500.out
/work/baker/boinc_analysis/2-22/HBLR_1.01dtj_top500.out
/work/baker/boinc_analysis/2-22/HBLR_1.01hz6_top500.out
/work/baker/boinc_analysis/2-22/HBLR_1.01mky_top500.out
/work/baker/boinc_analysis/2-22/HBLR_1.01n0u_top500.out
/work/baker/boinc_analysis/2-22/HBLR_1.01ogw_top500.out
/work/baker/boinc_analysis/2-22/HBLR_1.01r69_top500.out
/work/baker/boinc_analysis/2-22/HBLR_1.02reb_top500.out
/work/baker/boinc_analysis/2-22/HBLR_1.02tif_top500.out    
""".strip().split()

nstruct = {}

for o in outfiles:
    log = o[1:].replace('/','-')[:-4]+"_sasa.log"
    nstruct[log] = int(os.popen("grep SCORE %s | wc -l"%(o)).read())-1

#print nstruct

while True:
    count = 0
    for log in nstruct:
        if log in os.listdir('log'):
            done = int(os.popen("grep 'output basic features' log/%(log)s | wc -l"%vars()).read())
            if done == nstruct[log]:
                count += 1
            else:
                print "%0.3f"%(float(done)/nstruct[log]),log
    total = len(nstruct)
    print "%(count)i of %(total)i jobs done"%vars()
    time.sleep(10)
