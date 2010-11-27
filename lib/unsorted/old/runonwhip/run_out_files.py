import sys,os,random,time
#import pathsfile


dfjobs = ("distances","sasa")
if len(sys.argv) > 1:
    dfjobs = sys.argv[1].split()


#print outfiles
#print outnames

if len(sys.argv) > 2:
    outfiles = sys.argv[2:]

#print outfiles
#sys.exit()

outfiles = [ o.replace("/Users","/work") for o in outfiles ]
outnames = [ o.replace('/','-')[1:-4] for o in outfiles ]
workdir = "/scratch/sheffler/run_on_whip/" # os.popen('pwd').read().strip().replace("/Users",'/work')
#paths   = "/users/sheffler/paths/abspaths.txt"
paths   = "/work/sheffler/paths/abspaths.txt"
binary  = "/work/sheffler/rosetta_decoyfeatures.gcc"
#options_template = "-score -new_reader -try_both_his_tautomers -fa_input -fa_output -decoyfeatures %(dfjob)s -df_file df/%(dffile)s -paths %(paths)s "#-scorefile %(scorefile)s"
options_template = "-score -new_reader -try_both_his_tautomers -fa_input -fa_output -df_%(dfjob)s -paths %(paths)s "
cmd_template = "ssh -nx whip%(whip)02i 'cd %(workdir)s ; nice -n100 %(binary)s %(options)s %(input_options)s > log/%(logfile)s' &"
whips   = [1,2,3,4,5,6,7,8,9,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26]
whipjob = {}
for whip in whips:
    whipjob[whip] = "GARBAGE GARBAGE GARBAGE"

#pathsbasedir = ""
#paths = makepathsfile(location='./',sequence=)


for dfjob in dfjobs:
  for ii in range(len(outfiles)):
    inputfile = outfiles[ii]
    outname = outnames[ii]
    dffile  = outname+".df"
    logfile = outname+"_"+dfjob+".log"
    scorefile = outname+"_"+dfjob
    #paths   = "paths.txt"
    options = options_template%vars()
    started = False
    while not started:
#      random.shuffle(whips)
      for whip in whips:
          topcmd = "ssh whip%02i -nx 'top -bn1 | head -n10 | tail -n3'"%whip
          top = os.popen(topcmd).read()
          if top.count("No route to host") is 0:
              print 'checking whip%02i'%whip,
              load = float(top.split("\n")[1].split()[8])                    
              print "has load:",load,
              ps = os.popen("ssh -x whip%02i 'ps aux' | grep %s"%(whip,binary)).read().count(binary)/2
              print "running for sheffler:",ps
              if load < 20.0:
                  #print ps
                  if ps < 1:# or (random.random()<0.5 and ps < 2):
                      if inputfile.endswith('.out'):
                          input_options =  "-silent_input -all -s "+inputfile
                      else:
                          input_options = "-l "+inputfile
                      # print cmd_template%vars()
                      os.system(cmd_template%vars())
                      whipjob[whip] = inputfile
                      print "started job",inputfile,"on whip%2i"%whip
                      started = True
          time.sleep(1)
          if started:
              break
                    
## outfiles = """
## /work/baker/boinc_analysis/3-7/BARCODE_301a19_top500.out
## /work/baker/boinc_analysis/3-7/BARCODE_301a32_top500.out
## /work/baker/boinc_analysis/3-7/BARCODE_301a68_top200.out
## /work/baker/boinc_analysis/3-7/BARCODE_301acf_top500.out
## /work/baker/boinc_analysis/3-7/BARCODE_301aiu_top200.out
## /work/baker/boinc_analysis/3-7/BARCODE_301b3a_top200.out
## /work/baker/boinc_analysis/3-7/BARCODE_301bk2_top500.out
## /work/baker/boinc_analysis/3-7/BARCODE_301bm8_top200.out
## /work/baker/boinc_analysis/3-7/BARCODE_301bq9_top500.out
## /work/baker/boinc_analysis/3-7/BARCODE_301c8c_top500.out
## /work/baker/boinc_analysis/3-7/BARCODE_301c9o_top200.out
## /work/baker/boinc_analysis/3-7/BARCODE_301cc8_top500.out
## /work/baker/boinc_analysis/3-7/BARCODE_301ctf_top500.out
## /work/baker/boinc_analysis/3-7/BARCODE_301elw_top500.out
## /work/baker/boinc_analysis/3-7/BARCODE_301enh_top500.out
## /work/baker/boinc_analysis/3-7/BARCODE_301fna_top500.out
## /work/baker/boinc_analysis/3-7/BARCODE_301ig5_top500.out
## /work/baker/boinc_analysis/3-7/BARCODE_301iib_top500.out
## /work/baker/boinc_analysis/3-7/BARCODE_301opd_top500.out
## /work/baker/boinc_analysis/3-7/BARCODE_301pgx_top500.out
## /work/baker/boinc_analysis/3-7/BARCODE_301scj_top500.out
## /work/baker/boinc_analysis/3-7/BARCODE_301shf_top500.out
## /work/baker/boinc_analysis/3-7/BARCODE_301ten_top500.out
## /work/baker/boinc_analysis/3-7/BARCODE_301tig_top500.out
## /work/baker/boinc_analysis/3-7/BARCODE_301ubi_top500.out
## /work/baker/boinc_analysis/3-7/BARCODE_30256b_top500.out
## /work/baker/boinc_analysis/3-7/BARCODE_302chf_top500.out
## /work/baker/boinc_analysis/3-7/BARCODE_302ci2_top500.out
## /work/baker/boinc_analysis/3-7/BARCODE_304ubp_top500.out
## /work/baker/boinc_analysis/3-7/BARCODE_305cro_top500.out
## /work/baker/boinc_analysis/2-22/HBLR_1.01b72_top500.out
## /work/baker/boinc_analysis/2-22/HBLR_1.01dcj_top500.out
## /work/baker/boinc_analysis/2-22/HBLR_1.01di2_top500.out
## /work/baker/boinc_analysis/2-22/HBLR_1.01dtj_top500.out
## /work/baker/boinc_analysis/2-22/HBLR_1.01hz6_top500.out
## /work/baker/boinc_analysis/2-22/HBLR_1.01mky_top500.out
## /work/baker/boinc_analysis/2-22/HBLR_1.01n0u_top500.out
## /work/baker/boinc_analysis/2-22/HBLR_1.01ogw_top500.out
## /work/baker/boinc_analysis/2-22/HBLR_1.01r69_top500.out
## /work/baker/boinc_analysis/2-22/HBLR_1.02reb_top500.out
## /work/baker/boinc_analysis/2-22/HBLR_1.02tif_top500.out
##    """.strip().split() 



# /work/boinc/results/BARCODE_/BARCODE_30_1a19A_NATIVE_313_0.first.out
# /work/boinc/results/BARCODE_/BARCODE_30_1a32__NATIVE_313_0.first.out
# /work/boinc/results/BARCODE_/BARCODE_30_1a68__NATIVE_313_0.first.out
# /work/boinc/results/BARCODE_/BARCODE_30_1acf__NATIVE_313_0.first.out
# /work/boinc/results/BARCODE_/BARCODE_30_1aiu__NATIVE_313_0.first.out
# /work/boinc/results/BARCODE_/BARCODE_30_1b3aA_NATIVE_313_0.first.out
# /work/boinc/results/BARCODE_/BARCODE_30_1bk2__NATIVE_313_0.first.out
# /work/boinc/results/BARCODE_/BARCODE_30_1bm8__NATIVE_313_0.first.out
# /work/boinc/results/BARCODE_/BARCODE_30_1bq9A_NATIVE_313_0.first.out
# /work/boinc/results/BARCODE_/BARCODE_30_1c8cA_NATIVE_313_0.first.out
# /work/boinc/results/BARCODE_/BARCODE_30_1c9oA_NATIVE_313_0.first.out
# /work/boinc/results/BARCODE_/BARCODE_30_1cc8A_NATIVE_313_0.first.out
# /work/boinc/results/BARCODE_/BARCODE_30_1ctf__NATIVE_313_0.first.out
# /work/boinc/results/BARCODE_/BARCODE_30_1elwA_NATIVE_313_0.first.out
# /work/boinc/results/BARCODE_/BARCODE_30_1enh__NATIVE_313_0.first.out
# /work/boinc/results/BARCODE_/BARCODE_30_1fna__NATIVE_313_0.first.out
# /work/boinc/results/BARCODE_/BARCODE_30_1ig5A_NATIVE_313_0.first.out
# /work/boinc/results/BARCODE_/BARCODE_30_1iibA_NATIVE_313_0.first.out
# /work/boinc/results/BARCODE_/BARCODE_30_1opd__NATIVE_313_0.first.out
# /work/boinc/results/BARCODE_/BARCODE_30_1pgx__NATIVE_313_0.first.out
# /work/boinc/results/BARCODE_/BARCODE_30_1scjB_NATIVE_313_0.first.out
# /work/boinc/results/BARCODE_/BARCODE_30_1shfA_NATIVE_313_0.first.out
# /work/boinc/results/BARCODE_/BARCODE_30_1ten__NATIVE_313_0.first.out
# /work/boinc/results/BARCODE_/BARCODE_30_1tig__NATIVE_313_0.first.out
# /work/boinc/results/BARCODE_/BARCODE_30_1ubi__NATIVE_313_0.first.out
# /work/boinc/results/BARCODE_/BARCODE_30_256bA_NATIVE_313_0.first.out
# /work/boinc/results/BARCODE_/BARCODE_30_2chf__NATIVE_313_0.first.out
# /work/boinc/results/BARCODE_/BARCODE_30_2ci2I_NATIVE_313_0.first.out
# /work/boinc/results/BARCODE_/BARCODE_30_4ubpA_NATIVE_313_0.first.out
# /work/boinc/results/BARCODE_/BARCODE_30_5croA_NATIVE_313_0.first.out
# /work/boinc/results/HBLR_NAT/HBLR_NATIVE_1b72_335_0.first.out
# /work/boinc/results/HBLR_NAT/HBLR_NATIVE_1dcj_335_0.first.out
# /work/boinc/results/HBLR_NAT/HBLR_NATIVE_1di2_335_0.first.out
# /work/boinc/results/HBLR_NAT/HBLR_NATIVE_1dtj_335_0.first.out
# /work/boinc/results/HBLR_NAT/HBLR_NATIVE_1hz6_335_0.first.out
# /work/boinc/results/HBLR_NAT/HBLR_NATIVE_1mky_335_0.first.out
# /work/boinc/results/HBLR_NAT/HBLR_NATIVE_1n0u_335_0.first.out
# /work/boinc/results/HBLR_NAT/HBLR_NATIVE_1ogw_335_0.first.out
# /work/boinc/results/HBLR_NAT/HBLR_NATIVE_1r69_335_0.first.out
# /work/boinc/results/HBLR_NAT/HBLR_NATIVE_2reb_335_0.first.out
# /work/boinc/results/HBLR_NAT/HBLR_NATIVE_2tif_335_0.first.out





###/work/baker/boinc_analysis/3-7/BARCODE_301a19_rand1000.out
###/work/baker/boinc_analysis/3-7/BARCODE_301a32_rand1000.out
###/work/baker/boinc_analysis/3-7/BARCODE_301a68_rand1000.out
###/work/baker/boinc_analysis/3-7/BARCODE_301acf_rand1000.out
###/work/baker/boinc_analysis/3-7/BARCODE_301aiu_rand1000.out
###/work/baker/boinc_analysis/3-7/BARCODE_301b3a_rand1000.out
###/work/baker/boinc_analysis/3-7/BARCODE_301bk2_rand1000.out
###/work/baker/boinc_analysis/3-7/BARCODE_301bm8_rand1000.out

##outfiles += """
##/work/baker/boinc_analysis/3-7/BARCODE_301bq9_rand1000.out
##/work/baker/boinc_analysis/3-7/BARCODE_301c8c_rand1000.out
##/work/baker/boinc_analysis/3-7/BARCODE_301c9o_rand1000.out
##/work/baker/boinc_analysis/3-7/BARCODE_301cc8_rand1000.out
##/work/baker/boinc_analysis/3-7/BARCODE_301ctf_rand1000.out
##/work/baker/boinc_analysis/3-7/BARCODE_301elw_rand1000.out
##/work/baker/boinc_analysis/3-7/BARCODE_301enh_rand1000.out
##/work/baker/boinc_analysis/3-7/BARCODE_301fna_rand1000.out
##/work/baker/boinc_analysis/3-7/BARCODE_301ig5_rand1000.out
##/work/baker/boinc_analysis/3-7/BARCODE_301iib_rand1000.out
##/work/baker/boinc_analysis/3-7/BARCODE_301opd_rand1000.out
##/work/baker/boinc_analysis/3-7/BARCODE_301pgx_rand1000.out
##/work/baker/boinc_analysis/3-7/BARCODE_301scj_rand1000.out
##/work/baker/boinc_analysis/3-7/BARCODE_301shf_rand1000.out
##/work/baker/boinc_analysis/3-7/BARCODE_301ten_rand1000.out
##/work/baker/boinc_analysis/3-7/BARCODE_301tig_rand1000.out
##/work/baker/boinc_analysis/3-7/BARCODE_301ubi_rand1000.out
##/work/baker/boinc_analysis/3-7/BARCODE_30256b_rand1000.out
##/work/baker/boinc_analysis/3-7/BARCODE_302chf_rand1000.out
##/work/baker/boinc_analysis/3-7/BARCODE_302ci2_rand1000.out
##/work/baker/boinc_analysis/3-7/BARCODE_304ubp_rand1000.out
##/work/baker/boinc_analysis/3-7/BARCODE_305cro_rand1000.out
##/work/baker/boinc_analysis/2-22/HBLR_1.01b72_rand1000.out
##/work/baker/boinc_analysis/2-22/HBLR_1.01dcj_rand1000.out
##/work/baker/boinc_analysis/2-22/HBLR_1.01di2_rand1000.out
##/work/baker/boinc_analysis/2-22/HBLR_1.01dtj_rand1000.out
##/work/baker/boinc_analysis/2-22/HBLR_1.01hz6_rand1000.out
##/work/baker/boinc_analysis/2-22/HBLR_1.01mky_rand1000.out
##/work/baker/boinc_analysis/2-22/HBLR_1.01n0u_rand1000.out
##/work/baker/boinc_analysis/2-22/HBLR_1.01ogw_rand1000.out
##/work/baker/boinc_analysis/2-22/HBLR_1.01r69_rand1000.out
##/work/baker/boinc_analysis/2-22/HBLR_1.02reb_rand1000.out
##/work/baker/boinc_analysis/2-22/HBLR_1.02tif_rand1000.out
##""".strip().split()

##hhfix_outfiles = """
##/work/sheffler/data/decoys/score_tweaks/h-h-repulsion/target22301_1a19A.out
##/work/sheffler/data/decoys/score_tweaks/h-h-repulsion/target22302_1a32_.out
##/work/sheffler/data/decoys/score_tweaks/h-h-repulsion/target22303_1a68_.out
##/work/sheffler/data/decoys/score_tweaks/h-h-repulsion/target22304_1acf_.out
##/work/sheffler/data/decoys/score_tweaks/h-h-repulsion/target22305_1ail_.out
##/work/sheffler/data/decoys/score_tweaks/h-h-repulsion/target22306_1aiu_.out
##/work/sheffler/data/decoys/score_tweaks/h-h-repulsion/target22307_1b3aA.out
##/work/sheffler/data/decoys/score_tweaks/h-h-repulsion/target22308_1bgf_.out
##/work/sheffler/data/decoys/score_tweaks/h-h-repulsion/target22309_1bk2_.out
##/work/sheffler/data/decoys/score_tweaks/h-h-repulsion/target22310_1bkrA.out
##/work/sheffler/data/decoys/score_tweaks/h-h-repulsion/target22311_1bm8_.out
##/work/sheffler/data/decoys/score_tweaks/h-h-repulsion/target22312_1bq9A.out
##/work/sheffler/data/decoys/score_tweaks/h-h-repulsion/target22313_1c8cA.out
##/work/sheffler/data/decoys/score_tweaks/h-h-repulsion/target22314_1c9oA.out
##/work/sheffler/data/decoys/score_tweaks/h-h-repulsion/target22315_1cc8A.out
##/work/sheffler/data/decoys/score_tweaks/h-h-repulsion/target22316_1cei_.out
##/work/sheffler/data/decoys/score_tweaks/h-h-repulsion/target22317_1cg5B.out
##/work/sheffler/data/decoys/score_tweaks/h-h-repulsion/target22318_1ctf_.out
##/work/sheffler/data/decoys/score_tweaks/h-h-repulsion/target22319_1dhn_.out
##/work/sheffler/data/decoys/score_tweaks/h-h-repulsion/target22320_1e6iA.out
##/work/sheffler/data/decoys/score_tweaks/h-h-repulsion/target22321_1elwA.out
##/work/sheffler/data/decoys/score_tweaks/h-h-repulsion/target22322_1enh_.out
##/work/sheffler/data/decoys/score_tweaks/h-h-repulsion/target22323_1ew4A.out
##/work/sheffler/data/decoys/score_tweaks/h-h-repulsion/target22324_1eyvA.out
##/work/sheffler/data/decoys/score_tweaks/h-h-repulsion/target22325_1fkb_.out
##/work/sheffler/data/decoys/score_tweaks/h-h-repulsion/target22326_1fna_.out
##/work/sheffler/data/decoys/score_tweaks/h-h-repulsion/target22327_1gvp_.out
##/work/sheffler/data/decoys/score_tweaks/h-h-repulsion/target22328_1hz6A.out
##/work/sheffler/data/decoys/score_tweaks/h-h-repulsion/target22329_1ig5A.out
##/work/sheffler/data/decoys/score_tweaks/h-h-repulsion/target22330_1iibA.out
##/work/sheffler/data/decoys/score_tweaks/h-h-repulsion/target22331_1kpeA.out
##/work/sheffler/data/decoys/score_tweaks/h-h-repulsion/target22332_1lis_.out
##/work/sheffler/data/decoys/score_tweaks/h-h-repulsion/target22333_1louA.out
##/work/sheffler/data/decoys/score_tweaks/h-h-repulsion/target22334_1npsA.out
##/work/sheffler/data/decoys/score_tweaks/h-h-repulsion/target22335_1opd_.out
##/work/sheffler/data/decoys/score_tweaks/h-h-repulsion/target22336_1pgx_.out
##/work/sheffler/data/decoys/score_tweaks/h-h-repulsion/target22337_1ptq_.out
##/work/sheffler/data/decoys/score_tweaks/h-h-repulsion/target22338_1r69_.out
##/work/sheffler/data/decoys/score_tweaks/h-h-repulsion/target22339_1rnbA.out
##/work/sheffler/data/decoys/score_tweaks/h-h-repulsion/target22340_1scjB.out
##/work/sheffler/data/decoys/score_tweaks/h-h-repulsion/target22341_1shfA.out
##/work/sheffler/data/decoys/score_tweaks/h-h-repulsion/target22342_1ten_.out
##/work/sheffler/data/decoys/score_tweaks/h-h-repulsion/target22343_1tif_.out
##/work/sheffler/data/decoys/score_tweaks/h-h-repulsion/target22344_1tig_.out
##/work/sheffler/data/decoys/score_tweaks/h-h-repulsion/target22345_1tit_.out
##/work/sheffler/data/decoys/score_tweaks/h-h-repulsion/target22346_1tul_.out
##/work/sheffler/data/decoys/score_tweaks/h-h-repulsion/target22347_1ubi_.out
##/work/sheffler/data/decoys/score_tweaks/h-h-repulsion/target22348_1ughI.out
##/work/sheffler/data/decoys/score_tweaks/h-h-repulsion/target22350_1utg_.out
##/work/sheffler/data/decoys/score_tweaks/h-h-repulsion/target22351_1vcc_.out
##/work/sheffler/data/decoys/score_tweaks/h-h-repulsion/target22352_1vie_.out
##/work/sheffler/data/decoys/score_tweaks/h-h-repulsion/target22353_1vls_.out
##/work/sheffler/data/decoys/score_tweaks/h-h-repulsion/target22354_1who_.out
##/work/sheffler/data/decoys/score_tweaks/h-h-repulsion/target22355_1wit_.out
##/work/sheffler/data/decoys/score_tweaks/h-h-repulsion/target22356_256bA.out
##/work/sheffler/data/decoys/score_tweaks/h-h-repulsion/target22357_2acy_.out
##/work/sheffler/data/decoys/score_tweaks/h-h-repulsion/target22358_2chf_.out
##/work/sheffler/data/decoys/score_tweaks/h-h-repulsion/target22359_2ci2I.out
##/work/sheffler/data/decoys/score_tweaks/h-h-repulsion/target22360_2vik_.out
##/work/sheffler/data/decoys/score_tweaks/h-h-repulsion/target22361_4ubpA.out
##/work/sheffler/data/decoys/score_tweaks/h-h-repulsion/target22362_5croA.out
##""".strip().split()


hbfix_outfiles = """
/work/sheffler/data/decoys/score_tweaks/patrick-hb-fix/target22363_1a19A.out
/work/sheffler/data/decoys/score_tweaks/patrick-hb-fix/target22364_1a32_.out
/work/sheffler/data/decoys/score_tweaks/patrick-hb-fix/target22365_1a68_.out
/work/sheffler/data/decoys/score_tweaks/patrick-hb-fix/target22366_1acf_.out
/work/sheffler/data/decoys/score_tweaks/patrick-hb-fix/target22367_1ail_.out
/work/sheffler/data/decoys/score_tweaks/patrick-hb-fix/target22368_1aiu_.out
/work/sheffler/data/decoys/score_tweaks/patrick-hb-fix/target22369_1b3aA.out
/work/sheffler/data/decoys/score_tweaks/patrick-hb-fix/target22370_1bgf_.out
/work/sheffler/data/decoys/score_tweaks/patrick-hb-fix/target22371_1bk2_.out
/work/sheffler/data/decoys/score_tweaks/patrick-hb-fix/target22372_1bkrA.out
/work/sheffler/data/decoys/score_tweaks/patrick-hb-fix/target22373_1bm8_.out
/work/sheffler/data/decoys/score_tweaks/patrick-hb-fix/target22374_1bq9A.out
/work/sheffler/data/decoys/score_tweaks/patrick-hb-fix/target22375_1c8cA.out
/work/sheffler/data/decoys/score_tweaks/patrick-hb-fix/target22376_1c9oA.out
/work/sheffler/data/decoys/score_tweaks/patrick-hb-fix/target22377_1cc8A.out
/work/sheffler/data/decoys/score_tweaks/patrick-hb-fix/target22378_1cei_.out
/work/sheffler/data/decoys/score_tweaks/patrick-hb-fix/target22379_1cg5B.out
/work/sheffler/data/decoys/score_tweaks/patrick-hb-fix/target22380_1ctf_.out
/work/sheffler/data/decoys/score_tweaks/patrick-hb-fix/target22381_1dhn_.out
/work/sheffler/data/decoys/score_tweaks/patrick-hb-fix/target22382_1e6iA.out
/work/sheffler/data/decoys/score_tweaks/patrick-hb-fix/target22383_1elwA.out
/work/sheffler/data/decoys/score_tweaks/patrick-hb-fix/target22384_1enh_.out
/work/sheffler/data/decoys/score_tweaks/patrick-hb-fix/target22385_1ew4A.out
/work/sheffler/data/decoys/score_tweaks/patrick-hb-fix/target22386_1eyvA.out
/work/sheffler/data/decoys/score_tweaks/patrick-hb-fix/target22387_1fkb_.out
/work/sheffler/data/decoys/score_tweaks/patrick-hb-fix/target22388_1fna_.out
/work/sheffler/data/decoys/score_tweaks/patrick-hb-fix/target22389_1gvp_.out
/work/sheffler/data/decoys/score_tweaks/patrick-hb-fix/target22390_1hz6A.out
/work/sheffler/data/decoys/score_tweaks/patrick-hb-fix/target22391_1ig5A.out
/work/sheffler/data/decoys/score_tweaks/patrick-hb-fix/target22392_1iibA.out
/work/sheffler/data/decoys/score_tweaks/patrick-hb-fix/target22393_1kpeA.out
/work/sheffler/data/decoys/score_tweaks/patrick-hb-fix/target22394_1lis_.out
/work/sheffler/data/decoys/score_tweaks/patrick-hb-fix/target22395_1louA.out
/work/sheffler/data/decoys/score_tweaks/patrick-hb-fix/target22396_1npsA.out
/work/sheffler/data/decoys/score_tweaks/patrick-hb-fix/target22397_1opd_.out
/work/sheffler/data/decoys/score_tweaks/patrick-hb-fix/target22398_1pgx_.out
/work/sheffler/data/decoys/score_tweaks/patrick-hb-fix/target22399_1ptq_.out
/work/sheffler/data/decoys/score_tweaks/patrick-hb-fix/target22400_1r69_.out
/work/sheffler/data/decoys/score_tweaks/patrick-hb-fix/target22401_1rnbA.out
/work/sheffler/data/decoys/score_tweaks/patrick-hb-fix/target22402_1scjB.out
/work/sheffler/data/decoys/score_tweaks/patrick-hb-fix/target22403_1shfA.out
/work/sheffler/data/decoys/score_tweaks/patrick-hb-fix/target22405_1tif_.out
/work/sheffler/data/decoys/score_tweaks/patrick-hb-fix/target22406_1tig_.out
/work/sheffler/data/decoys/score_tweaks/patrick-hb-fix/target22407_1tit_.out
/work/sheffler/data/decoys/score_tweaks/patrick-hb-fix/target22408_1tul_.out
/work/sheffler/data/decoys/score_tweaks/patrick-hb-fix/target22409_1ubi_.out
/work/sheffler/data/decoys/score_tweaks/patrick-hb-fix/target22410_1ughI.out
/work/sheffler/data/decoys/score_tweaks/patrick-hb-fix/target22411_1urnA.out
/work/sheffler/data/decoys/score_tweaks/patrick-hb-fix/target22412_1utg_.out
/work/sheffler/data/decoys/score_tweaks/patrick-hb-fix/target22413_1vcc_.out
/work/sheffler/data/decoys/score_tweaks/patrick-hb-fix/target22414_1vie_.out
/work/sheffler/data/decoys/score_tweaks/patrick-hb-fix/target22415_1vls_.out
/work/sheffler/data/decoys/score_tweaks/patrick-hb-fix/target22416_1who_.out
/work/sheffler/data/decoys/score_tweaks/patrick-hb-fix/target22417_1wit_.out
/work/sheffler/data/decoys/score_tweaks/patrick-hb-fix/target22418_256bA.out
/work/sheffler/data/decoys/score_tweaks/patrick-hb-fix/target22419_2acy_.out
/work/sheffler/data/decoys/score_tweaks/patrick-hb-fix/target22420_2chf_.out
/work/sheffler/data/decoys/score_tweaks/patrick-hb-fix/target22422_2vik_.out
/work/sheffler/data/decoys/score_tweaks/patrick-hb-fix/target22423_4ubpA.out
/work/sheffler/data/decoys/score_tweaks/patrick-hb-fix/target22424_5croA.out
""".strip().split()

