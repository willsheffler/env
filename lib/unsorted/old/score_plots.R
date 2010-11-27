source('~/scripts/R-utils.R')

scorefiles <- list()

scorefiles[['hblrrot462']] <- list()
scorefiles[['hblrrot462']][['1b72']] <- "/users/sheffler/project/score_tweaks/HBLR_1.0_1b72_ROT_TRIALS_TRIE_462_0.rand10000.sc"
scorefiles[['hblrrot462']][['1di2']] <- "/users/sheffler/project/score_tweaks/HBLR_1.0_1di2_ROT_TRIALS_TRIE_462_0.rand10000.sc"
scorefiles[['hblrrot462']][['1dtj']] <- "/users/sheffler/project/score_tweaks/HBLR_1.0_1dtj_ROT_TRIALS_TRIE_462_0.rand10000.sc"
scorefiles[['hblrrot462']][['1hz6']] <- "/users/sheffler/project/score_tweaks/HBLR_1.0_1hz6_ROT_TRIALS_TRIE_462_0.rand10000.sc"
scorefiles[['hblrrot462']][['1mky']] <- "/users/sheffler/project/score_tweaks/HBLR_1.0_1mky_ROT_TRIALS_TRIE_462_0.rand10000.sc"
scorefiles[['hblrrot462']][['1n0u']] <- "/users/sheffler/project/score_tweaks/HBLR_1.0_1n0u_ROT_TRIALS_TRIE_462_0.rand10000.sc"
scorefiles[['hblrrot462']][['1ogw']] <- "/users/sheffler/project/score_tweaks/HBLR_1.0_1ogw_ROT_TRIALS_TRIE_462_0.rand10000.sc"
scorefiles[['hblrrot462']][['2tif']] <- "/users/sheffler/project/score_tweaks/HBLR_1.0_2tif_ROT_TRIALS_TRIE_462_0.rand10000.sc"

scorefiles[['hblrsmet']] <- list()
scorefiles[['hblrsmet1']][['1b72']] <- "/work/boinc/results_ralph/HBLR_SM_/HBLR_SM_ET_1b72_458_0.sc"
scorefiles[['hblrsmet1']][['1di2']] <- "/work/boinc/results_ralph/HBLR_SM_/HBLR_SM_ET_1di2_458_0.sc"
scorefiles[['hblrsmet1']][['1dtj']] <- "/work/boinc/results_ralph/HBLR_SM_/HBLR_SM_ET_1dtj_458_0.sc"
scorefiles[['hblrsmet1']][['1hz6']] <- "/work/boinc/results_ralph/HBLR_SM_/HBLR_SM_ET_1hz6_458_0.sc"
scorefiles[['hblrsmet1']][['1mky']] <- "/work/boinc/results_ralph/HBLR_SM_/HBLR_SM_ET_1mky_458_0.sc"
scorefiles[['hblrsmet1']][['1n0u']] <- "/work/boinc/results_ralph/HBLR_SM_/HBLR_SM_ET_1n0u_458_0.sc"
scorefiles[['hblrsmet1']][['1ogw']] <- "/work/boinc/results_ralph/HBLR_SM_/HBLR_SM_ET_1ogw_458_0.sc"
scorefiles[['hblrsmet1']][['2tif']] <- "/work/boinc/results_ralph/HBLR_SM_/HBLR_SM_ET_2tif_458_0.sc"

scorefiles[['hblrsmet2']] <- list()
scorefiles[['hblrsmet2']][['1b72']] <- "/work/boinc/results_ralph/HBLR_SM_/HBLR_SM_ET_W1_1b72_466_0.sc"
scorefiles[['hblrsmet2']][['1di2']] <- "/work/boinc/results_ralph/HBLR_SM_/HBLR_SM_ET_W1_1di2_466_0.sc"
scorefiles[['hblrsmet2']][['1dtj']] <- "/work/boinc/results_ralph/HBLR_SM_/HBLR_SM_ET_W1_1dtj_466_0.sc"
scorefiles[['hblrsmet2']][['1hz6']] <- "/work/boinc/results_ralph/HBLR_SM_/HBLR_SM_ET_W1_1hz6_466_0.sc"
scorefiles[['hblrsmet2']][['1mky']] <- "/work/boinc/results_ralph/HBLR_SM_/HBLR_SM_ET_W1_1mky_466_0.sc"
scorefiles[['hblrsmet2']][['1n0u']] <- "/work/boinc/results_ralph/HBLR_SM_/HBLR_SM_ET_W1_1n0u_466_0.sc"
scorefiles[['hblrsmet2']][['1ogw']] <- "/work/boinc/results_ralph/HBLR_SM_/HBLR_SM_ET_W1_1ogw_466_0.sc"
scorefiles[['hblrsmet2']][['2tif']] <- "/work/boinc/results_ralph/HBLR_SM_/HBLR_SM_ET_W1_2tif_466_0.sc"


scorefiles <- list()
scorefiles[['smrelax1']] <- list()
scorefiles[['smrelax1']][['1b72']] <- "/users/sheffler/data/clusters/boinc_relax/results/dd1b72.out"
scorefiles[['smrelax1']][['1di2']] <- "/users/sheffler/data/clusters/boinc_relax/results/dd1di2.out"
scorefiles[['smrelax1']][['1dtj']] <- "/users/sheffler/data/clusters/boinc_relax/results/dd1dtj.out"
scorefiles[['smrelax1']][['1hz6']] <- "/users/sheffler/data/clusters/boinc_relax/results/dd1hz6.out"
scorefiles[['smrelax1']][['1mky']] <- "/users/sheffler/data/clusters/boinc_relax/results/dd1mky.out"
scorefiles[['smrelax1']][['1n0u']] <- "/users/sheffler/data/clusters/boinc_relax/results/dd1n0u.out"
scorefiles[['smrelax1']][['1ogw']] <- "/users/sheffler/data/clusters/boinc_relax/results/dd1ogw.out"
scorefiles[['smrelax1']][['2tif']] <- "/users/sheffler/data/clusters/boinc_relax/results/dd2tif.out"

scorefiles[['boinctop1']] <- list()
scorefiles[['boinctop1']][['1b72']] <- "/users/sheffler/data/clusters/boinc_relax/results/HBLR_1.01b72_rottrials_top2000.out"
scorefiles[['boinctop1']][['1di2']] <- "/users/sheffler/data/clusters/boinc_relax/results/HBLR_1.01di2_rottrials_top2000.out"
scorefiles[['boinctop1']][['1dtj']] <- "/users/sheffler/data/clusters/boinc_relax/results/HBLR_1.01dtj_rottrials_top2000.out" 
scorefiles[['boinctop1']][['1hz6']] <- "/users/sheffler/data/clusters/boinc_relax/results/HBLR_1.01hz6_rottrials_top2000.out" 
scorefiles[['boinctop1']][['1mky']] <- "/users/sheffler/data/clusters/boinc_relax/results/HBLR_1.01mky_rottrials_top2000.out" 
scorefiles[['boinctop1']][['1n0u']] <- "/users/sheffler/data/clusters/boinc_relax/results/HBLR_1.01n0u_rottrials_top2000.out"  
scorefiles[['boinctop1']][['1ogw']] <- "/users/sheffler/data/clusters/boinc_relax/results/HBLR_1.01ogw_rottrials_top2000.out" 
scorefiles[['boinctop1']][['2tif']] <- "/users/sheffler/data/clusters/boinc_relax/results/HBLR_1.02tif_rottrials_top2000.out"



scorefiles <- list()
scorefiles[['t4rlxsm']] <- list()
scorefiles[['t4rlxsm']][['1b72']] <- "/users/sheffler/data/clusters/boinc_relax/results/ee1b72.out"
scorefiles[['t4rlxsm']][['1di2']] <- "/users/sheffler/data/clusters/boinc_relax/results/ee1di2.out"
scorefiles[['t4rlxsm']][['1dtj']] <- "/users/sheffler/data/clusters/boinc_relax/results/ee1dtj.out"
#scorefiles[['t4rlxsm']][['1hz6']] <- "/users/sheffler/data/clusters/boinc_relax/results/ee1hz6.out"
#scorefiles[['t4rlxsm']][['1mky']] <- "/users/sheffler/data/clusters/boinc_relax/results/ee1mky.out"
scorefiles[['t4rlxsm']][['1n0u']] <- "/users/sheffler/data/clusters/boinc_relax/results/ee1n0u.out"
scorefiles[['t4rlxsm']][['1ogw']] <- "/users/sheffler/data/clusters/boinc_relax/results/ee1ogw.out"

scorefiles[['t4rlxsm']][['1dcj']] <- "/users/sheffler/data/clusters/boinc_relax/results/gg1dcj.out"
scorefiles[['t4rlxsm']][['1r69']] <- "/users/sheffler/data/clusters/boinc_relax/results/gg1r69.out"
scorefiles[['t4rlxsm']][['2reb']] <- "/users/sheffler/data/clusters/boinc_relax/results/gg2reb.out"

scorefiles[['t4rlxstd']] <- list()
scorefiles[['t4rlxstd']][['1b72']] <- "/users/sheffler/data/clusters/boinc_relax/results/ff1b72.out"
scorefiles[['t4rlxstd']][['1di2']] <- "/users/sheffler/data/clusters/boinc_relax/results/ff1di2.out"
scorefiles[['t4rlxstd']][['1dtj']] <- "/users/sheffler/data/clusters/boinc_relax/results/ff1dtj.out"
#scorefiles[['t4rlxstd']][['1hz6']] <- "/users/sheffler/data/clusters/boinc_relax/results/ff1hz6.out"
#scorefiles[['t4rlxstd']][['1mky']] <- "/users/sheffler/data/clusters/boinc_relax/results/ff1mky.out"
scorefiles[['t4rlxstd']][['1n0u']] <- "/users/sheffler/data/clusters/boinc_relax/results/ff1n0u.out"
scorefiles[['t4rlxstd']][['1ogw']] <- "/users/sheffler/data/clusters/boinc_relax/results/ff1ogw.out"

scorefiles[['t4rlxstd']][['1dcj']] <- "/users/sheffler/data/clusters/boinc_relax/results/hh1dcj.out"
scorefiles[['t4rlxstd']][['1r69']] <- "/users/sheffler/data/clusters/boinc_relax/results/hh1r69.out"
scorefiles[['t4rlxstd']][['2reb']] <- "/users/sheffler/data/clusters/boinc_relax/results/hh2reb.out"



possufiles <- list()
possufiles[['top400_relax_std']] <- list()
possufiles[['top400_relax_std']][['1b72']] <- "/users/sheffler/data/clusters/boinc_relax/results/ra1b72.out"
possufiles[['top400_relax_std']][['1di2']] <- "/users/sheffler/data/clusters/boinc_relax/results/ra1di2.out"
possufiles[['top400_relax_std']][['1dtj']] <- "/users/sheffler/data/clusters/boinc_relax/results/ra1dtj.out"
possufiles[['top400_relax_std']][['1hz6']] <- "/users/sheffler/data/clusters/boinc_relax/results/ra1hz6.out"
possufiles[['top400_relax_std']][['1mky']] <- "/users/sheffler/data/clusters/boinc_relax/results/ra1mky.out"
possufiles[['top400_relax_std']][['1n0u']] <- "/users/sheffler/data/clusters/boinc_relax/results/ra1n0u.out"
possufiles[['top400_relax_std']][['1ogw']] <- "/users/sheffler/data/clusters/boinc_relax/results/ra1ogw.out"
possufiles[['top400_relax_std']][['2tif']] <- "/users/sheffler/data/clusters/boinc_relax/results/ra2tif.out"
possufiles[['top400_relax_std']][['1dcj']] <- "/users/sheffler/data/clusters/boinc_relax/results/hh1dcj.out"
possufiles[['top400_relax_std']][['1r69']] <- "/users/sheffler/data/clusters/boinc_relax/results/hh1r69.out"
possufiles[['top400_relax_std']][['2reb']] <- "/users/sheffler/data/clusters/boinc_relax/results/hh2reb.out"

possufiles[['t4rlxsm']] <- list()
possufiles[['t4rlxsm']][['1b72']] <- "/users/sheffler/data/clusters/boinc_relax/results/ee1b72.out"
possufiles[['t4rlxsm']][['1di2']] <- "/users/sheffler/data/clusters/boinc_relax/results/ee1di2.out"
possufiles[['t4rlxsm']][['1dtj']] <- "/users/sheffler/data/clusters/boinc_relax/results/ee1dtj.out"
#possufiles[['t4rlxsm']][['1hz6']] <- "/users/sheffler/data/clusters/boinc_relax/results/ee1hz6.out"
#possufiles[['t4rlxsm']][['1mky']] <- "/users/sheffler/data/clusters/boinc_relax/results/ee1mky.out"
possufiles[['t4rlxsm']][['1n0u']] <- "/users/sheffler/data/clusters/boinc_relax/results/ee1n0u.out"
possufiles[['t4rlxsm']][['1ogw']] <- "/users/sheffler/data/clusters/boinc_relax/results/ee1ogw.out"

possufiles[['t4rlxsm']][['1dcj']] <- "/users/sheffler/data/clusters/boinc_relax/results/gg1dcj.out"
possufiles[['t4rlxsm']][['1r69']] <- "/users/sheffler/data/clusters/boinc_relax/results/gg1r69.out"
possufiles[['t4rlxsm']][['2reb']] <- "/users/sheffler/data/clusters/boinc_relax/results/gg2reb.out"


possufiles[['top400_relax_sm_possu']] <- list()
possufiles[['top400_relax_sm_possu']][['1b72']] <- "/users/sheffler/data/clusters/boinc_relax/results/rb1b72.out"
possufiles[['top400_relax_sm_possu']][['1di2']] <- "/users/sheffler/data/clusters/boinc_relax/results/rb1di2.out"
possufiles[['top400_relax_sm_possu']][['1dtj']] <- "/users/sheffler/data/clusters/boinc_relax/results/rb1dtj.out"
possufiles[['top400_relax_sm_possu']][['1hz6']] <- "/users/sheffler/data/clusters/boinc_relax/results/rb1hz6.out"
possufiles[['top400_relax_sm_possu']][['1mky']] <- "/users/sheffler/data/clusters/boinc_relax/results/rb1mky.out"
possufiles[['top400_relax_sm_possu']][['1n0u']] <- "/users/sheffler/data/clusters/boinc_relax/results/rb1n0u.out"
possufiles[['top400_relax_sm_possu']][['1ogw']] <- "/users/sheffler/data/clusters/boinc_relax/results/rb1ogw.out"
possufiles[['top400_relax_sm_possu']][['2tif']] <- "/users/sheffler/data/clusters/boinc_relax/results/rb2tif.out"
possufiles[['top400_relax_sm_possu']][['1dcj']] <- "/users/sheffler/data/clusters/boinc_relax/results/ii1dcj.out"
possufiles[['top400_relax_sm_possu']][['1r69']] <- "/users/sheffler/data/clusters/boinc_relax/results/ii1r69.out"
possufiles[['top400_relax_sm_possu']][['2reb']] <- "/users/sheffler/data/clusters/boinc_relax/results/ii2reb.out"

#possufiles[['top400_relax_sm_possu_sameatrrep']] <- list()
#possufiles[['top400_relax_sm_possu_sameatrrep']][['1b72']] <- "/users/sheffler/data/clusters/boinc_relax/results/rc1b72.out"
#possufiles[['top400_relax_sm_possu_sameatrrep']][['1di2']] <- "/users/sheffler/data/clusters/boinc_relax/results/rc1di2.out"
#possufiles[['top400_relax_sm_possu_sameatrrep']][['1dtj']] <- "/users/sheffler/data/clusters/boinc_relax/results/rc1dtj.out"
#possufiles[['top400_relax_sm_possu_sameatrrep']][['1hz6']] <- "/users/sheffler/data/clusters/boinc_relax/results/rc1hz6.out"
#possufiles[['top400_relax_sm_possu_sameatrrep']][['1mky']] <- "/users/sheffler/data/clusters/boinc_relax/results/rc1mky.out"
#possufiles[['top400_relax_sm_possu_sameatrrep']][['1n0u']] <- "/users/sheffler/data/clusters/boinc_relax/results/rc1n0u.out"
#possufiles[['top400_relax_sm_possu_sameatrrep']][['1ogw']] <- "/users/sheffler/data/clusters/boinc_relax/results/rc1ogw.out"
#possufiles[['top400_relax_sm_possu_sameatrrep']][['2tif']] <- "/users/sheffler/data/clusters/boinc_relax/results/rc2tif.out"


abscfiles <- list()

abscfiles[['hblrrot462']] <- list()
abscfiles[['hblrrot462']][['1b72']] <- "/users/sheffler/project/score_tweaks/HBLR_1.0_1b72_ROT_TRIALS_TRIE_462_0.rand10000.sc"
abscfiles[['hblrrot462']][['1di2']] <- "/users/sheffler/project/score_tweaks/HBLR_1.0_1di2_ROT_TRIALS_TRIE_462_0.rand10000.sc"
abscfiles[['hblrrot462']][['1dtj']] <- "/users/sheffler/project/score_tweaks/HBLR_1.0_1dtj_ROT_TRIALS_TRIE_462_0.rand10000.sc"
abscfiles[['hblrrot462']][['1hz6']] <- "/users/sheffler/project/score_tweaks/HBLR_1.0_1hz6_ROT_TRIALS_TRIE_462_0.rand10000.sc"
abscfiles[['hblrrot462']][['1mky']] <- "/users/sheffler/project/score_tweaks/HBLR_1.0_1mky_ROT_TRIALS_TRIE_462_0.rand10000.sc"
abscfiles[['hblrrot462']][['1n0u']] <- "/users/sheffler/project/score_tweaks/HBLR_1.0_1n0u_ROT_TRIALS_TRIE_462_0.rand10000.sc"
abscfiles[['hblrrot462']][['1ogw']] <- "/users/sheffler/project/score_tweaks/HBLR_1.0_1ogw_ROT_TRIALS_TRIE_462_0.rand10000.sc"
abscfiles[['hblrrot462']][['2tif']] <- "/users/sheffler/project/score_tweaks/HBLR_1.0_2tif_ROT_TRIALS_TRIE_462_0.rand10000.sc"

abscfiles[['hblrrot462']][['1dcj']] <- "/users/sheffler/project/score_tweaks/HBLR_1.01dcj_rand1000.sc"
abscfiles[['hblrrot462']][['1r69']] <- "/users/sheffler/project/score_tweaks/HBLR_1.01r69_rand1000.sc"
abscfiles[['hblrrot462']][['2reb']] <- "/users/sheffler/project/score_tweaks/HBLR_1.02reb_rand1000.sc"


abscfiles[['abrlx_sol1.6_ralph']] <- list()
abscfiles[['abrlx_sol1.6_ralph']][['1b72']] <- "/users/sheffler/project/score_tweaks/HBLR_SM_ET_W2_1b72_525_0.sc"
abscfiles[['abrlx_sol1.6_ralph']][['1di2']] <- "/users/sheffler/project/score_tweaks/HBLR_SM_ET_W2_1di2_525_0.sc"
abscfiles[['abrlx_sol1.6_ralph']][['1dtj']] <- "/users/sheffler/project/score_tweaks/HBLR_SM_ET_W2_1dtj_525_0.sc"
abscfiles[['abrlx_sol1.6_ralph']][['1hz6']] <- "/users/sheffler/project/score_tweaks/HBLR_SM_ET_W2_1hz6_525_0.sc"
abscfiles[['abrlx_sol1.6_ralph']][['1mky']] <- "/users/sheffler/project/score_tweaks/HBLR_SM_ET_W2_1mky_525_0.sc"
abscfiles[['abrlx_sol1.6_ralph']][['1n0u']] <- "/users/sheffler/project/score_tweaks/HBLR_SM_ET_W2_1n0u_525_0.sc"
abscfiles[['abrlx_sol1.6_ralph']][['1ogw']] <- "/users/sheffler/project/score_tweaks/HBLR_SM_ET_W2_1ogw_525_0.sc"
abscfiles[['abrlx_sol1.6_ralph']][['2tif']] <- "/users/sheffler/project/score_tweaks/HBLR_SM_ET_W2_2tif_525_0.sc"

abscfiles[['abrlx_sol1.6_ralph']][['1dcj']] <- "/users/sheffler/project/score_tweaks/HBLR_SM_ET_W2_1dcj_562_0.out"
abscfiles[['abrlx_sol1.6_ralph']][['1r69']] <- "/users/sheffler/project/score_tweaks/HBLR_SM_ET_W2_1r69_562_0.out"
abscfiles[['abrlx_sol1.6_ralph']][['2reb']] <- "/users/sheffler/project/score_tweaks/HBLR_SM_ET_W2_2reb_562_0.out"



read.score.files <- function( s=list() , scorefiles ) {
  for(sn in names(scorefiles) ) {
    for(prot in names(scorefiles[[sn]]) ) {
      if( ! prot %in% names(s) )
        s[[prot]] <- list()
      print(paste('reading scores',sn,prot))
      scfile <- scorefiles[[sn]][[prot]]      
      #print("grep'ing")
      system(paste("grep SCORE",scfile,"> /tmp/tmp.sc"))
      #print("cleaning")
      system("python ~/python/clean_outfile.py /tmp/tmp.sc /tmp/tmp2.sc")      
      #print("reading")
      s[[prot]][[sn]] <- read.table("/tmp/tmp2.sc",header=T,row.names=NULL)
      print(dim(s[[prot]][[sn]]))
      print(sum(is.na( s[[prot]][[sn]] )))
      #print("done")
    }
  }
  return(s)
}

#sc <- read.score.files(list(),scorefiles)


plot.scores <- function(s,
                        xval="rms",yval="score",
                        Qmn=0.0,
                        Qmx=0.9,
                        dsetnames=NULL,
                        pch='.',col=1:7,
                        Sxmn=999999,Sxmx=-999999,Symn=999999,Symx=-9999999,
                        NORM=F,
                        NBEST=1,
                        BSRNUM=99999,
                        ...) {
  par(mfrow=hw(length(s)+2))
  par(mar=c(2,2,2,1))
  par(mgp=c(1,0,0))
  par(tck=0)
  if(is.null(dsetnames))
    dsetnames <- names(s[[1]])
  bestscorerms <-list()
  for(prot in names(s)) {
    print(paste('starting loop1',prot))
    bestscorerms[[prot]] <- list()
    smean <- list()
    ssd   <- list()    
    for(dset in names(s[[prot]])){
      bestscorerms[[prot]][[dset]] <- list()
      if(NORM) {
        tmp <- s[[prot]][[dset]][[yval]]
        tmp <- tmp[ tmp > mean(tmp) - 10*sd(tmp) ]
        smean[[dset]] <- min(tmp)
        ssd  [[dset]] <- median(tmp)-min(tmp)
      } else {
        smean[[dset]] <- 0
        ssd  [[dset]] <- 1
      }
    }
    xmn <- Sxmn
    xmx <- Sxmx
        names(tmp) <- names(bestscorerms)
    ymn <- Symn
    ymx <- Symx
    for(ii in 1:length(dsetnames)) {
      print(paste('starting loop2',prot,dset))
      dset <- dsetnames[[ii]]
      tmp <- s[[prot]][[dset]][ 1:min(BSRNUM,dim(s[[prot]][[dset]])[1] ), ]
      print(dim(tmp))
      bestscorerms[[prot]][[dset]] <- mean(
              tmp$rms[ order( tmp$score )[1:NBEST] ]  )
      print(paste(prot,dset,tmp$rms[ order( tmp$score )[1:NBEST] ] ))
      xmn <- min(xmn,quantile(s[[prot]][[dset]][[xval]],Qmn))
      print('a')
      xmx <- max(xmx,quantile(s[[prot]][[dset]][[xval]],Qmx))
      print('b')
      tmp <- s[[prot]][[dset]][[yval]]
      print('c')
      tmp <- tmp[ tmp > mean(tmp) - 10*sd(tmp) ]
      print('d')
      ymn <- min(ymn,quantile((tmp-smean[[dset]])/ssd[[dset]],Qmn))
      print('e')
      ymx <- max(ymx,quantile((tmp-smean[[dset]])/ssd[[dset]],Qmx))
    }
    print('f')
    #browser()
    plot(1,type='n',
         xlim=c(xmn,xmx),
         ylim=c(ymn,ymx),
         main=prot,xlab=xval,ylab=yval)
    print('f2')
    for(ii in 1:length(dsetnames)) {
      print('g')
      dset <- dsetnames[[ii]]
      print(dset)
      points(s[[prot]][[dset]][[xval]],
             (s[[prot]][[dset]][[yval]]-smean[[dset]])/ssd[[dset]],
             pch=pch,col=col[ii],...)
      print('h')
    }
  }

  #browser()
  tmp <- list()
  for(n in names(bestscorerms)) {
    tmp[[n]] <- list()
    for(d in 1:length(dsetnames)) {
      tmp[[n]][[d]] <- bestscorerms[[n]][[dsetnames[d]]]
    }
  }
  browser()
  par(mar=c(6,1,1,0))
  barplot(as.matrix(as.data.frame(lapply( tmp ,as.numeric))),beside=T,
          legend.text=names(tmp[[1]]),col=1:length(dsetnames) )

  par(cex=1)
  plot(0:1,0:1,'n',main="legend")
  legend(0.0,1,dsetnames,col=col[1:length(dsetnames)],lwd=3)

  return(bestscorerms)
}



if(F) {

s <- read.score.files(list(),scorefiles)

plot.scores(s)






}
