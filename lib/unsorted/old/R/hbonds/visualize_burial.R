gvkburial <- read.table('gvkb_polar_burial.data',header=F)[,-1]
names(gvkburial) <- c('prot','pdb','aa','res','atom',
                      'hbe','hbet','dist','cosh','cosa',
                      'atomtype','sasa14','sasa10','sasa7',
                      'lk','lk2','lk4','lke','lk2e','lk4e',
                      'nb','sasafrac','anb5','anb75','anb10')
save(gvkburial,file="gvkburial.rdata")

#load('pdbpolar.rdata')
load('gvkburial.rdata')

# polarisH <- pdbpolar$atomtype==22 | pdbpolar$atomtype==25
# lkHthresh  <- quantile(pdbpolar$lk[polarisH]    , 0.5 )
# lkAthresh  <- quantile(pdbpolar$lk[!polarisH]   , 0.5 )
# lk2Hthresh  <- quantile(pdbpolar$lk2[polarisH]  , 0.5 )
# lk2Athresh  <- quantile(pdbpolar$lk2[!polarisH] , 0.5 )
# lk4Hthresh  <- quantile(pdbpolar$lk4[polarisH]  , 0.5 )
# lk4Athresh  <- quantile(pdbpolar$lk4[!polarisH] , 0.5 )
# nbthresh  <- quantile(pdbpolar$nb , 0.5 )

isH    <- gvkburial$atomtype==22 | gvkburial$atomtype==25
lkHthresh   <- quantile(gvkburial$lk[gvkisH]   , 0.5 )
lkAthresh   <- quantile(gvkburial$lk[!gvkisH]  , 0.5 )
lk2Hthresh  <- quantile(gvkburial$lk2[gvkisH]  , 0.5 )
lk2Athresh  <- quantile(gvkburial$lk2[!gvkisH] , 0.5 )
lk4Hthresh  <- quantile(gvkburial$lk4[gvkisH]  , 0.5 )
lk4Athresh  <- quantile(gvkburial$lk4[!gvkisH] , 0.5 )
lkeHthresh   <- quantile(gvkburial$lke[gvkisH]   , 0.5 )
lkeAthresh   <- quantile(gvkburial$lke[!gvkisH]  , 0.5 )
lk2eHthresh  <- quantile(gvkburial$lk2e[gvkisH]  , 0.5 )
lk2eAthresh  <- quantile(gvkburial$lk2e[!gvkisH] , 0.5 )
lk4eHthresh  <- quantile(gvkburial$lk4e[gvkisH]  , 0.5 )
lk4eAthresh  <- quantile(gvkburial$lk4e[!gvkisH] , 0.5 )
nbthresh    <- quantile(gvkburial$nb     , 0.5 )
anb5thresh  <- quantile(gvkburial$anb5   , 0.5 )
anb75thresh  <- quantile(gvkburial$anb75 , 0.5 )
anb10thresh  <- quantile(gvkburial$anb10 , 0.5 )



rasmol.whole.chain.vdw <- function(prot,file,pdbdir,res,atom) {
  # write rasmol script
  write('zap',file=file)
  write("stereo -5",file=file,append=T)
  write("set picking center",file=file,append=T)
  write(paste('load pdb "',pdbdir,prot,'.pdb"',sep=""),file=file,append=T)
  write("restrict none",file=file,append=T)
  write("select backbone",file=file,append=T)
  write("strands",file=file,append=T)
  write("select all",file=file,append=T)
  write("wireframe",file=file,append=T)
  for(ii in 1:length(res)) {
    sel <- paste(paste("(",res[ii],' and ','*.',atom[ii],')',sep=''))
    write(paste('select',sel),file=file,append=T)
    write("spacefill",file=file,append=T)
    #write(paste("select within(4.0,(",sel,"))"),file=file,append=T)
    #write("wireframe 0.2",file=file,append=T)
    #write(paste("select within(3.0,(",sel,")) and not",rosres[ii]),file=file,append=T)
    #write("wireframe 0.4",file=file,append=T)
  }
  
}

rasmol.whole.chain.vdw.comp <- function(prot,file,pdbdir,res,atom,I1,I2) {
  # write rasmol script
  write('zap',file=file)
  write("stereo -5",file=file,append=T)
  write("set picking center",file=file,append=T)
  write(paste('load pdb "',pdbdir,prot,'.pdb"',sep=""),file=file,append=T)
  write("restrict none",file=file,append=T)
  write("select backbone",file=file,append=T)
  write("strands",file=file,append=T)
  write("select all",file=file,append=T)
  write("wireframe",file=file,append=T)
  #write("spacefill 1.0",file=file,append=T)
  write("color white",file=file,append=T)
  res1not2 <- res[I1 & !I2]
  res2not1 <- res[I2 & !I1]
  resboth  <- res[I1 & I2]
  atom1not2 <- atom[I1 & !I2]
  atom2not1 <- atom[I2 & !I1]
  atomboth  <- atom[I1 & I2]
  if(length(res1not2)){
    for(ii in 1:length(res1not2)) {
      sel <- paste(paste("(",res1not2[ii],' and ','*.',atom1not2[ii],')',sep=''))
      write(paste('select',sel),file=file,append=T)
      write("spacefill",file=file,append=T)
      write("color blue",file=file,append=T)
    }
  }
  if(length(res2not1)){
    for(ii in 1:length(res2not1)) {
      sel <- paste(paste("(",res2not1[ii],' and ','*.',atom2not1[ii],')',sep=''))
      write(paste('select',sel),file=file,append=T)
      write("spacefill",file=file,append=T)
      write("color yellow",file=file,append=T)
    }
  }
  if(length(resboth)){
    for(ii in 1:length(resboth)) {
      sel <- paste(paste("(",resboth[ii],' and ','*.',atomboth[ii],')',sep=''))
      write(paste('select',sel),file=file,append=T)
      write("spacefill",file=file,append=T)
      write("color green",file=file,append=T)
    }
  }
}

find.buried.examples <- function() {
  load('gvkburial.rdata')
  prot <- "1gvkB"
  pdbdir <- "/data/sheffler/pdb/rosetta_pdb_chains/"
  filepath  <- "~/project/features/buried_polar/visualizations/"

  attach(gvkburial)
  
  file <- paste(filepath,prot,'_sasa14.rasmol',sep='')
  I <- sasa14==0
  print(sum(I))
  rasmol.whole.chain.vdw(prot,file,pdbdir, res[I] , atom[I] )
  
  file <- paste(filepath,prot,'_sasa10.rasmol',sep='')
  I <- sasa10==0
  print(sum(I))
  rasmol.whole.chain.vdw(prot,file,pdbdir, res[I] , atom[I] )
  
  file <- paste(filepath,prot,'_lk.rasmol',sep='')
  I <- (isH & lk >= lkHthresh) | (!isH & lk >= lkAthresh)
  print(sum(I))
  rasmol.whole.chain.vdw(prot,file,pdbdir, res[I] , atom[I] )
  
  file <- paste(filepath,prot,'_lk2.rasmol',sep='')
  I <- (isH & lk2 >= lk2Hthresh) | (!isH & lk2 >= lk2Athresh)
  print(sum(I))
  rasmol.whole.chain.vdw(prot,file,pdbdir, res[I] , atom[I] )
  
  file <- paste(filepath,prot,'_lk4.rasmol',sep='')
  I <- (isH & lk4 >= lk4Hthresh) | (!isH & lk4 >= lk4Athresh)
  print(sum(I))
  rasmol.whole.chain.vdw(prot,file,pdbdir, res[I] , atom[I] )

  file <- paste(filepath,prot,'_nb.rasmol',sep='')
  I <- nb >= nbthresh
  print(sum(I))
  rasmol.whole.chain.vdw(prot,file,pdbdir, res[I] , atom[I] )

  file <- paste(filepath,prot,'_sasa14_sasa10.rasmol',sep='')
  I1 <- sasa14==0
  I2 <- sasa10==0
  print(paste(sum(I1),sum(I2),sum(I1&I2)))
  rasmol.whole.chain.vdw.comp(prot,file,pdbdir,res,atom,I1,I2)
  
  file <- paste(filepath,prot,'_sasa14_lk.rasmol',sep='')
  I1 <- sasa14==0
  I2 <- (isH & lk >= lkHthresh) | (!isH & lk >= lkAthresh)
  print(paste(sum(I1),sum(I2),sum(I1&I2)))
  rasmol.whole.chain.vdw.comp(prot,file,pdbdir,res,atom,I1,I2)
  
  file <- paste(filepath,prot,'_sasa14_lk2.rasmol',sep='')
  I1 <- sasa14==0
  I2 <- (isH & lk2 >= lk2Hthresh) | (!isH & lk2 >= lk2Athresh)
  print(paste(sum(I1),sum(I2),sum(I1&I2)))
  rasmol.whole.chain.vdw.comp(prot,file,pdbdir,res,atom,I1,I2)
  
  file <- paste(filepath,prot,'_sasa14_lk4.rasmol',sep='')
  I1 <- sasa14==0
  I2 <- (isH & lk4 >= lk4Hthresh) | (!isH & lk4 >= lk4Athresh)
  print(paste(sum(I1),sum(I2),sum(I1&I2)))
  rasmol.whole.chain.vdw.comp(prot,file,pdbdir,res,atom,I1,I2)
  
  file <- paste(filepath,prot,'_lk_lk2.rasmol',sep='')
  I1 <- (isH & lk >= lkHthresh) | (!isH & lk >= lkAthresh)
  I2 <- (isH & lk2 >= lk2Hthresh) | (!isH & lk2 >= lk2Athresh)
  print(paste(sum(I1),sum(I2),sum(I1&I2)))
  rasmol.whole.chain.vdw.comp(prot,file,pdbdir,res,atom,I1,I2)
   
  file <- paste(filepath,prot,'_lk_lk4.rasmol',sep='')
  I1 <- (isH & lk >= lkHthresh) | (!isH & lk >= lkAthresh)
  I2 <- (isH & lk4 >= lk4Hthresh) | (!isH & lk4 >= lk4Athresh)
  print(paste(sum(I1),sum(I2),sum(I1&I2)))
  rasmol.whole.chain.vdw.comp(prot,file,pdbdir,res,atom,I1,I2)
  
  file <- paste(filepath,prot,'_lk2_lk4.rasmol',sep='')
  I1 <- (isH & lk2 >= lk2Hthresh) | (!isH & lk2 >= lk2Athresh)
  I2 <- (isH & lk4 >= lk4Hthresh) | (!isH & lk4 >= lk4Athresh)
  print(paste(sum(I1),sum(I2),sum(I1&I2)))
  rasmol.whole.chain.vdw.comp(prot,file,pdbdir,res,atom,I1,I2)
  
  file <- paste(filepath,prot,'_nb_sasa14.rasmol',sep='')
  I1 <- nb > nbthresh
  I2 <- sasa14==0
  print(paste(sum(I1),sum(I2),sum(I1&I2)))
  rasmol.whole.chain.vdw.comp(prot,file,pdbdir,res,atom,I1,I2)
  
  file <- paste(filepath,prot,'_nb_lk.rasmol',sep='')
  I1 <- nb > nbthresh
  I2 <- (isH & lk >= lkHthresh) | (!isH & lk >= lkAthresh)
  print(paste(sum(I1),sum(I2),sum(I1&I2)))
  rasmol.whole.chain.vdw.comp(prot,file,pdbdir,res,atom,I1,I2)

##############3

  file <- paste(filepath,prot,'_lke.rasmol',sep='')
  I <- (isH & lke >= lkeHthresh) | (!isH & lke >= lkeAthresh)
  print(sum(I))
  rasmol.whole.chain.vdw(prot,file,pdbdir, res[I] , atom[I] )
  
  file <- paste(filepath,prot,'_lk2e.rasmol',sep='')
  I <- (isH & lk2e >= lk2eHthresh) | (!isH & lk2e >= lk2eAthresh)
  print(sum(I))
  rasmol.whole.chain.vdw(prot,file,pdbdir, res[I] , atom[I] )
  
  file <- paste(filepath,prot,'_lk4e.rasmol',sep='')
  I <- (isH & lk4e >= lk4eHthresh) | (!isH & lk4e >= lk4eAthresh)
  print(sum(I))
  rasmol.whole.chain.vdw(prot,file,pdbdir, res[I] , atom[I] )

  
  file <- paste(filepath,prot,'_lk_lke.rasmol',sep='')
  I1 <- (isH & lk >= lkHthresh) | (!isH & lk >= lkAthresh)
  I2 <- (isH & lke >= lkeHthresh) | (!isH & lke >= lkeAthresh)
  print(paste(sum(I1),sum(I2),sum(I1&I2)))
  rasmol.whole.chain.vdw.comp(prot,file,pdbdir,res,atom,I1,I2)

  file <- paste(filepath,prot,'_lk2_lk2e.rasmol',sep='')
  I1 <- (isH & lk2 >= lk2Hthresh) | (!isH & lk2 >= lk2Athresh)
  I2 <- (isH & lk2e >= lk2eHthresh) | (!isH & lk2e >= lk2eAthresh)
  print(paste(sum(I1),sum(I2),sum(I1&I2)))
  rasmol.whole.chain.vdw.comp(prot,file,pdbdir,res,atom,I1,I2)

  file <- paste(filepath,prot,'_lk4_lk4e.rasmol',sep='')
  I1 <- (isH & lk4 >= lk4Hthresh) | (!isH & lk4 >= lk4Athresh)
  I2 <- (isH & lk4e >= lk4eHthresh) | (!isH & lk4e >= lk4eAthresh)
  print(paste(sum(I1),sum(I2),sum(I1&I2)))
  rasmol.whole.chain.vdw.comp(prot,file,pdbdir,res,atom,I1,I2)

  
  file <- paste(filepath,prot,'_anb5.rasmol',sep='')
  I <- anb5 >= anb5thresh
  print(sum(I))
  rasmol.whole.chain.vdw(prot,file,pdbdir, res[I] , atom[I] )
  
  file <- paste(filepath,prot,'_anb75.rasmol',sep='')
  I <- anb75 >= anb75thresh
  print(sum(I))
  rasmol.whole.chain.vdw(prot,file,pdbdir, res[I] , atom[I] )
  
  file <- paste(filepath,prot,'_anb10.rasmol',sep='')
  I <- anb10 >= anb10thresh
  print(sum(I))
  rasmol.whole.chain.vdw(prot,file,pdbdir, res[I] , atom[I] )
  
  file <- paste(filepath,prot,'_anb5_anb75.rasmol',sep='')
  I1 <- anb5 >= anb5thresh
  I2 <- anb75 >= anb75thresh
  print(paste(sum(I1),sum(I2),sum(I1&I2)))
  rasmol.whole.chain.vdw.comp(prot,file,pdbdir,res,atom,I1,I2)

  file <- paste(filepath,prot,'_anb5_anb10.rasmol',sep='')
  I1 <- anb5 >= anb5thresh
  I2 <- anb10 >= anb10thresh
  print(paste(sum(I1),sum(I2),sum(I1&I2)))
  rasmol.whole.chain.vdw.comp(prot,file,pdbdir,res,atom,I1,I2)

  file <- paste(filepath,prot,'_anb75_anb10.rasmol',sep='')
  I1 <- anb75 >= anb75thresh
  I2 <- anb10 >= anb10thresh
  print(paste(sum(I1),sum(I2),sum(I1&I2)))
  rasmol.whole.chain.vdw.comp(prot,file,pdbdir,res,atom,I1,I2)

  
  
  file <- paste(filepath,prot,'_sasa14_anb5.rasmol',sep='')
  I1 <- sasa14==0
  I2 <- anb5 >= anb5thresh
  print(paste(sum(I1),sum(I2),sum(I1&I2)))
  rasmol.whole.chain.vdw.comp(prot,file,pdbdir,res,atom,I1,I2)

  file <- paste(filepath,prot,'_sasa14_anb75.rasmol',sep='')
  I1 <- sasa14==0
  I2 <- anb75 >= anb75thresh
  print(paste(sum(I1),sum(I2),sum(I1&I2)))
  rasmol.whole.chain.vdw.comp(prot,file,pdbdir,res,atom,I1,I2)

  file <- paste(filepath,prot,'_sasa14_anb10.rasmol',sep='')
  I1 <- sasa14==0
  I2 <- anb10 >= anb10thresh
  print(paste(sum(I1),sum(I2),sum(I1&I2)))
  rasmol.whole.chain.vdw.comp(prot,file,pdbdir,res,atom,I1,I2)

  file <- paste(filepath,prot,'_sasa14nb_anb5.rasmol',sep='')
  I1 <- sasa14==0 & nb > nbthresh
  I2 <- anb5 >= anb5thresh
  print(paste(sum(I1),sum(I2),sum(I1&I2)))
  rasmol.whole.chain.vdw.comp(prot,file,pdbdir,res,atom,I1,I2)

####################
  
  file <- paste(filepath,prot,'_sasa14nb_sasa14anb5.rasmol',sep='')
  I1 <- sasa14==0 & nb > nbthresh
  I2 <- sasa14==0 & anb5 >= anb5thresh
  print(paste(sum(I1),sum(I2),sum(I1&I2)))
  rasmol.whole.chain.vdw.comp(prot,file,pdbdir,res,atom,I1,I2)
  
  file <- paste(filepath,prot,'_sasa14nb_sasa14anb75.rasmol',sep='')
  I1 <- sasa14==0 & nb > nbthresh
  I2 <- sasa14==0 & anb75 >= anb75thresh
  print(paste(sum(I1),sum(I2),sum(I1&I2)))
  rasmol.whole.chain.vdw.comp(prot,file,pdbdir,res,atom,I1,I2)

  file <- paste(filepath,prot,'_sasa14nb_sasa14anb10.rasmol',sep='')
  I1 <- sasa14==0 & nb > nbthresh
  I2 <- sasa14==0 & anb10 >= anb10thresh
  print(paste(sum(I1),sum(I2),sum(I1&I2)))
  rasmol.whole.chain.vdw.comp(prot,file,pdbdir,res,atom,I1,I2)

  
  file <- paste(filepath,prot,'_sasa14anb10_sasa14lk4.rasmol',sep='')
  I1 <- sasa14==0 & anb10 >= anb10thresh
  I2 <- sasa14==0 & (isH & lk4 >= lk4Hthresh) | (!isH & lk4 >= lk4Athresh)
  print(paste(sum(I1),sum(I2),sum(I1&I2)))
  rasmol.whole.chain.vdw.comp(prot,file,pdbdir,res,atom,I1,I2)

  
}
