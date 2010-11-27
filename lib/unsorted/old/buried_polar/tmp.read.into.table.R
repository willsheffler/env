pdbpolar <- read.table("pdb_polar.data",header=F)
pdbpolar <- pdbpolar[,-c(1,14)]
names(pdbpolar) <- c("pdb",
                     "aa",
                     "res",
                     "atom",
                     "hbe",
                     "hbetot",
                     "dis",
                     "cosh",
                     "cosa",
                     "atomtype",
                     "sasa14",
                     "sasa10",
                     "lk",
                     "lk2",
                     "lk4",
                     "nb",
                     "sasafrac" )


sasa14max <- quantile(pdbpolar$sasa14,0.999)
sasa10max <- quantile(pdbpolar$sasa10,0.999)
pdbpolar$sasa14[ pdbpolar$sasa14 > sasa14max ] <- sasa14max
pdbpolar$sasa10[ pdbpolar$sasa10 > sasa10max ] <- sasa10max

pdbpolar$lk [ is.na(pdbpolar$lk ) ] <- 0
pdbpolar$lk2[ is.na(pdbpolar$lk2) ] <- 0
pdbpolar$lk4[ is.na(pdbpolar$lk4) ] <- 0
lkmax     <- quantile(pdbpolar$lk ,0.999)
lk2max    <- quantile(pdbpolar$lk2,0.999)
lk4max    <- quantile(pdbpolar$lk4,0.999)
pdbpolar$lk [ pdbpolar$lk  > lkmax ]  <- lkmax
pdbpolar$lk2[ pdbpolar$lk2 > lk2max ] <- lk2max
pdbpolar$lk4[ pdbpolar$lk4 > lk4max ] <- lk4max

save(pdbpolar,file="tmprdata/pdbpolar.rdata")


pdbhbond <- read.table("pdb_hbond.data",header=F)
pdbhbond <- pdbhbond[,-c(1)]
names(pdbhbond) <- c("pdb",
                     "donres",
                     "donaa",
                     "donatom",
                     "actres",
                     "actaa",
                     "actatom",
                     "dis",
                     "cosh",
                     "cosa",
                     "dise",
                     "cosae",
                     "coshe"   )
save(pdbhbond,file="tmprdata/pdbhond.rdata")


gvkburial <- read.table('1gvkB.polar.data',header=F)[,-1]
names(gvkburial) <- c("pdb",
                      "aa",
                      "res",
                      "atom",
                      "hbe",
                      "hbetot",
                      "dis",
                      "cosh",
                      "cosa",
                      "atomtype",
                      "sasa14",
                      "sasa10",
                      "sasa7",
                      "lk",
                      "lk2",
                      "lk4",
                      "nb",
                      "sasafrac",
                      "atomnb5",
                      "atomnb75",
                      "atomnb10",
                      "lke",
                      "lk2e",
                      "lk4e"                     
                      )
save(gvkburial,file='~/project/packing_index/rdata/gvkburial.rdata')


pdbburial <- read.table('pdb_burial.data0',header=F)[,-1]
names(pdbburial) <- c("pdb",
                      "aa",
                      "res",
                      "atom",
                      "hbe",
                      "hbetot",
                      "dis",
                      "cosh",
                      "cosa",
                      "atomtype",
                      "sasa14",
                      "sasa10",
                      "sasa7",
                      "lk",
                      "lk2",
                      "lk4",
                      "nb",
                      "sasafrac",
                      "atomnb5",
                      "atomnb75",
                      "atomnb10",
                      "lke",
                      "lk2e",
                      "lk4e"                     
                      )
save(pdbburial,file='pdbburial.rdata')
