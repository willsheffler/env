source('~/scripts/R-utils.R')

#loadshit <- function() {
  ppoc <- read.table("~/w/project/pocket/pdb.pocket.dat",header=F)
  dpoc <- read.table("~/w/project/pocket/ds2.pocket.dat",header=F)
  n <- c('pdb','natom','psurf','pvol','cid','csurf','cvol')
  names(ppoc) <- n
  names(dpoc) <- n
  
  ptotvol <- tapply(ppoc$pvol,ppoc$pdb,mean)
  dtotvol <- tapply(dpoc$pvol,dpoc$pdb,mean)
#}

par(mfrow=c(2,3))

for ( ii in 1:6) {
  pcavvol <- tapply(ppoc$cvol**ii,ppoc$pdb,sum)
  dcavvol <- tapply(dpoc$cvol**ii,dpoc$pdb,sum)
  pvolfrac <- pcavvol/ptotvol
  dvolfrac <- dcavvol/dtotvol
  pvolfrac <- pvolfrac[!is.infinite(pvolfrac)]
  dvolfrac <- dvolfrac[!is.infinite(dvolfrac)]
  plotden(pvolfrac[pvolfrac<0.1], dvolfrac[pvolfrac<0.1] , xlim=c(0.02),
          title="PDB vs. Decoy Set 2, Fraction Cavity 1.4A Volume",lwd=2)
}

pcavvol <- tapply(ppoc$cvol**1,ppoc$pdb,sum)
dcavvol <- tapply(dpoc$cvol**1,dpoc$pdb,sum)
pvolfrac <- pcavvol/ptotvol
dvolfrac <- dcavvol/dtotvol
plotden(pvolfrac[pvolfrac<0.1], dvolfrac[dvolfrac<0.1] , xlim=c(0,0.02),
        title="PDB vs. Decoy Set 2, Fraction Cavity 1.4A Volume",lwd=2)

pcavvol2 <- tapply(ppoc$cvol**2,ppoc$pdb,sum)
dcavvol2 <- tapply(dpoc$cvol**2,dpoc$pdb,sum)
pvolfrac2 <- log(pcavvol2/ptotvol**2)
dvolfrac2 <- log(dcavvol2/dtotvol**2)
plotden(pvolfrac2[pvolfrac<0.1], dvolfrac2[dvolfrac<0.1] , xlim=c(0,0.0001),
        title="PDB vs. Decoy Set 2, Fraction Cavity 1.4A Volume",lwd=2)

pcavvol3 <- tapply(ppoc$cvol**3,ppoc$pdb,sum)
dcavvol3 <- tapply(dpoc$cvol**3,dpoc$pdb,sum)
pvolfrac3 <- pcavvol3/ptotvol
dvolfrac3 <- dcavvol3/dtotvol
plotden(pvolfrac3[pvolfrac<0.1], dvolfrac3[dvolfrac<0.1] ,
        title="PDB vs. Decoy Set 2, Fraction Cavity 1.4A Volume",lwd=2)

pcavvol4 <- tapply(ppoc$cvol**4,ppoc$pdb,sum)
dcavvol4 <- tapply(dpoc$cvol**4,dpoc$pdb,sum)
pvolfrac4 <- pcavvol4/ptotvol
dvolfrac4 <- dcavvol4/dtotvol
plotden(pvolfrac4[pvolfrac<0.1], dvolfrac4[dvolfrac<0.1] ,
        title="PDB vs. Decoy Set 2, Fraction Cavity 1.4A Volume",lwd=2)




#dev.copy2eps(file="~/figures/gm060925/cavity_volume_pdb_v_decoy")
