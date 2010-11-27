#################################################################3
#        with SS and more decoys
#################################################################


##### start analysis

# decided not to use these because I don't quite get what ROCR is doing
#library(gtools,lib.loc='~/software/Rlib')
#library(gdata, lib.loc='~/software/Rlib')
#library(gplots,lib.loc='~/software/Rlib')
#library(ROCR,  lib.loc='~/software/Rlib')


source('~/scripts/buried_polar/buried_polar_util.R')



#############################################3
#############################################
##        change to gsgr vs. gsbr decoys!
#############################################
##############################################3

if(0){
  
setwd('/users/sheffler/project/features/buried_polar/')
#load('rdata/gsbr.decoyset2.dssp.polar.rdata')
#load('rdata/gsgr.pdb.polar.rdata')

#test.calc.functions(gsgr,gsbr)

bth <- 340
gsgr.br <- c();
gsgr.br[gsgr$anb >bth] <- "BR";
gsgr.br[gsgr$anb<=bth] <- "SF"; 
gsgr.br <- as.factor(gsgr.br)
gsbr.br <- c();
gsbr.br[gsbr$anb >bth] <- "BR";
gsbr.br[gsbr$anb<=bth] <- "SF"; 
gsbr.br <- as.factor(gsbr.br)

factors <- list(c('type'),c('type','ss'),c('type','br'),c('type','ss','br'))
gsbr$br <- gsbr.br
gsgr$br <- gsgr.br
dbin  <- calc.joint.bins(gsbr,factors,list(c('type','aa')))
pbin  <- calc.joint.bins(gsgr,factors,list(c('type','aa')))
gsgruf  <- calc.frac.unsatisfied( unsatisfied=gsgr$hbe==0 , bins=pbin )
gsbruf  <- calc.frac.unsatisfied( unsatisfied=gsbr$hbe==0 , bins=dbin )
#gsbrred <- calc.pred(dbin,gsgruf,satisfied=gsbr$hbe!=0,rm.satisfied=T)
#gsgrred <- calc.pred(pbin,gsgruf,satisfied=gsgr$hbe!=0,rm.satisfied=T)
gsbrred <- calc.pred(dbin,gsgruf,satisfied=gsbr$hbe!=0,rm.satisfied=F)
gsgrred <- calc.pred(pbin,gsgruf,satisfied=gsgr$hbe!=0,rm.satisfied=F)
jointuf <- list()
jointuf$aa     <- gsgruf$aa/gsbruf$aa
jointuf$aabr   <- gsgruf$aa/gsbruf$aabr
jointuf$aass   <- gsgruf$aa/gsbruf$aass
jointuf$aassbr <- gsgruf$aa/gsbruf$aassbr
gsbrred2 <- calc.pred(dbin,jointuf,satisfied=gsbr$hbe!=0,rm.satisfied=F)
gsgrred2 <- calc.pred(pbin,jointuf,satisfied=gsgr$hbe!=0,rm.satisfied=F)

gsgr.pdb.pred          <- tapply(gsgr$hbe!=0,    gsgr$pdb, mean)
gsgr.pdb.pred.aa       <- tapply(gsgrred$aa,     gsgr$pdb, mean)
gsgr.pdb.pred.aa.br    <- tapply(gsgrred$aabr,   gsgr$pdb, mean)
gsgr.pdb.pred.aa.ss    <- tapply(gsgrred$aass,   gsgr$pdb, mean)
gsgr.pdb.pred.aa.ss.br <- tapply(gsgrred$aassbr, gsgr$pdb, mean)
gsbr.pdb.pred          <- tapply(gsbr$hbe!=0,    gsbr$pdb, mean)
gsbr.pdb.pred.aa       <- tapply(gsbrred$aa,     gsbr$pdb, mean)
gsbr.pdb.pred.aa.br    <- tapply(gsbrred$aabr,   gsbr$pdb, mean)
gsbr.pdb.pred.aa.ss    <- tapply(gsbrred$aass,   gsbr$pdb, mean)
gsbr.pdb.pred.aa.ss.br <- tapply(gsbrred$aassbr, gsbr$pdb, mean)

gsgr.pdb.pred2          <- tapply(gsgr$hbe!=0,     gsgr$pdb, mean)
gsgr.pdb.pred2.aa       <- tapply(gsgrred2$aa,     gsgr$pdb, mean)
gsgr.pdb.pred2.aa.br    <- tapply(gsgrred2$aabr,   gsgr$pdb, mean)
gsgr.pdb.pred2.aa.ss    <- tapply(gsgrred2$aass,   gsgr$pdb, mean)
gsgr.pdb.pred2.aa.ss.br <- tapply(gsgrred2$aassbr, gsgr$pdb, mean)
gsbr.pdb.pred2          <- tapply(gsbr$hbe!=0,     gsbr$pdb, mean)
gsbr.pdb.pred2.aa       <- tapply(gsbrred2$aa,     gsbr$pdb, mean)
gsbr.pdb.pred2.aa.br    <- tapply(gsbrred2$aabr,   gsbr$pdb, mean)
gsbr.pdb.pred2.aa.ss    <- tapply(gsbrred2$aass,   gsbr$pdb, mean)
gsbr.pdb.pred2.aa.ss.br <- tapply(gsbrred2$aassbr, gsbr$pdb, mean)

pdb.roc          <- calc.roc(gsgr.pdb.pred,         gsbr.pdb.pred          )
pdb.roc.aa       <- calc.roc(gsgr.pdb.pred.aa,      gsbr.pdb.pred.aa       )
pdb.roc.aa.br    <- calc.roc(gsgr.pdb.pred.aa.br,   gsbr.pdb.pred.aa.br    )
pdb.roc.aa.ss    <- calc.roc(gsgr.pdb.pred.aa.ss,   gsbr.pdb.pred.aa.ss.br )
pdb.roc.aa.ss.br <- calc.roc(gsgr.pdb.pred.aa.ss.br,gsbr.pdb.pred.aa.ss.br )
pdb.roc2          <- calc.roc(gsgr.pdb.pred2,         gsbr.pdb.pred2          )
pdb.roc2.aa       <- calc.roc(gsgr.pdb.pred2.aa,      gsbr.pdb.pred2.aa       )
pdb.roc2.aa.br    <- calc.roc(gsgr.pdb.pred2.aa.br,   gsbr.pdb.pred2.aa.br    )
pdb.roc2.aa.ss    <- calc.roc(gsgr.pdb.pred2.aa.ss,   gsbr.pdb.pred2.aa.ss.br )
pdb.roc2.aa.ss.br <- calc.roc(gsgr.pdb.pred2.aa.ss.br,gsbr.pdb.pred2.aa.ss.br )

par(cex=0.5)
plot (pdb.roc         , col=1 , type='l',
      main="PDB/DS2 Buried Polar Protein Level Discrimination ",
      xlab="Fraction False Positives",ylab="Fraction True Positives")
legend(.50,.25,c("ROC",paste(format(round(c(pdb.roc$score,
                                            pdb.aa.roc$score,
                                            pdb.aa.br.roc$score,
                                            pdb.aa.ss.roc$score,
                                            pdb.aa.ss.br.roc$score),3)),
               c('No Bins','Type Bins','Type & Burial Bins',
                 'Type & SS Bins','Type & SS & Burial Bins'))),
       col=c('white','black','black','red','blue','purple'),
       lwd=3)
lines(pdb.aa.roc      , col=1)
lines(pdb.aa.br.roc   , col='red')
lines(pdb.aa.ss.roc   , col='blue')
lines(pdb.aa.ss.br.roc, col='purple')
#dev.print(file=paste('fig/050921/DS2_v_PDB_protein_discrimination_dssp.ps',sep=''),
#          device=postscript,horizontal=F)

par(cex=0.5)
plot (pdb.roc2         , col=1 , type='l',
      main="PDB/DS2 Buried Polar Protein Level Discrimination ",
      xlab="Fraction False Positives",ylab="Fraction True Positives")
legend(.50,.25,c("ROC",paste(format(round(c(pdb.roc2$score,
                                            pdb.aa.roc2$score,
                                            pdb.aa.br.roc2$score,
                                            pdb.aa.ss.roc2$score,
                                            pdb.aa.ss.br.roc2$score),3)),
               c('No Bins','Type Bins','Type & Burial Bins',
                 'Type & SS Bins','Type & SS & Burial Bins'))),
       col=c('white','black','black','red','blue','purple'),
       lwd=3)
lines(pdb.aa.roc2      , col=1)
lines(pdb.aa.br.roc2   , col='red')
lines(pdb.aa.ss.roc2   , col='blue')
lines(pdb.aa.ss.br.roc2, col='purple')
#dev.print(file=paste('fig/050921/DS2_v_PDB_protein_discrimination_dssp.ps',sep=''),
#          device=postscript,horizontal=F)


#uf.cor       <- cor(gsbr.uf,gsgr.uf)
#uf.br.cor    <- cor(gsbr.uf.br,gsgr.uf.br)
#uf.ss.cor    <- cor(gsbr.uf.ss,gsgr.uf.ss)
#uf.ss.br.cor <- cor(gsbr.uf.ss.br,gsgr.uf.ss.br)

#par(mfrow=c(2,1))
#par(cex=0.5)
#plot (gsgr.roc         , col=1 , type='l',
#      main="PDB/DS2 Buried Polar Group Level Discrimination ",
#      xlab="Fraction False Positives",ylab="Fraction True Positives")
#legend(.45,.25,c("ROC  Cor",paste(format(round(c(gsgr.roc$score,
#                                                 aaroc$score,
#                                                 aabrroc$score,
#                                                 aassroc$score,
#                                                 aassbrbr.roc$score),3)),
#                                format(round(c(0,
#                                               uf.cor,
#                                               uf.br.cor,
#                                               uf.ss.cor,
#                                               uf.ss.br.cor),3)),
#               c('No Bins','Type Bins','Type & Burial Bins',
#                 'Type & SS Bins','Type & SS & Burial Bins'))),
#       col=c('white','black','black','red','blue','purple'),
#       lwd=3)
#lines(gsgr.aa.roc      , col=1)
#lines(gsgr.aa.br.roc   , col='red')
#lines(gsgr.aa.ss.roc   , col='blue')
#lines(gsgr.aa.ss.br.roc, col='purple')





}






