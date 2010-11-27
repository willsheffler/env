#################################################################3
#        with SS and more decoys
#################################################################


##### start analysis

# decided not to use these because I don't quite get what ROCR is doing
#library(gtools,lib.loc='~/software/Rlib')
#library(gdata, lib.loc='~/software/Rlib')
#library(gplots,lib.loc='~/software/Rlib')
#library(ROCR,  lib.loc='~/software/Rlib')



calc.joint.bins <- function(bp,fieldlist,namesubs,fsep=".",nsep="") {
  if(missing(namesubs))
    sub <- list()
  bins <- list()
  for(ii in 1:length(fieldlist)) {
    fields <- fieldlist[[ii]]
    name <- paste(fieldlist[[ii]],collapse=nsep)
    for(sub in namesubs)
      name <- sub(sub[1],sub[2],name)
    print(paste(c("calculating joint factor on fields:",fields),collapse=" "))
    jointfield <- as.character(bp[[fields[1]]])
    if(length(fields)>1) {
      for(jj in 2:length(fields)) {
        field <- fields[jj]
        jointfield <- paste(jointfield,as.character(bp[[field]]),sep=fsep)
      }
    }
    bins[[name]] <- as.factor(jointfield)
  }
  return(bins)  
}

calc.frac.unsatisfied <- function(unsatisfied,bins) {
  uf <- list()
  for(bin in names(bins)) {
    ct <- tapply(unsatisfied, bins[[bin]], length)
    uc <- tapply(unsatisfied, bins[[bin]], sum)
    uf[[bin]] <- uc/ct
  }
  return(uf)
}

calc.pred <- function(bins,predictor,satisfied=F,rm.satisfied=F) {
  pred <- list()
  for(bin in names(bins)) {
    tmp <- bins[[bin]]
    levels(tmp) <- predictor[[bin]]
    tmp <- as.numeric(as.character(tmp))
    if(rm.satisfied)
      tmp <- tmp[!satisfied]
    else if(any(satisfied))
      tmp[satisfied] <- 1
    pred[[bin]] <- tmp
  }
  return(pred)
}


# function to plot roc curves
# labels must be true and false
calc.roc <- function(posval,negval,n=10,highgood=T) {
  pred  <- c(posval,negval)
  label <- c(rep(T,length(posval)),rep(F,length(negval)))
  tpr <- 0:(n-1)/(n-1)
  fpr <- 0:(n-1)/(n-1)
  q <- quantile(pred,0:(n-1)/(n-1))  
  for(ii in 2:(n-1)) {
    x <- sum(pred<=q[ii] & !label)/sum( !label )
    y <- sum(pred<=q[ii] &  label)/sum(  label )
    fpr[ii] <- x
    tpr[ii] <- y
  }
  if(highgood) {
    tpr <- rev(1-tpr)
    fpr <- rev(1-fpr)
  }
  score <- 0
  for(ii in 2:n) {
    #print(paste(ii,score, mean(c(tpr[ii-1],tpr[ii])),abs(fpr[ii]-fpr[ii-1])))
    score <- score + mean(c(tpr[ii-1],tpr[ii]))*abs(fpr[ii]-fpr[ii-1])
  }
  return(list(y=tpr,x=fpr,score=score))
}


if(0){
  
setwd('/users/sheffler/project/features/buried_polar/')
load('rdata/dp.decoyset2.dssp.polar.rdata')
load('rdata/pp.pdb.polar.rdata')

#test.calc.functions(pp,dp)

bth <- 340
pp.br <- c();
pp.br[pp$anb >bth] <- "BR";
pp.br[pp$anb<=bth] <- "SF"; 
pp.br <- as.factor(pp.br)
dp.br <- c();
dp.br[dp$anb >bth] <- "BR";
dp.br[dp$anb<=bth] <- "SF"; 
dp.br <- as.factor(dp.br)

factors <- list(c('type'),c('type','ss'),c('type','br'),c('type','ss','br'))
dp$br <- dp.br
pp$br <- pp.br
dbin  <- calc.joint.bins(dp,factors,list(c('type','aa')))
pbin  <- calc.joint.bins(pp,factors,list(c('type','aa')))
ppuf  <- calc.frac.unsatisfied( unsatisfied=pp$hbe==0 , bins=pbin )
dpuf  <- calc.frac.unsatisfied( unsatisfied=dp$hbe==0 , bins=dbin )
#dpred <- calc.pred(dbin,ppuf,satisfied=dp$hbe!=0,rm.satisfied=T)
#ppred <- calc.pred(pbin,ppuf,satisfied=pp$hbe!=0,rm.satisfied=T)
dpred <- calc.pred(dbin,ppuf,satisfied=dp$hbe!=0,rm.satisfied=F)
ppred <- calc.pred(pbin,ppuf,satisfied=pp$hbe!=0,rm.satisfied=F)
jointuf <- list()
jointuf$aa     <- ppuf$aa/dpuf$aa
jointuf$aabr   <- ppuf$aa/dpuf$aabr
jointuf$aass   <- ppuf$aa/dpuf$aass
jointuf$aassbr <- ppuf$aa/dpuf$aassbr
dpred2 <- calc.pred(dbin,jointuf,satisfied=dp$hbe!=0,rm.satisfied=F)
ppred2 <- calc.pred(pbin,jointuf,satisfied=pp$hbe!=0,rm.satisfied=F)

pp.pdb.pred          <- tapply(pp$hbe!=0,    pp$pdb, mean)
pp.pdb.pred.aa       <- tapply(ppred$aa,     pp$pdb, mean)
pp.pdb.pred.aa.br    <- tapply(ppred$aabr,   pp$pdb, mean)
pp.pdb.pred.aa.ss    <- tapply(ppred$aass,   pp$pdb, mean)
pp.pdb.pred.aa.ss.br <- tapply(ppred$aassbr, pp$pdb, mean)
dp.pdb.pred          <- tapply(dp$hbe!=0,    dp$pdb, mean)
dp.pdb.pred.aa       <- tapply(dpred$aa,     dp$pdb, mean)
dp.pdb.pred.aa.br    <- tapply(dpred$aabr,   dp$pdb, mean)
dp.pdb.pred.aa.ss    <- tapply(dpred$aass,   dp$pdb, mean)
dp.pdb.pred.aa.ss.br <- tapply(dpred$aassbr, dp$pdb, mean)

pp.pdb.pred2          <- tapply(pp$hbe!=0,     pp$pdb, mean)
pp.pdb.pred2.aa       <- tapply(ppred2$aa,     pp$pdb, mean)
pp.pdb.pred2.aa.br    <- tapply(ppred2$aabr,   pp$pdb, mean)
pp.pdb.pred2.aa.ss    <- tapply(ppred2$aass,   pp$pdb, mean)
pp.pdb.pred2.aa.ss.br <- tapply(ppred2$aassbr, pp$pdb, mean)
dp.pdb.pred2          <- tapply(dp$hbe!=0,     dp$pdb, mean)
dp.pdb.pred2.aa       <- tapply(dpred2$aa,     dp$pdb, mean)
dp.pdb.pred2.aa.br    <- tapply(dpred2$aabr,   dp$pdb, mean)
dp.pdb.pred2.aa.ss    <- tapply(dpred2$aass,   dp$pdb, mean)
dp.pdb.pred2.aa.ss.br <- tapply(dpred2$aassbr, dp$pdb, mean)

pdb.roc          <- calc.roc(pp.pdb.pred,         dp.pdb.pred          )
pdb.roc.aa       <- calc.roc(pp.pdb.pred.aa,      dp.pdb.pred.aa       )
pdb.roc.aa.br    <- calc.roc(pp.pdb.pred.aa.br,   dp.pdb.pred.aa.br    )
pdb.roc.aa.ss    <- calc.roc(pp.pdb.pred.aa.ss,   dp.pdb.pred.aa.ss.br )
pdb.roc.aa.ss.br <- calc.roc(pp.pdb.pred.aa.ss.br,dp.pdb.pred.aa.ss.br )
pdb.roc2          <- calc.roc(pp.pdb.pred2,         dp.pdb.pred2          )
pdb.roc2.aa       <- calc.roc(pp.pdb.pred2.aa,      dp.pdb.pred2.aa       )
pdb.roc2.aa.br    <- calc.roc(pp.pdb.pred2.aa.br,   dp.pdb.pred2.aa.br    )
pdb.roc2.aa.ss    <- calc.roc(pp.pdb.pred2.aa.ss,   dp.pdb.pred2.aa.ss.br )
pdb.roc2.aa.ss.br <- calc.roc(pp.pdb.pred2.aa.ss.br,dp.pdb.pred2.aa.ss.br )

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


#uf.cor       <- cor(dp.uf,pp.uf)
#uf.br.cor    <- cor(dp.uf.br,pp.uf.br)
#uf.ss.cor    <- cor(dp.uf.ss,pp.uf.ss)
#uf.ss.br.cor <- cor(dp.uf.ss.br,pp.uf.ss.br)

#par(mfrow=c(2,1))
#par(cex=0.5)
#plot (pp.roc         , col=1 , type='l',
#      main="PDB/DS2 Buried Polar Group Level Discrimination ",
#      xlab="Fraction False Positives",ylab="Fraction True Positives")
#legend(.45,.25,c("ROC  Cor",paste(format(round(c(pp.roc$score,
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
#lines(pp.aa.roc      , col=1)
#lines(pp.aa.br.roc   , col='red')
#lines(pp.aa.ss.roc   , col='blue')
#lines(pp.aa.ss.br.roc, col='purple')





}













############### tests, utils, plots ##############

# function to plot roc curves
# labels must be true and false
calcroc <- function(pred,label,n=10,highgood=T) {
  #pred  <- c(posval,negval)
  #label <- c(rep(T,length(posval)),rep(F,length(negval)))
  tpr <- 0:(n-1)/(n-1)
  fpr <- 0:(n-1)/(n-1)
  q <- quantile(pred,0:(n-1)/(n-1))  
  for(ii in 2:(n-1)) {
    x <- sum(pred<=q[ii] & !label)/sum( !label )
    y <- sum(pred<=q[ii] &  label)/sum(  label )
    fpr[ii] <- x
    tpr[ii] <- y
  }
  if(highgood) {
    tpr <- rev(1-tpr)
    fpr <- rev(1-fpr)
  }
  score <- 0
  for(ii in 2:n) {
    #print(paste(ii,score, mean(c(tpr[ii-1],tpr[ii])),abs(fpr[ii]-fpr[ii-1])))
    score <- score + mean(c(tpr[ii-1],tpr[ii]))*abs(fpr[ii]-fpr[ii-1])
  }
  return(list(y=tpr,x=fpr,score=score))
}


read.raw.data <- function(pp,dp) {

#dp1 <- read.table('ds2_polar_dssp_0.data',header=F)[,-c(1,10:13,16)]
#dp2 <- read.table('ds2_polar_dssp_1.data',header=F)[,-c(1,10:13,16)]
#dp3 <- read.table('ds2_polar_dssp_2.data',header=F)[,-c(1,10:13,16)]
#dp4 <- read.table('ds2_polar_dssp_3.data',header=F)[,-c(1,10:13,16)]
#dp <- rbind(dp1,dp2,dp3,dp4)
#names(dp) <- c('prot','pdb','aa','res','ss','atom','type',
#               'hbe','nb','anb')
#save(dp,file='~/project/features/buried_polar/rdata/dp.decoyset2.dssp.polar.rdata')

  atomtype <- dp$type
  type <- c()
  type[dp$type== 8] <- "NH"
  type[dp$type==13] <- "OH"
  type[dp$type==14] <- "OA"
  type[dp$type==15] <- "OC"
  type[dp$type==20] <- "OB"
  type[dp$type==22] <- "HP"
  type[dp$type==25] <- "HB"
  isbb <- rep(F,length(type))
  isbb[dp$type==20] <- T
  isbb[dp$type==25] <- T
  type <- paste(as.character(dp$aa),as.character(type),sep='.')
  dp$type <- as.factor(type)
  dp <- cbind(dp,isbb)
  dp <- cbind(dp,atomtype)
  rm(atomtype,isbb)
  pdb <- levels(dp$pdb)
  pdb <- strsplit(gsub('//','/',as.character(pdb)),'/',fixed=T)
  pdb <- lapply(pdb,function(x) paste(x[length(x)-2],x[length(x)-3],sep='.'))
  pdb <- as.character(pdb)
  dp$prot <- dp$pdb
  levels(dp$prot) <- pdb
  rm(pdb)
  save(dp,file='rdata/dp.decoyset2.dssp.polar.rdata')
  
  atomtype <- pp$type
  type <- c()
  type[pp$type== 8] <- "NH"
  type[pp$type==13] <- "OH"
  type[pp$type==14] <- "OA"
  type[pp$type==15] <- "OC"
  type[pp$type==20] <- "OB"
  type[pp$type==22] <- "HP"
  type[pp$type==25] <- "HB"
  isbb <- rep(F,length(type))
  isbb[pp$type==20] <- T
  isbb[pp$type==25] <- T
  type <- paste(as.character(pp$aa),as.character(type),sep='.')
  pp$type <- as.factor(type)
  pp <- cbind(pp,isbb)
  pp <- cbind(pp,atomtype)
  pp$prot <- pp$pdb
  save(pp,file='rdata/pp.pdb.polar.rdata')
}



test.calc.functions <- function(pdbpolar,decoypolar) {
  pp.bt <- tapply(pdbpolar$anb, list(pdbpolar$type,pdbpolar$ss), median)
  dp.bt <- tapply(decoypolar$anb, list(decoypolar$type,decoypolar$ss), median)
  bth <- 340
  pp.br <- c();
  pp.br[pdbpolar$anb >bth] <- "BR";
  pp.br[pdbpolar$anb<=bth] <- "SF"; 
  pp.br <- as.factor(pp.br)
  dp.br <- c();
  dp.br[decoypolar$anb >bth] <- "BR";
  dp.br[decoypolar$anb<=bth] <- "SF"; 
  dp.br <- as.factor(dp.br)
  
  
                                        # computing UNSFRAC statistics
  dp.aa <- as.factor(decoypolar$type)
  pp.aa <- as.factor(pdbpolar$type)
  pp.ct <- tapply(pdbpolar$hbe,   pp.aa, length)
  pp.uc <- tapply(pdbpolar$hbe==0,pp.aa, sum)
  pp.uf <- pp.uc/pp.ct
  dp.ct <- tapply(decoypolar$hbe,   dp.aa, length)
  dp.uc <- tapply(decoypolar$hbe==0,dp.aa, sum)
  dp.uf <- dp.uc/dp.ct
  
  dp.aa.ss <- as.factor(paste(as.character(decoypolar$type),as.character(decoypolar$ss),sep='.'))
  pp.aa.ss <- as.factor(paste(as.character(pdbpolar$type),as.character(pdbpolar$ss),sep='.'))
  pp.ct.ss <- tapply(pdbpolar$hbe,    pp.aa.ss, length)
  pp.uc.ss <- tapply(pdbpolar$hbe==0, pp.aa.ss, sum)
  pp.uf.ss <- pp.uc.ss / pp.ct.ss
  dp.ct.ss <- tapply(decoypolar$hbe,    dp.aa.ss, length)
  dp.uc.ss <- tapply(decoypolar$hbe==0, dp.aa.ss, sum)
  dp.uf.ss <- dp.uc.ss / dp.ct.ss
  
  dp.aa.br <- as.factor(paste(as.character(decoypolar$type),as.character(dp.br),sep='.'))
  pp.aa.br <- as.factor(paste(as.character(pdbpolar$type),as.character(pp.br),sep='.'))
  pp.ct.br <- tapply(pdbpolar$hbe,    pp.aa.br, length)
  pp.uc.br <- tapply(pdbpolar$hbe==0, pp.aa.br, sum)
  pp.uf.br <- pp.uc.br / pp.ct.br
  dp.ct.br <- tapply(decoypolar$hbe,    dp.aa.br, length)
  dp.uc.br <- tapply(decoypolar$hbe==0, dp.aa.br, sum)
  dp.uf.br <- dp.uc.br / dp.ct.br

  dp.aa.ss.br <- as.factor(paste(as.character(decoypolar$type),
                                 as.character(decoypolar$ss),
                                 as.character(dp.br),sep='.'))
  pp.aa.ss.br <- as.factor(paste(as.character(pdbpolar$type),
                                 as.character(pdbpolar$ss),
                                 as.character(pp.br),sep='.'))
  pp.ct.ss.br <- tapply(pdbpolar$hbe,    pp.aa.ss.br, length)
  pp.uc.ss.br <- tapply(pdbpolar$hbe==0, pp.aa.ss.br, sum)
  pp.uf.ss.br <- pp.uc.ss.br / pp.ct.ss.br
  dp.ct.ss.br <- tapply(decoypolar$hbe,    dp.aa.ss.br, length)
  dp.uc.ss.br <- tapply(decoypolar$hbe==0, dp.aa.ss.br, sum)
  dp.uf.ss.br <- dp.uc.ss.br / dp.ct.ss.br

  # computing predictions based on unsfrac bins
  # for all the buried polar groups (only buried!)
  Idp <- dp$hbe==0
  Ipp <- pdbpolar$hbe==0
  pp.pred <- rep(sum(Ipp)/length(Ipp),sum(Ipp))
  dp.pred <- rep(sum(Ipp)/length(Ipp),sum(Idp))

  pptmp <- pp.aa
  dptmp <- dp.aa
  levels(pptmp) <- pp.uf
  levels(dptmp) <- pp.uf
  pp.aa.pred    <- as.numeric(as.character(pptmp))[Ipp]
  dp.aa.pred    <- as.numeric(as.character(dptmp))[Idp]

  pptmp <- pp.aa.ss
  dptmp <- dp.aa.ss
  levels(pptmp) <- pp.uf.ss
  levels(dptmp) <- pp.uf.ss
  pp.aa.ss.pred    <- as.numeric(as.character(pptmp))[Ipp]
  dp.aa.ss.pred    <- as.numeric(as.character(dptmp))[Idp]

  pptmp <- pp.aa.br
  dptmp <- dp.aa.br
  levels(pptmp) <- pp.uf.br
  levels(dptmp) <- pp.uf.br
  pp.aa.br.pred    <- as.numeric(as.character(pptmp))[Ipp]
  dp.aa.br.pred    <- as.numeric(as.character(dptmp))[Idp]

  pptmp <- pp.aa.ss.br
  dptmp <- dp.aa.ss.br
  levels(pptmp) <- pp.uf.ss.br
  levels(dptmp) <- pp.uf.ss.br
  pp.aa.ss.br.pred    <- as.numeric(as.character(pptmp))[Ipp]
  dp.aa.ss.br.pred    <- as.numeric(as.character(dptmp))[Idp]

  factors <- list(c('type'),c('type','ss'),c('type','br'),c('type','ss','br'))
  dp$br <- dp.br
  pp$br <- pp.br
  dbin  <- calc.joint.bins(dp,factors,list(c('type','aa')))
  pbin  <- calc.joint.bins(pp,factors,list(c('type','aa')))
  ppuf  <- calc.frac.unsatisfied( unsatisfied=pp$hbe==0 , bins=pbin )
  dpuf  <- calc.frac.unsatisfied( unsatisfied=dp$hbe==0 , bins=dbin )
  dpred <- calc.pred(dbin,ppuf,satisfied=dp$hbe!=0,rm.satisfied=T)
  ppred <- calc.pred(pbin,ppuf,satisfied=pp$hbe!=0,rm.satisfied=T)

  if((any(dbin$aa     != dp.aa)       |
      any(dbin$aabr   != dp.aa.br)    |
      any(dbin$aass   != dp.aa.ss)    |
      any(dbin$aassbr != dp.aa.ss.br) |
      any(pbin$aa     != pp.aa)       |
      any(pbin$aabr   != pp.aa.br)    |
      any(pbin$aass   != pp.aa.ss)    |
      any(pbin$aassbr != pp.aa.ss.br)   ))
    {
      print("error in calc.joint.bins!!!!!")
    } else {
      print("bins all agree")
    }


  if((any(pp.uf       != ppuf$aa    ) |
      any(pp.uf.ss    != ppuf$aass  ) |
      any(pp.uf.br    != ppuf$aabr  ) |
      any(pp.uf.ss.br != ppuf$aassbr) |
      any(dp.uf       != dpuf$aa    ) |
      any(dp.uf.ss    != dpuf$aass  ) |
      any(dp.uf.br    != dpuf$aabr  ) |
      any(dp.uf.ss.br != dpuf$aassbr)   ))
     {
       print("error in calc.unsfrac!!!!")
     } else {
      print("unsfracs all agree")
    }

    
  if((any(dpred$aa     != dp.aa.pred)       |
      any(dpred$aabr   != dp.aa.br.pred)    |
      any(dpred$aass   != dp.aa.ss.pred)    |
      any(dpred$aassbr != dp.aa.ss.br.pred) |
      any(ppred$aa     != pp.aa.pred)       |
      any(ppred$aabr   != pp.aa.br.pred)    |
      any(ppred$aass   != pp.aa.ss.pred)    |
      any(ppred$aassbr != pp.aa.ss.br.pred)   ))
    {
      print("error in calc.pred!!!!!")
    } else {
      print("predictions all agree")
    }

}


  
make.pp.dp.plots <- function(pp,dp) {
  setwd('/users/sheffler/project/features/buried_polar/')

  aanames3 <- c("ALA","ARG","ASN","ASP","CYS","GLN","GLU","GLY","HIS","ILE",
                "LEU","LYS","MET","PHE","PRO","SER","THR","TRP","TYR","VAL")

  otypes.nosp <- unique(pp$type[!pp$isbb])[
                                           order(substr(unique(as.character(pp$type[!pp$isbb])),5,6))]

  sc <- rev(c("ARG.HP","ASN.HP","GLN.HP","THR.HP","SER.HP","LYS.HP",
              "TYR.HP","TRP.HP","HIS.HP","CYS.HP",NA,"HIS.NH",NA,"ASN.OA",
              "GLN.OA",NA,"ASP.OC","GLU.OC",NA,"SER.OH","THR.OH","TYR.OH"))
  bbh <- rev(c(paste(aanames3[-15],'HB',sep='.')))
  bbo <- rev(c(paste(aanames3,'OB',sep='.')))
  bb  <- c(bbh,NA,bbo)
  otypes    <- c(bb,NA,sc)


############ graphics params ############

  dev.off()
  x11(height=10,width=7.5)
  par(mar=c(5,6,3,1))
  par(las=2)
  par(mfrow=c(1,2))
  par(cex=0.6)

########## burial thresholds #############

  pp.bt <- tapply(pp$anb, list(pp$type,pp$ss), median)
  dp.bt <- tapply(dp$anb, list(dp$type,dp$ss), median)
  bth <- 340
  pp.br <- c();
  pp.br[pp$anb >bth] <- "BR";
  pp.br[pp$anb<=bth] <- "SF"; 
  pp.br <- as.factor(pp.br)
  dp.br <- c();
  dp.br[dp$anb >bth] <- "BR";
  dp.br[dp$anb<=bth] <- "SF"; 
  dp.br <- as.factor(dp.br)


########## uns frac stuff ################

  pp.pc <- tapply(pp$prot,  pp$type, function(x) length(unique(x)))
  pp.ct <- tapply(pp$hbe,   pp$type, length)
  pp.uc <- tapply(pp$hbe==0,pp$type, sum)
  pp.uf <- pp.uc/pp.ct
  dp.pc <- tapply(dp$prot,  dp$type, function(x) length(unique(x)))
  dp.ct <- tapply(dp$hbe,   dp$type, length)
  dp.uc <- tapply(dp$hbe==0,dp$type, sum)
  dp.uf <- dp.uc/dp.ct
  ur    <- dp.uf/pp.uf

  pp.pc.ss <- tapply(pp$prot,   list(pp$type,pp$ss), function(x) length(unique(x)))
  pp.ct.ss <- tapply(pp$hbe,    list(pp$type,pp$ss), length)
  pp.uc.ss <- tapply(pp$hbe==0, list(pp$type,pp$ss), sum)
  pp.uf.ss <- pp.uc.ss / pp.ct.ss
  dp.pc.ss <- tapply(dp$prot,   list(dp$type,dp$ss), function(x) length(unique(x)))
  dp.ct.ss <- tapply(dp$hbe,    list(dp$type,dp$ss), length)
  dp.uc.ss <- tapply(dp$hbe==0, list(dp$type,dp$ss), sum)
  dp.uf.ss <- dp.uc.ss / dp.ct.ss
  ur.ss    <- dp.uf.ss / pp.uf.ss

  pp.pc.br <- tapply(pp$prot,   list(pp$type,pp.br), function(x) length(unique(x)))
  pp.ct.br <- tapply(pp$hbe,    list(pp$type,pp.br), length)
  pp.uc.br <- tapply(pp$hbe==0, list(pp$type,pp.br), sum)
  pp.uf.br <- pp.uc.br / pp.ct.br
  dp.pc.br <- tapply(dp$prot,   list(dp$type,dp.br), function(x) length(unique(x)))
  dp.ct.br <- tapply(dp$hbe,    list(dp$type,dp.br), length)
  dp.uc.br <- tapply(dp$hbe==0, list(dp$type,dp.br), sum)
  dp.uf.br <- dp.uc.br / dp.ct.br
  ur.br    <- dp.uf.br / pp.uf.br

  pp.ssbr    <- as.factor(paste(as.character(pp$ss),as.character(pp.br),sep='.'))
  pp.pc.ssbr <- tapply(pp$prot,   list(pp$type,pp.ssbr), function(x) length(unique(x)))
  pp.ct.ssbr <- tapply(pp$hbe,    list(pp$type,pp.ssbr), length)
  pp.uc.ssbr <- tapply(pp$hbe==0, list(pp$type,pp.ssbr), sum)
  pp.uf.ssbr <- pp.uc.ssbr / pp.ct.ssbr
  dp.ssbr    <- as.factor(paste(as.character(dp$ss),as.character(dp.br),sep='.'))
  dp.pc.ssbr <- tapply(dp$prot,   list(dp$type,dp.ssbr), function(x) length(unique(x)))
  dp.ct.ssbr <- tapply(dp$hbe,    list(dp$type,dp.ssbr), length)
  dp.uc.ssbr <- tapply(dp$hbe==0, list(dp$type,dp.ssbr), sum)
  dp.uf.ssbr <- dp.uc.ssbr / dp.ct.ssbr
  ur.ssbr    <- dp.uf.ssbr / pp.uf.ssbr

############ median burial thresholds

  barplot(rbind(pp.bt[,'E'][otypes],pp.bt[,'H'][otypes],pp.bt[,'L'][otypes]),
          beside=T,col=c('yellow','red','blue'),space=c(0.2,1),horiz=T,xlim=c(200,400),xpd=F,
          main="PDB Median ANB Burial by SS\n(blue=loop,red=helix,yellow=sheet)")
  barplot(rbind(dp.bt[,'E'][otypes],dp.bt[,'H'][otypes],dp.bt[,'L'][otypes]),
          beside=T,col=c('yellow','red','blue'),space=c(0.2,1),horiz=T,xlim=c(200,400),xpd=F,
          main="DS2 Median ANB Burial by SS\n(blue=loop,red=helix,yellow=sheet)")
                                        #dev.print(file='fig/050921/anb_median_burial_pdb_dssp.ps',device=postscript,horizontal=F)

############ just aa bins #########

                                        #dev.print(file='fig/050921/type_counts_pdb_dssp.ps',device=postscript,horizontal=F)
  barplot(rbind(dp.uc[otypes],pp.uc[otypes]),
          beside=T,col=c('grey','black'),space=c(0.2,1),horiz=T,xlim=c(0,5000),xpd=F,
          main="AA bin DS2 Buried Polar Counts")
  barplot(rbind(dp.pc[otypes],pp.pc[otypes]),
          beside=T,col=c(1,'grey'),space=c(0.2,1),horiz=T,xlim=c(0,max(dp.pc)),xpd=F,
          main="AA bin DS2 Buried Polar Protein Counts")
                                        #dev.print(file='fig/050921/count_aa_bins_dssp.ps',device=postscript,horizontal=F)

  barplot(rbind(pp.uf[otypes],dp.uf[otypes]),
          beside=T,col=c('grey','black'),space=c(0.2,1),horiz=T,xlim=c(0,1),xpd=F,
          main="Buried Polar UNS Frac PDB=Black, DS2=Red")
  barplot(ur[otypes],
          beside=T,col=c(1),space=c(0.2,1),horiz=T,xpd=F,
          main="Buried Polar UNS Frac Ratio DS2/PDB")
  lines(c(1,1),c(0,10000),col='purple')
                                        #dev.print(file='fig/050921/unsfrac_ratio_ds2_pdb_aa_dssp.ps',device=postscript,horizontal=F)



########## burial bins

  barplot(rbind(pp.uc.br[,'BR'][otypes],pp.uc.br[,'SF'][otypes]),
          beside=T,col=c('black','grey'),space=c(0.2,1),horiz=T,xlim=c(0,1000),xpd=F,
          main="PDB Buried Polar UNS Counts by Burialn\n(buried=black,surface=grey)")
  barplot(rbind(dp.uc.br[,'BR'][otypes],dp.uc.br[,'SF'][otypes]),
          beside=T,col=c('black','grey'),space=c(0.2,1),horiz=T,xlim=c(0,1000),xpd=F,
          main="DS2 Buried Polar UNS Counts by Burial\n(buried=black,surface=grey)")
                                        #dev.print(file='fig/050921/type_counts_pdb_ds2_aa_br_dssp.ps',device=postscript,horizontal=F)

  barplot(rbind(pp.uf.br[,'BR'][otypes],pp.uf.br[,'SF'][otypes]),
          beside=T,col=c('black','grey'),space=c(0.2,1),horiz=T,xpd=F,
          main="PDB Buried Polar UNS Frac by Burial\n(buried=black,surface=grey)")
  barplot(rbind(dp.uf.br[,'BR'][otypes],dp.uf.br[,'SF'][otypes]),
          beside=T,col=c('black','grey'),space=c(0.2,1),horiz=T,xpd=F,
          main="DS2 Buried Polar UNS Frac by Burial\n(buried=black,surface=grey)")
                                        #dev.print(file='fig/050921/unsfrac_pdb_ds2_aa_br_dssp.ps',device=postscript,horizontal=F)

 
  barplot(rbind(ur.br[,'BR'][otypes],ur.br[,'SF'][otypes]),
          beside=T,col=c('black','grey'),space=c(0.2,1),horiz=T,xpd=F,
          main="DS2 Buried Polar UNS Frac Ratio DS2/PDB by Burial
         buried=black,surface=grey)")
  lines(c(1,1),c(0,10000),col='purple')
  pc.br <- pmin(dp.pc.br,pp.pc.br)
  barplot(rbind(pc.br[,'BR'][otypes],pc.br[,'SF'][otypes]),
          beside=T,col=c('black','grey'),space=c(0.2,1),horiz=T,xpd=F,
          main="DS2 Buried Polar UNS Counts by Burial\n(buried=black,surface=grey)")
                                        #dev.print(file='fig/050921/unsfrac_ratio_ds2_pdb_aa_br_dssp.ps',device=postscript,horizontal=F)



######### ss bins ###########
                                        #for(types in c('sc','bbo','bbh')) {
  barplot(rbind(pp.uc.ss[,'E'][otypes],pp.uc.ss[,'H'][otypes],pp.uc.ss[,'L'][otypes]),
          beside=T,col=c('yellow','red','blue'),space=c(0.2,1),horiz=T,xlim=c(0,1000),xpd=F,
          main="PDB Buried Polar UNS Counts by SS
                (blue=loop,red=helix,yellow=sheet)")
  barplot(rbind(dp.uc.ss[,'E'][otypes],dp.uc.ss[,'H'][otypes],dp.uc.ss[,'L'][otypes]),
          beside=T,col=c('yellow','red','blue'),space=c(0.2,1),horiz=T,xlim=c(0,1000),xpd=F,
          main="DS2 Buried Polar UNS Counts by SS
                (blue=loop,red=helix,yellow=sheet)")
  dev.print(file='fig/050921/type_counts_pdb_ds2_aa_ss_dssp.ps',device=postscript,horizontal=F)

  # for R.R.
  otypes <- c("ALA.HB","LEU.HB","PHE.HB","TYR.HB","ARG.HB","ASP.HB",NA,"ALA.OB","LEU.OB","PHE.OB","PRO.OB","TYR.OB","ARG.OB","ASP.OB",NA,"TYR.HP","TYR.OH","CYS.HP","TRP.HP","ASN.HP","ASN.OA","ARG.HP")
 

  barplot(rbind(pp.uf.ss[,'E'][otypes],pp.uf.ss[,'H'][otypes],pp.uf.ss[,'L'][otypes]),
          beside=T,col=c('yellow','red','blue'),space=c(0.2,1),horiz=T,xpd=F,
          main="PDB Buried Polar UNS Frac by SS
                (blue=loop,red=helix,yellow=sheet)")
  barplot(rbind(dp.uf.ss[,'E'][otypes],dp.uf.ss[,'H'][otypes],dp.uf.ss[,'L'][otypes]),
          beside=T,col=c('yellow','red','blue'),space=c(0.2,1),horiz=T,xpd=F,
          main="DS2 Buried Polar UNS Frac by SS
                (blue=loop,red=helix,yellow=sheet)")
  dev.print(file='fig/050921/unsfrac_pdb_ds2_aa_ss_dssp.ps',device=postscript,horizontal=F)

  barplot(rbind(ur.ss[,'E'][otypes],ur.ss[,'H'][otypes],ur.ss[,'L'][otypes]),
          beside=T,col=c('yellow','red','blue'),space=c(0.2,1),horiz=T,xpd=F,
          main="DS2 Buried Polar UNS Frac Ratio DS2/PDB by SS
                (blue=loop,red=helix,yellow=sheet)")
  lines(c(1,1),c(0,10000),col='purple')
  pc.ss <- pmin(pp.pc.ss,dp.pc.ss)
  barplot(rbind(pc.ss[,'E'][otypes],pc.ss[,'H'][otypes],pc.ss[,'L'][otypes]),
          beside=T,col=c('yellow','red','blue'),space=c(0.2,1),horiz=T,xpd=F,
          main="PDB Buried Polar UNS Counts by SS
                (blue=loop,red=helix,yellow=sheet)")
  dev.print(file='fig/050921/unsfrac_ratio_ds2_pdb_aa_ss_dssp.ps',device=postscript,horizontal=F)
                                        #}

  ############ for R.R.
  par(mar=c(5,6,3,1))
  par(mfrow=c(1,3))
  par(las=1)
  par(cex=0.75)
  otypes <- c("ALA.HB","LEU.HB","PHE.HB","TYR.HB","ARG.HB","ASP.HB",NA,"ALA.OB","LEU.OB","PHE.OB","PRO.OB","TYR.OB","ARG.OB","ASP.OB",NA,"TYR.HP","TYR.OH","CYS.HP","TRP.HP","ASN.HP","ASN.OA","ARG.HP")
 
  barplot(rbind(pp.uf.ss[,'E'][otypes],pp.uf.ss[,'H'][otypes],pp.uf.ss[,'L'][otypes]),
          beside=T,col=c('green','red','blue'),space=c(0.2,1),horiz=T,xpd=F,xlim=c(0,0.5),
          main="Native Proteins",xlab="Fraction Unsatisfied")
  barplot(rbind(dp.uf.ss[,'E'][otypes],dp.uf.ss[,'H'][otypes],dp.uf.ss[,'L'][otypes]),
          beside=T,col=c('green','red','blue'),space=c(0.2,1),horiz=T,xpd=F,xlim=c(0,0.5),
          main="Rosetta Decoys",xlab="Fraction Unsatisfied")
  barplot(rbind(ur.ss[,'E'][otypes],ur.ss[,'H'][otypes],ur.ss[,'L'][otypes]),
          beside=T,col=c('green','red','blue'),space=c(0.2,1),horiz=T,xpd=F,
          main="Ratio: Decoy / Native",xlab="Fraction Unsatisfied")
  lines(c(1,1),c(0,10000),col='black',lwd=2)
  dev.copy2eps(file="/users/sheffler/project/research_reports/fig/unsat_polar_by_ss.eps")


  uf <- tapply(pp$hbe>-0.01,pp$pdb,sum) / tapply(pp$pdb,pp$pdb,length)

  uf <- 1:length(levels(pp$pdb))
  names(uf) <- levels(pp$pdb)
  for(ii in 1:length(levels(pp$pdb))) {
    print( ii )
    tmp <- pp$hbe[pp$pdb==levels(pp$pdb)[ii]] >= -0.01
    uf[ levels(pp$pdb)[ii] ] <- sum(tmp)/length(tmp)
    print ( uf[ levels(pp$pdb)[ii] ] ) 
  }
  sol <- tapply(pdbe$sole,pdbe$pdb,mean)
  n <- names(sol)[ names(uf) %in% names(sol) ]

  par(mfrow=c(1,1))  
  plot(sol[n],1-uf[n],ylim=c(0.7,1),xlim=c(1,2.5),
       xlab="Lazaridis-Karplus Solvation Energy",
       ylab="Buried Unsatisfied Polar Score",
       main="Bur. Unsat. Polar Score vs. Solvation Energy")
  dev.copy2eps(file="/users/sheffler/project/research_reports/fig/bup_vs_solvation.eps")

  
######### ss nad burialbins ###########

  for(types in c('sc','bbo','bbh')) {

                                        #par(mfrow=c(1,2))
                                        #par(cex=0.75)
    
    barplot(rbind(pp.uc.ssbr[,'E.BR'][get(types)],pp.uc.ssbr[,'E.SF'][get(types)],
                  pp.uc.ssbr[,'H.BR'][get(types)],pp.uc.ssbr[,'H.SF'][get(types)],
                  pp.uc.ssbr[,'L.BR'][get(types)],pp.uc.ssbr[,'L.SF'][get(types)]),
            col=c('yellow','lightyellow','red','pink','blue','lightblue'),
            space=c(0.4,2),horiz=T,xlim=c(0,500),beside=T,xpd=F,
            main=paste("PDB",types,"Buried Polar UNS Counts by SSBR
           (blue=loop,red=helix,yellow=sheet)"))
    barplot(rbind(dp.uc.ssbr[,'E.BR'][get(types)],dp.uc.ssbr[,'E.BR'][get(types)],
                  dp.uc.ssbr[,'H.BR'][get(types)],dp.uc.ssbr[,'H.SF'][get(types)],
                  dp.uc.ssbr[,'L.BR'][get(types)],dp.uc.ssbr[,'L.SF'][get(types)]),
            col=c('yellow','lightyellow','red','pink','blue','lightblue'),
            space=c(0.4,2),horiz=T,xlim=c(0,500),beside=T,xpd=F,
            main=paste("DS2",types,"Buried Polar UNS Counts by SSBR
           (blue=loop,red=helix,yellow=sheet)"))
    dev.print(file=paste('fig/050921/type_counts_pdb_ds2_aa_ssbr_',types,'_dssp.ps',sep=''),
              device=postscript,horizontal=F)
    
    barplot(rbind(pp.uf.ssbr[,'E.BR'][get(types)],pp.uf.ssbr[,'E.SF'][get(types)],
                  pp.uf.ssbr[,'H.BR'][get(types)],pp.uf.ssbr[,'H.SF'][get(types)],
                  pp.uf.ssbr[,'L.BR'][get(types)],pp.uf.ssbr[,'L.SF'][get(types)]),
            col=c('yellow','lightyellow','red','pink','blue','lightblue'),
            space=c(0.4,2),horiz=T,beside=T,xpd=F,
            main=paste("DS2",types,"Buried Polar UNS Frac by SSBR
           (blue=loop,red=helix,yellow=sheet)"))
    barplot(rbind(dp.uf.ssbr[,'E.BR'][get(types)],dp.uf.ssbr[,'E.SF'][get(types)],
                  dp.uf.ssbr[,'H.BR'][get(types)],dp.uf.ssbr[,'H.SF'][get(types)],
                  dp.uf.ssbr[,'L.BR'][get(types)],dp.uf.ssbr[,'L.SF'][get(types)]),
            col=c('yellow','lightyellow','red','pink','blue','lightblue'),
            space=c(0.4,2),horiz=T,beside=T,xpd=F,
            main=paste("DS2",types," Buried Polar UNS Frac by SSBR
           (blue=loop,red=helix,yellow=sheet)"))
    dev.print(file=paste('fig/050921/unsfrac_pdb_ds2_aa_ssbr_',types,'_dssp.ps',sep=''), 
              device=postscript,horizontal=F)

    
    barplot(rbind(ur.ssbr[,'E.BR'][get(types)],ur.ssbr[,'E.SF'][get(types)],
                  ur.ssbr[,'H.BR'][get(types)],ur.ssbr[,'H.SF'][get(types)],
                  ur.ssbr[,'L.BR'][get(types)],ur.ssbr[,'L.SF'][get(types)]),
            col=c('yellow','lightyellow','red','pink','blue','lightblue'),
            space=c(0.4,2),horiz=T,beside=T,xlim=c(0,5),xpd=F,
            main=paste("DS2/PDB",types,"Buried Polar UNS Frac Ratio DS2/PDB by SSBR
           (blue=loop,red=helix,yellow=sheet)"))
    lines(c(1,1),c(0,10000),col='purple')
    pc.ssbr <- pmin(pp.pc.ssbr,dp.pc.ssbr)
    barplot(rbind(pc.ssbr[,'E.BR'][get(types)],pc.ssbr[,'E.SF'][get(types)],
                  pc.ssbr[,'H.BR'][get(types)],pc.ssbr[,'H.SF'][get(types)],
                  pc.ssbr[,'L.BR'][get(types)],pc.ssbr[,'L.SF'][get(types)]),
            col=c('yellow','lightyellow','red','pink','blue','lightblue'),
            space=c(0.4,2),horiz=T,beside=T,xpd=F,
            main=paste("min(PDB,DS2)",types,"Buried Polar UNS Counts by SSBR
           (blue=loop,red=helix,yellow=sheet)"))
    dev.print(file=paste('fig/050921/unsfrac_ratio_ds2_pdb_aa_ssbr_',types,'_dssp.ps',sep=''),
              device=postscript,horizontal=F)
    
  }
}
