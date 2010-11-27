load('~/project/features/buried_polar/rdata/dp.decoyset2.dssp.polar.rdata')
load('~/project/features/buried_polar/rdata/pp.pdb.polar.rdata')

source("~/scripts/buried_polar/buried_polar_util.R")


######### Buried UNS protein level discrimination

if(0){
factors <- list(c('type'),c('type','ss'),c('type','br'),c('type','ss','br'))
bth <- 340
pp.br <- c();
pp.br[pp$anb >bth] <- "BR";
pp.br[pp$anb<=bth] <- "SF"; 
pp.br <- as.factor(pp.br)
dp.br <- c();
dp.br[dp$anb >bth] <- "BR";
dp.br[dp$anb<=bth] <- "SF"; 
dp.br <- as.factor(dp.br)
dp$br <- dp.br
pp$br <- pp.br
rm(dp.br,pp.br,bth)

dbin  <- calc.joint.bins(dp,factors,list(c('type','aa')))
pbin  <- calc.joint.bins(pp,factors,list(c('type','aa')))

ncross <- 3
pp.train <- list()
dp.train <- list()
if(ncross==1){
  pp.train[[1]] <- T
  dp.train[[1]] <- T
} else {
  ppr <- runif(dim(pp)[1])
  dpr <- runif(dim(dp)[1])  
  for( ii in 1:ncross ) {
    pp.train[[ii]] <- !( (ii-1)/ncross < ppr & ppr <= (ii)/ncross )
    dp.train[[ii]] <- !( (ii-1)/ncross < dpr & dpr <= (ii)/ncross )
  }
}
pp.pred.pp <- list(); dp.pred.pp <- list()
pp.pred.dp <- list(); dp.pred.dp <- list()
pp.pred.lr <- list(); dp.pred.lr <- list()
for(bin in names(dbin)) {
  pp.pred.pp[[bin]] <- rep(NA,dim(pp)[1]); dp.pred.pp[[bin]] <- rep(NA,dim(dpb)[1]);
  pp.pred.dp[[bin]] <- rep(NA,dim(pp)[1]); dp.pred.dp[[bin]] <- rep(NA,dim(dp)[1]);
  pp.pred.lr[[bin]] <- rep(NA,dim(pp)[1]); dp.pred.lr[[bin]] <- rep(NA,dim(dp)[1]);
}
ppuf <- list(); dpuf <- list(); lruf <- list()
for( ii in 1:ncross ) {
  ppuf[[ii]] <- calc.frac.unsatisfied( unsatisfied=pp$hbe==0 , bins=pbin , pp.train[[ii]] )
  dpuf[[ii]] <- calc.frac.unsatisfied( unsatisfied=dp$hbe==0 , bins=dbin , dp.train[[ii]] )
  lruf[[ii]] <- list()
  for(bin in names(dbin))
    lruf[[ii]][[bin]] <- ppuf[[ii]][[bin]]/dpuf[[ii]][[bin]]
  for(bin in names(dbin)) {
    print(paste("computing preds for",ii,bin))
    pp.test <- !pp.train[[ii]]
    dp.test <- !dp.train[[ii]]
    pp.pred.pp[[bin]][pp.test] <- calc.pred(pbin[[bin]],ppuf[[ii]][[bin]],satisfied=pp$hbe!=0,test=(pp.test))
    pp.pred.dp[[bin]][pp.test] <- calc.pred(pbin[[bin]],dpuf[[ii]][[bin]],satisfied=pp$hbe!=0,test=(pp.test))
    pp.pred.lr[[bin]][pp.test] <- calc.pred(pbin[[bin]],lruf[[ii]][[bin]],satisfied=pp$hbe!=0,test=(pp.test))
    dp.pred.pp[[bin]][dp.test] <- calc.pred(dbin[[bin]],ppuf[[ii]][[bin]],satisfied=dp$hbe!=0,test=(dp.test))
    dp.pred.dp[[bin]][dp.test] <- calc.pred(dbin[[bin]],dpuf[[ii]][[bin]],satisfied=dp$hbe!=0,test=(dp.test))
    dp.pred.lr[[bin]][dp.test] <- calc.pred(dbin[[bin]],lruf[[ii]][[bin]],satisfied=dp$hbe!=0,test=(dp.test))
  }
}
save(pp.pred.pp,pp.pred.dp,pp.pred.lr,dp.pred.pp,dp.pred.dp,dp.pred.lr,
     file='~/project/features/discrimination/rdata/predictions.rdata')

pp.pdb.pred.pp <- list() ; dp.pdb.pred.pp <- list()
pp.pdb.pred.dp <- list() ; dp.pdb.pred.dp <- list()
pp.pdb.pred.lr <- list() ; dp.pdb.pred.lr <- list()
for(bin in c('nobins','aa','aabr','aass','aassbr')) {
  print(paste("computing pdb pred for",bin))
  pp.pdb.pred.pp[[bin]] <- tapply(pp.pred.pp[[bin]],pp$pdb,mean)
  pp.pdb.pred.dp[[bin]] <- tapply(pp.pred.dp[[bin]],pp$pdb,mean)
  pp.pdb.pred.lr[[bin]] <- tapply(pp.pred.lr[[bin]],pp$pdb,mean)
  dp.pdb.pred.pp[[bin]] <- tapply(dp.pred.pp[[bin]],dp$pdb,mean)
  dp.pdb.pred.dp[[bin]] <- tapply(dp.pred.dp[[bin]],dp$pdb,mean)
  dp.pdb.pred.lr[[bin]] <- tapply(dp.pred.lr[[bin]],dp$pdb,mean)
}

pdb.roc.pp <- calc.roc.categories(pp.pdb.pred.pp, dp.pdb.pred.pp,100,T)
pdb.roc.dp <- calc.roc.categories(pp.pdb.pred.dp, dp.pdb.pred.dp,100,T)
pdb.roc.lr <- calc.roc.categories(pp.pdb.pred.lr, dp.pdb.pred.lr,100,T)


########## polar group level discrimination
pp.pred.unsat.pp <- list(); dp.pred.unsat.pp <- list();
pp.pred.unsat.dp <- list(); dp.pred.unsat.dp <- list();
pp.pred.unsat.lr <- list(); dp.pred.unsat.lr <- list(); 
pp.unsat <- pp$hbe==0
dp.unsat <- dp$hbe==0
for(bin in names(dbin)) {
  pp.pred.unsat.pp[[bin]] <- pp.pred.pp[[bin]][pp.unsat]
  pp.pred.unsat.dp[[bin]] <- pp.pred.dp[[bin]][pp.unsat]
  pp.pred.unsat.lr[[bin]] <- pp.pred.lr[[bin]][pp.unsat]
  dp.pred.unsat.pp[[bin]] <- dp.pred.pp[[bin]][dp.unsat]
  dp.pred.unsat.dp[[bin]] <- dp.pred.dp[[bin]][dp.unsat]
  dp.pred.unsat.lr[[bin]] <- dp.pred.lr[[bin]][dp.unsat]
}

roc.pp <- calc.roc.categories(pp.pred.unsat.pp, dp.pred.unsat.pp,100,T)
roc.dp <- calc.roc.categories(pp.pred.unsat.dp, dp.pred.unsat.dp,100,T)
roc.lr <- calc.roc.categories(pp.pred.unsat.lr, dp.pred.unsat.lr,100,T)

par(cex=0.5)
par(mar=c(3,3,3,1))
par(mfrow=c(2,3))
plot.roc(pdb.roc.pp,"PP Buried UNS Protein Discrimination ")
plot.roc(pdb.roc.dp,"DP Buried UNS Protein Discrimination ")
plot.roc(pdb.roc.lr,"PDB/DS2 Ratio Buried UNS Protein Discrimination")

plot.roc(roc.pp,"PP Buried UNS Atom Discrimination ")
plot.roc(roc.dp,"DP Buried UNS Atom Discrimination ")
plot.roc(roc.lr,"PDB/DS2 Ratio Buried UNS Atom Discrimination")

file <- "~/project/features/discrimination/fig/05sept26_roc_pdb_v_decoy_buried_uns_bins.ps"
dev.print(file=file,device=postscript,horizontal=T)       


}

######### residue level discrimination

if(0) {
source("~/scripts/buried_polar/buried_polar_util.R")
setwd("~/project/features/discrimination/")

load('~/project/features/buried_polar/rdata/dp.decoyset2.dssp.polar.rdata')
load('~/project/features/buried_polar/rdata/pp.pdb.polar.rdata')
pp.pdb  <- pp$pdb
pp.prot <- pp$prot
pp.res  <- pp$res
dp.pdb  <- dp$pdb
dp.prot <- dp$prot
dp.res  <- dp$res
save(pp.pdb,file='rdata/pp.pdb.rdata')
save(dp.pdb,file='rdata/dp.pdb.rdata')
rm(pp,dp)

#load('~/project/features/discrimination/rdata/predictions.rdata')
#save(pp.pred.lr,dp.pred.lr,file='rdata/pred.lr.rdata')

load('rdata/pred.lr.rdata')

factors <- "aa" #names(dp.pred.pp)
#pp.res.pred.pp <- list(); dp.res.pred.pp <- list()
#pp.res.pred.dp <- list(); dp.res.pred.dp <- list()
pp.res.pred.lr <- list(); dp.res.pred.lr <- list()
for(bin in factors) {
  pp.res.pred.pp[[bin]] <- c()
}
for(ii in 1:length(levels(dp.prot))) {
  I <- dp.prot==levels(dp.prot)[ii]
  lbl <- paste(dp.pdb[I],dp.res[I])
  for(bin in factors) {
    print(paste("making bp.res.pred",ii,bin))
    #tmp1 <- tapply(dp.pred.pp[[bin]][I] , lbl , mean )
    #tmp2 <- tapply(dp.pred.dp[[bin]][I] , lbl , mean )
    tmp3 <- tapply(dp.pred.lr[[bin]][I] , lbl , mean )
    #dp.res.pred.pp[[bin]] <- c(dp.res.pred.pp[[bin]],tmp1)
    #dp.res.pred.dp[[bin]] <- c(dp.res.pred.dp[[bin]],tmp2)
    dp.res.pred.lr[[bin]] <- c(dp.res.pred.lr[[bin]],tmp3)
  }
}
     for(ii in 1:length(levels(pp.pdb))) {
  I <- pp.prot==levels(pp.prot)[ii]
  lbl <- paste(pp.pdb[I],pp.res[I])
  for(bin in factors) {
    print(paste("making bpp.res.pred",ii,bin))
    #tmp1 <- tapply(pp.pred.pp[[bin]][I] , lbl , mean )
    #tmp2 <- tapply(pp.pred.dp[[bin]][I] , lbl , mean )
    tmp3 <- tapply(pp.pred.lr[[bin]][I] , lbl , mean )
    #pp.res.pred.pp[[bin]] <- c(pp.res.pred.pp[[bin]],tmp1)
    #pp.res.pred.dp[[bin]] <- c(pp.res.pred.dp[[bin]],tmp2)
    pp.res.pred.lr[[bin]] <- c(pp.res.pred.lr[[bin]],tmp3)
  }
}
save(pp.res.pred.pp,dp.res.pred.pp,
     file='~/project/features/discrimination/rdata/res.pred.pp.rdata')
save(pp.res.pred.dp,dp.res.pred.dp,
     file='~/project/features/discrimination/rdata/res.pred.dp.rdata')
save(pp.res.pred.lr,dp.res.pred.lr,
     file='~/project/features/discrimination/rdata/res.pred.lr.rdata')

}

if(0) {

source("~/scripts/buried_polar/buried_polar_util.R")
setwd("~/project/features/discrimination/")

load('~/decoy_energy_polar/pe.pdb.energies.rdata')
load('~/decoy_energy_polar/de.decoyset2.energy.rdata')
load('~/project/features/discrimination/rdata/res.pred.lr.rdata')
load('~/project/features/discrimination/rdata/pred.lr.rdata')
load('rdata/pp.pdb.rdata')
load('rdata/dp.pdb.rdata')

pp.res.pred <- pp.res.pred.lr$aa
pp.names <- names(pp.res.pred)
pe.names <- paste(pe$pdb,pe$res)
print(paste("len pp", length(pp.names),
            "len pe", length(pe.names),
            "pp in pe",sum(pp.names %in% pe.names),
            "pe in pp",sum(pe.names %in% pp.names)
))
einp  <- pe.names %in% pp.names
extra <- rep(1,sum(!einp))
names(extra) <- pe.names[!einp]
pp.res.pred  <- c(pp.res.pred,extra)
pp.names <- names(pp.res.pred)
opp <- order(pp.names)
ope <- order(pe.names)
pp.res.pred <- pp.res.pred[opp][order(ope)]
if(any(names(pp.res.pred)!=pe.names))
  print("PP names wrong!")


dp.res.pred <- dp.res.pred.lr$aa
dp.names <- names(dp.res.pred)
de.names <- paste(de$pdb,de$res)
print(paste("len dp", length(dp.names),
            "len de", length(de.names),
            "dp in de",sum(dp.names %in% de.names),
            "de in dp",sum(de.names %in% dp.names)
))
einp  <- de.names %in% dp.names
extra <- rep(1,sum(!einp))
names(extra) <- de.names[!einp]
dp.res.pred  <- c(dp.res.pred,extra)
dp.names <- names(dp.res.pred)
odp <- order(dp.names)
ode <- order(de.names)
dp.res.pred <- dp.res.pred[odp][order(ode)]
if(any(names(dp.res.pred)!=de.names))
  print("DP names wrong!")


pe.bpe.lr.aa <- pp.res.pred 
de.bpe.lr.aa <- dp.res.pred 
save(pe.bpe.lr.aa,file='rdata/pe.bpe.lr.aa')
save(de.bpe.lr.aa,file='rdata/de.bpe.lr.aa')

load('rdata/pe.bpe.lr.aa')
load('rdata/de.bpe.lr.aa')

}

if(0) {

enames <- c("atre","repe","sole","probe","dune","intrae","hbe","paire","rese","tlje")

roc <- list()
#for(e in enames) {
#  roc[[e]] <- calc.roc(-pe[[e]],-de[[e]])
#}
roc$sasapack <- calc.roc(-pe$sasapack,             -de$sasapack)
roc$sasaprob <- calc.roc( pe$sasaprob,              de$sasaprob)
roc$bpe  <-     calc.roc( pe.bpe.lr.aa,             de.bpe.lr.aa)
roc$bpe.times.sasa <- calc.roc( pe$sasaprob*  pe.bpe.lr.aa, de$sasaprob*  de.bpe.lr.aa)
roc$bpe.plus.sasa  <- calc.roc( pe$sasaprob+3*pe.bpe.lr.aa, de$sasaprob+3*de.bpe.lr.aa)
pdbscore <- pe$tlje + pe$sole + pe$probe + pe$intrae + pe$hbe + pe$paire
decscore <- de$tlje + de$sole + de$probe + de$intrae + de$hbe + de$paire
roc$score <- calc.roc( -pdbscore, -decscore )
roc$everything <- calc.roc( -pdbscore+pe$sasaprob+3*pe.bpe.lr.aa, -decscore+de$sasaprob+3*de.bpe.lr.aa)

pe.pdb <- list()
pp.pdb.names <- unique(pp.pdb)
#roc <- list()
for(e in enames) {
  pe.pdb[[e]]      <- tapply(-pe[[e]],      pe$pdb,mean)[pp.pdb.names]
}
pe.pdb$sasapack <- tapply(-pe$sasapack,  pe$pdb,mean)[pp.pdb.names]
pe.pdb$sasaprob <- tapply(pe$sasaprob,   pe$pdb,mean)[pp.pdb.names]
pe.pdb$bpe      <- tapply(pp.pred.lr$aa, pp.pdb,mean)[pp.pdb.names]
pe.pdb$bpe.times.sasa  <-   pe.pdb$bpe * pe.pdb$sasaprob
pe.pdb$bpe.plus.sasa   <- 3*pe.pdb$bpe + pe.pdb$sasaprob
pe.pdb$score <- pe.pdb$tlje + pe.pdb$sole + pe.pdb$probe + pe.pdb$intrae + pe.pdb$hbe + pe.pdb$paire
pe.pdb$everything   <- pe.pdb$score + pe.pdb$bpe.plus.sasa

depe.pdb$score <- pe.pdb$tlje + pe.pdb$sole + pe.pdb$probe + pe.pdb$intrae + pe.pdb$hbe + pe.pdb$paire
.pdb <- list()
dp.pdb.names <- unique(dp.pdb)
#roc <- list()
for(e in enames) {
  de.pdb[[e]]      <- tapply(-de[[e]],de$pdb,mean)[dp.pdb.names]
}
de.pdb$sasapack <- tapply(-de$sasapack,  de$pdb,mean)[dp.pdb.names]
de.pdb$sasaprob <- tapply(de$sasaprob,   de$pdb,mean)[dp.pdb.names]
de.pdb$bpe      <- tapply(dp.pred.lr$aa, dp.pdb,mean)[dp.pdb.names]
de.pdb$bpe.times.sasa  <-   de.pdb$bpe * de.pdb$sasaprob
de.pdb$bpe.plus.sasa   <- 3*de.pdb$bpe + de.pdb$sasaprob
de.pdb$score <- de.pdb$tlje + de.pdb$sole + de.pdb$probe + de.pdb$intrae + de.pdb$hbe + de.pdb$paire
de.pdb$everything   <- de.pdb$score + de.pdb$bpe.plus.sasa

pdbroc <- calc.roc.categories(pe.pdb,de.pdb,100)


par(mfrow=c(2,2))
par(mar=c(3,3,3,1))
par(cex=0.5)

I <- runif(length(pe$hbe)) < 10000/length(pe$hbe)
plot(pe$hbe[I],     pe.bpe.lr.aa[I],pch='.',main="Buried Polar E vs HBe")
plot(pe.pdb$sasaprob,pe.pdb$bpe       , xlim=c(0.2,0.8),ylim=c(0.8,1.1), pch=1,
     main="SASAprob vs Buried Polar E",xlab="Sasa Prob",ylab="Buried Polar Score")
points(de.pdb$sasaprob,de.pdb$bpe   ,pch='.'    ,col=2)

par(mfrow=c(2,1))
par(cex=0.5)
plot.roc(roc,    title="Residue level discrimination",legx=0.6,legy=0.4)
plot.roc(pdbroc, title="Protein level discrimination",legx=0.6,legy=0.4)
#file <- "fig/05sept28_res_and_pdb_discrimination_on_energies.ps"
#dev.print(file=file,device=postscript,horizontal=F)       


#file <- "fig/05sept28_res_and_pdb_discrimination_on_score_sasa_bpe.ps"
#file <- "test.ps"
#dev.print(file=file,device=postscript,horizontal=T)       


tmp <- list()
for(n in c("bpe.plus.sasa","everything","score","bpe","sasapack","tlje"))
  tmp[[n]] <- pdbroc[[n]]

plot.roc(tmp,title="PDB vs. Decoy Library Discrimination",legx=0.22,legy=0.35,
         colors=c("purple","black","darkgreen","blue","red","green"),
         xlim=c(0,0.5))
dev.copy2eps(file="/users/sheffler/project/research_reports/fig/ROC_pdb_v_decoy_whole_prot.eps")




}
