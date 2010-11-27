load('~/project/features/buried_polar/rdata/dp.decoyset2.dssp.polar.rdata')
load('~/project/features/buried_polar/rdata/pp.pdb.polar.rdata')

source("~/scripts/buried_polar/buried_polar_util.R")


plot.roc <- function(roc,title="Some ROC Curves",colors=NULL,file=NULL) {
  par(cex=0.5)
  rocscores <- format(round(as.numeric(lapply(roc,function(x) x$score)),3))
  if(is.null(colors))
    colors = rainbow(length(roc))
  plot (c(0,1),c(0,1) , col=1 , type='n',main=title,
        xlab="Fraction False Positives",ylab="Fraction True Positives")
  for(ii in 1:length(roc)) 
    lines(roc[[ii]], col=colors[ii])
  legend(.50,.25,c("ROC", paste(rocscores,'bin on:',names(roc))), col=c(1,colors),
         lwd=c(0,rep(3,length(roc))))
  if(!is.null(file))
    dev.print(file=file,device=postscript,horizontal=F)       
}


######### Buried UNS protein level discrimination

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
  pp.pred.pp[[bin]] <- rep(NA,dim(pp)[1]); dp.pred.pp[[bin]] <- rep(NA,dim(dp)[1]);
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

pdb.roc.pp <- calc.roc.categories(pp.pdb.pred.pp, dp.pdb.pred.pp,10,T)
pdb.roc.dp <- calc.roc.categories(pp.pdb.pred.dp, dp.pdb.pred.dp,10,T)
pdb.roc.lr <- calc.roc.categories(pp.pdb.pred.lr, dp.pdb.pred.lr,10,T)


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

roc.pp <- calc.roc.categories(pp.pred.unsat.pp, dp.pred.unsat.pp,10,T)
roc.dp <- calc.roc.categories(pp.pred.unsat.dp, dp.pred.unsat.dp,10,T)
roc.lr <- calc.roc.categories(pp.pred.unsat.lr, dp.pred.unsat.lr,10,T)

par(mfrow=c(2,3))
plot.roc(pdb.roc.pp,"PP Buried UNS Protein Discrimination ")
plot.roc(pdb.roc.dp,"DP Buried UNS Protein Discrimination ")
plot.roc(pdb.roc.lr,"PDB/DS2 Ratio Buried UNS Protein Discrimination")

plot.roc(roc.pp,"PP Buried UNS Atom Discrimination ")
plot.roc(roc.dp,"DP Buried UNS Atom Discrimination ")
plot.roc(roc.lr,"PDB/DS2 Ratio Buried UNS Atom Discrimination")

file <- "~/project/features/05sept26_roc_pdb_v_decoy_buried_uns_bins.ps"
dev.print(file=file,device=postscript,horizontal=T)       




######### residue level discrimination

load('~/decoy_energy_polar/pe.pdb.energies.rdata')
load('~/decoy_energy_polar/de.decoyset2.energy.rdata')

I <- 1:10000
pe.pdbres <- as.factor(paste(as.character(pp$pdb[I]),as.character(pp$res[I])))
de.pdbres <- as.factor(paste(as.character(dp$pdb[I]),as.character(dp$res[I])))
pp.pdbres <- as.factor(paste(as.character(pp$pdb[I]),as.character(pp$res[I])))
dp.pdbres <- as.factor(paste(as.character(dp$pdb[I]),as.character(dp$res[I])))

pp.res <- tapply(pp.pred.lr$aassbr, pp.pdbres, mean)
dp.res <- tapply(dp.pred.lr$aassbr, dp.pdbres, mean)



