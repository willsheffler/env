#load(file="~/project/features/discrimination/rdata/bpgsbr.gsbr.buried.polar.rdata")
#load(file="~/project/features/discrimination/rdata/bpgsgr.gsgr.buried.polar.rdata")
#load('~/project/features/buried_polar/rdata/dp.decoyset2.dssp.polar.rdata')
load('~/project/features/buried_polar/rdata/pp.pdb.polar.rdata')

source('~/scripts/R-utils.R')

setwd("~/project/features/buried_polar/")

write.uns.data.file <- function(uf,file) {
  n <- names(uf)
  aa <- substr(n,1,3)
  aai <- res2num[aa]
  at <- substr(n,6,7)
  ss <- substr(n,10,10)
  d <- as.data.frame(list(index=1:length(n),aa=aa,aai=aai,atomtype=at,ss=ss))
  d <- cbind(d,uf)
  write.data.frame(d,file=file,pretty=T,meta=paste("META Nbins",length(uf)))
}

aaat       <- as.factor(paste(pp$aa,'at',format(pp$atomtype),'ss.',sep=''))
aaatss     <- as.factor(paste(pp$aa,'at',format(pp$atomtype),'ss',pp$ss,sep=''))
ppufaaat   <- tapply(pp$hbe==0,aaat,sum)/tapply(pp$type,aaat,length)
ppufaaatss <- tapply(pp$hbe==0,aaatss,sum)/tapply(pp$type,aaatss,length)
write.uns.data.file(ppufaaat,file='unsatisfied_buried_polar__pdb__aa_at.data')
write.uns.data.file(ppufaaatss,file='unsatisfied_buried_polar__pdb__aa_at_ss.data')

if(0){
gsgruf <- tapply(bpgsgr$hbe==0,bpgsgr$type,sum)/tapply(bpgsgr$type,bpgsgr$type,length)
gsbruf <- tapply(bpgsbr$hbe==0,bpgsbr$type,sum)/tapply(bpgsbr$type,bpgsbr$type,length)
gsuf <- gsgruf/gsbruf
save(gsuf,file="~/project/features/discrimination/rdata/gsuf.buried.uns.frac.rdata")


ppuf <- tapply(pp$hbe==0,pp$type,sum)/tapply(pp$type,pp$type,length)
dpuf <- tapply(dp$hbe==0,dp$type,sum)/tapply(dp$type,dp$type,length)
rtuf <- ppuf/dpuf
save(rtuf,file="~/project/features/discrimination/rdata/uf.buried.uxns.frac.rdata")

}
