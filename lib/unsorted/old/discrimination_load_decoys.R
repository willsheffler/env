decoydir <- "/users/sheffler/data/decoys"
decoys <- list()
decoys$bqian_caspbench <- c("t196","t199","t205","t206","t223","t224","t232","t234","t249")
decoys$phil_homolog <-    c("1af7","1csp","1di2","1mky","1n0u","1ogw","1shf","1tig",
                            "1b72","1dcj","1dtj","1mla","1o2f","1r69","1tif","2reb" )
decoys$kira_hom <-        c("1a4p","1agi","1ar0","1b07","1bmb","1cfy","1d2n","1erv","1pva",
                            "1acf","1ahn","1aw2","1be9","1btn","1crb","1dol","1ig5","1qav" )
decoys$sraman_nmr <-      c("1fvk","1hb6","1hoe","1i27","1who","2int","3il8","4rnt")
decoys$kira_bench <-      c("1acp","1b72","1ig5","1r69","1tig","1vii",
                            "1aa3","1ail","1bf4","1pgx","1tif","1utg")#,"256b","1a32")



get.bpe <- function(bpfile,uf) {
  bp <- read.table(bpfile,header=F)[,-1]
  names(bp) <- c('prot','pdb','aa','res','ss','atom','atomtype','hbe','anb')
  type <- c()
  type[bp$atomtype== 8] <- "NH"
  type[bp$atomtype==13] <- "OH"
  type[bp$atomtype==14] <- "OA"
  type[bp$atomtype==15] <- "OC"
  type[bp$atomtype==20] <- "OB"
  type[bp$atomtype==22] <- "HP"
  type[bp$atomtype==25] <- "HB"
  type <- paste(as.character(bp$aa),as.character(type),sep='.')
  type <- as.factor(type)
  bp$type <- type
  levels(type) <- uf
  type <- as.numeric(as.character(type))
  type[bp$hbe!=0] <- 1
  bpe <- tapply(type,bp$pdb,mean)[sort(bp$pdb)]
  #n <- as.character(lapply(strsplit(gsub('//','/',names(bpe)),'/'),function(x) x[length(x)]))
  #names(bpe) <- n
  bpenaive <- tapply(bp$hbe!=0,bp$pdb,mean)[sort(bp$pdb)]
  #names(bpenaive) <- n
  #bp$pdb <- n
  #return(list(bpe=bpe,bpenaive=bpenaive))
  bp <- bp[order(bp$pdb),]
  return(list(bpe=bpe,bpenaive=bpenaive,bpdata=bp))
}

get.energies <- function(efile) {
  e <- read.table(efile,header=F)
  sasaprob <- tapply(e[,22],e[,3],mean)[sort(e[,3])]
  sasapack <- tapply(e[,21],e[,3],mean)[sort(e[,3])]
  #names(sasaprob) <- as.character(lapply(strsplit(gsub('//','/',names(sasaprob)),'/'),function(x) x[length(x)]))
  #names(sasapack) <- as.character(lapply(strsplit(gsub('//','/',names(sasapack)),'/'),function(x) x[length(x)]))
  H <- sum(e[,6]=='H')/length(e[,6])
  E <- sum(e[,6]=='E')/length(e[,6])
  L <- sum(e[,6]=='L')/length(e[,6])
  #return(list(sasaprob=sasaprob,sasapack=sasapack,len=sum(e[,3]==e[1,3]),H=H,E=E,L=L))
  return(list(sasaprob=sasaprob,sasapack=sasapack,len=sum(e[,3]==e[1,3]),H=H,E=E,L=L,edata=e))
}

calc.decoy.scores <- function(uf,categories,get.full.tables=F) {
  decoyscores <- data.frame()
  bpgsbr <- data.frame()
  bpgsgr <- data.frame()
  egsbr <- data.frame()
  egsgr <- data.frame()
  for(set in names(decoys)) {
    for(prot in decoys[[set]]) {
      if(missing(categories)) {
        categories <- dir(paste(decoydir,set,prot,sep='/'),pattern='\.list$')
        categories <- as.character(lapply(categories,function(x)substring(x,1,nchar(x)-5)))
      }
      for(category in categories) {      
        print(paste(set,prot,category))
        file <- paste(decoydir,'/',set,'/',prot,'/data/',category,'_score.fasc',sep='')
        if(!file.exists(file)){
          print(paste('no fasc',set,prot,category))
        } else {
          print('    reading fasc')
          fasc <- read.table(file,header=T,comment.char='')
          if(fasc$filename[1] == paste(prot,'.pdb',sep=''))
            fasc <- fasc[-1,]
          if(dim(fasc)[1]){
            if('rms' %in% names(fasc)) {
              #if(length(bpe)!=dim(fasc)[1])
              #  print(paste("truncated fasc",length(bpe),dim(fasc)[1],set,prot,category))             
              pfile <- paste(decoydir,'/',set,'/',prot,'/data/',category,'_polar.data',sep='')
              efile <- paste(decoydir,'/',set,'/',prot,'/data/',category,'_energy.data',sep='')
              if(!file.exists(pfile)) {
                print(paste("data files don't exist for",set,prot,category))
              } else {
                print('    reading polar')
                bp  <- get.bpe(pfile,uf=uf)
                bpe <- bp$bpe[1:dim(fasc)[1]]
                bpenaive <- bp$bpenaive[1:dim(fasc)[1]]
                print('    reading energy')
                energies  <- get.energies(efile)
                sasaprob <- energies$sasaprob[1:dim(fasc)[1]]
                sasapack <- energies$sasapack[1:dim(fasc)[1]]
                H <- energies$H
                L <- energies$L
                E <- energies$E
                len <- energies$len              
                fasc$pdb <- names(bpe)
                fasc$bpe <- bpe
                fasc$bpenaive <- bpenaive
                fasc$sasaprob <- sasaprob
                fasc$sasapack <- sasapack
                fasc$length <- len
                fasc$set <- set
                fasc$prot <- prot
                fasc$category <- category
                fasc$H <- H
                fasc$E <- E
                fasc$L <- L
                fasc$setprot <- as.factor(paste(set,prot))
                decoyscores <- rbind(decoyscores,fasc)
                if(get.full.tables) {
                  if(category=='gsgr') {
                    bpgsgr <- rbind(bpgsgr,bp$bpdata)
                    egsgr  <- rbind(egsgr,energies$edata)
                  } else if(category=='gsbr') {
                    bpgsbr <- rbind(bpgsbr,bp$bpdata)
                    egsbr  <- rbind(egsbr,energies$edata)
                  }
                }
              }
            } else {
              print(paste('no rms in fasc',set,prot,category))            
            }
          } else {        
            print(paste('no fasc',set,prot,category))            
          }
        }
      }
    }
  }
  decoyscores
  enames <- c("score","env","pair","vdw","hs","ss","sheet","cb","rsigma",
              "hb_srbb","hb_lrbb","rg","co","contact","rama","bk_tot",
              "fa_atr","fa_rep","fa_sol","h2o_sol","hbsc","fa_dun",
              "fa_intra","fa_pair","fa_plane","fa_prob","fa_h2o","h2o_hb",
              "gsolt","sasa")
  for(e in enames)
    decoyscores[[e]] <- decoyscores[[e]]/decoyscores$length

  #save(bpgsbr,file="~/project/features/discrimination/rdata/bpgsbr.gsbr.buried.polar.rdata")
  #save(bpgsgr,file="~/project/features/discrimination/rdata/bpgsgr.gsgr.buried.polar.rdata")
  return(list(ds=decoyscores,bpgsgr=bpgsgr,bpgsbr=bpgsbr,egsgr=rgsgr,egsbr=egsbr))
}

load("~/project/features/discrimination/rdata/gsuf.type.gsgr_over_gsbr.rdata")
load("~/project/features/discrimination/rdata/rtuf.type.pp_over_dp.rdata")

uf <- gsuf
get.full.tables <- T
l <- calc.decoy.scores(uf=rtuf,
                       categories=c('gsbr','gsgr','relax_native'),
                       get.full.tables=get.full.tables)
ds <- l$ds
save(ds,file="~/project/features/discrimination/rdata/ds.ds3.rtuf.decoyscores.rdata")
save(ds,file="~/project/features/discrimination/rdata/ds.ds3.gsuf.decoyscores.rdata")
if(get.full.tables) {
  bpgsgr <- l$bpgsgr
  bpgsbr <- l$bpgsbr
  egsgr  <- l$egsgr
  egsbr  <- l$egsbr
  save(bpgsbr,file="~/project/features/discrimination/rdata/bpgsbr.gsbr.buried.polar.rdata")
  save(bpgsgr,file="~/project/features/discrimination/rdata/bpgsgr.gsgr.buried.polar.rdata")
  save(bpgsbr,file="~/project/features/discrimination/rdata/egsbr.gsbr.energies.rdata")
  save(bpgsgr,file="~/project/features/discrimination/rdata/egsgr.gsgr.energies.rdata")
}
   
#load("~/project/features/discrimination/rdata/gsuf.buried.uns.frac.rdata")
#ds <- calc.decoy.scores(uf=gsuf,categories=c('gsbr','gsgr','relax_native'))
#save(ds,file="~/project/features/discrimination/rdata/ds.ds3.gsuf.decoyscores.rdata")


#load("~/project/features/discrimination/rdata/ds.ds3.gsuf.decoyscores.rdata")
#ds.gs <- ds
#load("~/project/features/discrimination/rdata/ds.ds3.rtuf.decoyscores.rdata")
#ds.rt <- ds

#roc <- list()
#roc$gs <- calc.roc(ds.gs$bpe[ds.gs$category=='gsgr'],ds.gs$bpe[ds.gs$category=='gsbr'])
#roc$rt <- calc.roc(ds.rt$bpe[ds.rt$category=='gsgr'],ds.rt$bpe[ds.rt$category=='gsbr'])
#plot.roc(roc)

if(0) {

load('~/project/features/buried_polar/rdata/dp.decoyset2.dssp.polar.rdata')
load('~/project/features/buried_polar/rdata/pp.pdb.polar.rdata')
source('~/scripts/buried_polar/buried_polar_util.R')

ppuf <- tapply(pp$hbe==0,pp$type,sum)/tapply(pp$type,pp$type,length)
dpuf <- tapply(dp$hbe==0,dp$type,sum)/tapply(dp$type,dp$type,length)
rtuf <- ppuf/dpuf
save(rtuf,file="~/project/features/discrimination/rdata/uf.buried.uns.frac.rdata")

load(file="~/project/features/discrimination/rdata/uf.buried.uns.frac.rdata")
ppred <- pp$type
levels(ppred) <- rtuf
ppred <- as.numeric(as.character(ppred))
ppred[pp$hbe!=0] <- 1
dpred <- dp$type
levels(dpred) <- rtuf
dpred <- as.numeric(as.character(dpred))
dpred[dp$hbe!=0] <- 1

pppdbpred <- tapply(ppred,pp$pdb,mean)
dppdbpred <- tapply(dpred,dp$pdb,mean)

roc <- calc.roc(pppdbpred,dppdbpred)
plot.roc(roc)

}


if(0) {

load(file="~/project/features/discrimination/rdata/bpgsbr.gsbr.buried.polar.rdata")
load(file="~/project/features/discrimination/rdata/bpgsgr.gsgr.buried.polar.rdata")
source('~/scripts/buried_polar/buried_polar_util.R')

#prots <- paste(as.charac)

gsgruf <- tapply(bpgsgr$hbe==0,bpgsgr$type,sum)/tapply(bpgsgr$type,bpgsgr$type,length)
gsbruf <- tapply(bpgsbr$hbe==0,bpgsbr$type,sum)/tapply(bpgsbr$type,bpgsbr$type,length)
gsuf <- gsgruf/gsbruf
save(gsuf,file="~/project/features/discrimination/rdata/gsuf.buried.uns.frac.rdata")

uf <- gsuf
#uf <- rtuf
               
gsgrred <- bpgsgr$type
levels(gsgrred) <- uf
gsgrred <- as.numeric(as.character(gsgrred))
gsgrred[bpgsgr$hbe!=0] <- 1
gsbrred <- bpgsbr$type
levels(gsbrred) <- uf
gsbrred <- as.numeric(as.character(gsbrred))
gsbrred[bpgsbr$hbe!=0] <- 1

gsgrpdbpred <- tapply(gsgrred,bpgsgr$pdb,mean)
gsbrpdbpred <- tapply(gsbrred,bpgsbr$pdb,mean)

roc <- list(calc.roc(gsgrpdbpred,gsbrpdbpred))
plot.roc(roc)

}

d1 <- density(gsgrpdbpred)
d2 <- density(gsgr.gs.bpe)
plot (d1)
lines(d2,col=2)

gsbr <- ds.gs$bpe[ds$category=='gsbr']
gsgr <- ds.gs$bpe[ds$category=='gsgr']
dgsbr1 <- density(gsbr)
dgsgr1 <- density(gsgr)
dgsbr2 <- density(gsbrpdbpred)
dgsgr2 <- density(gsgrpdbpred)

plot (dgsbr1)
lines(dgsgr1,col=2)
lines(dgsbr2,col=3)
lines(dgsgr2,col=4)

gsbr <- ds.rt$bpe[ds$category=='gsbr']
gsgr <- ds.rt$bpe[ds$category=='gsgr']
drtbr1 <- density(gsbr)
drtgr1 <- density(gsgr)
drtbr2 <- density(gsbrpdbpred)
drtgr2 <- density(gsgrpdbpred)

plot (drtbr1)
lines(drtgr1,col=2)
lines(drtbr2,col=3)
lines(drtgr2,col=4)


n <- as.character(lapply(strsplit(gsub('//','/',as.character(bpgsbr$pdb)),'/'),function(x) x[length(x)]))
