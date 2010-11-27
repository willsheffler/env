load('nativebp.rdata')
load('decoybp.rdata')
load('pdbbp.rdata')

Ndecoy  <- dim(decoybp)[1]
Nnative <- dim(nativebp)[1]
Npdb    <- dim(pdbbp)[1]


dsets <- c('decoybp','nativebp','pdbbp')

anbthresh <- quantile(decoybp$anb10,0:19/20)
nbthresh  <- 1:30
hbethresh <- quantile(decoybp$hbe,1:20/20)

burieduns <- array(0,c(2,length(dsets),30,20,20))
dimnames(burieduns) <- list(c('frac','total'),
                            dsets,
                            paste('nb',1:30),
                            paste('anb',names(anbthresh)),
                            paste('hbe',names(hbethresh))
                            )


nbburied  <- function(dset,nbth)  get(dset)$nb    >= nbth
anbburied <- function(dset,anbth) get(dset)$anb10 >= anbth
fhbeuns   <- function(dset,hbeth) get(dset)$hbe   >= hbeth
hbeth <- -0.01

for(dset in dsets) {  
  for(nbth in 5:30) {
    Nburieduns <- sum(nbburied(dset,nbth)&fhbeuns(dset,hbeth))
    burieduns['frac', dset,nbth,1] <- Nburieduns/sum(nbburied(dset,nbth))
    burieduns['total',dset,nbth,1] <- Nburieduns
  }
  for(anbth in 1:20) {
    Nburieduns <- sum(anbburied(dset,anbthresh[anbth])&fhbeuns(dset,hbeth))
    burieduns['frac', dset,1,1]     <- Nburieduns/sum(anbburied(dset,anbthresh[anbth]))
    burieduns['total',dset,1,anbth] <- Nburieduns
  }
}

##################3

dsets <- c('decoybp','nativebp','pdbbp')
types <- sort(unique(decoybp$atomtype))

nbthresh  <- 1:30
anbthresh <- quantile(decoybp$anb10,0:19/20)
hbethresh <- quantile(decoybp$hbe,1:20/20)

nballfrac  <- array(0,c(length(dsets),30))
anballfrac <- array(0,c(length(dsets),20))
hbeallfrac <- array(0,c(length(dsets),20))
nballtot   <- array(0,c(length(dsets),30))
anballtot  <- array(0,c(length(dsets),20))
hbealltot  <- array(0,c(length(dsets),20))
dimnames(nballfrac)  <- list(dsets,paste('nb',1:30))
dimnames(anballfrac) <- list(dsets,paste('anb',names(anbthresh)))
dimnames(hbeallfrac) <- list(dsets,paste('hbe',names(hbethresh)))
dimnames(nballtot)   <- list(dsets,paste('nb',1:30))
dimnames(anballtot)  <- list(dsets,paste('anb',names(anbthresh)))
dimnames(hbealltot)  <- list(dsets,paste('hbe',names(hbethresh)))

sasaburied <- function(dset)       rep(T,length(get(dset)$nb))
nbburied   <- function(dset,nbth)  get(dset)$nb    >= nbth
anbburied  <- function(dset,anbth) get(dset)$anb10 >= anbth
fhbeuns    <- function(dset,hbeth) get(dset)$hbe   >= hbeth


for(dset in dsets) {
  hbeth <- -0.01
  for(ii in 5:30) {    
    nbth <- ii
    Nbu <- sum( nbburied(dset,nbth) & fhbeuns(dset,hbeth) )
    nballfrac[ dset,ii] <- Nbu/sum( nbburied(dset,nbth) )
    nballtot[dset,ii] <- Nbu
  }
  for(ii in 1:20) {
    anbth <- anbthresh[ii]
    hbeth <- -0.01
    Nbu <- sum( anbburied(dset,anbth) & fhbeuns(dset,hbeth) )
    anballfrac[ dset,ii] <- Nbu/sum( anbburied(dset,anbth) )
    anballtot[dset,ii] <- Nbu
  }
  for(ii in 1:20) {
    hbeth <- hbethresh[ii]
    Nbu <- sum( sasaburied(dset) & fhbeuns(dset,hbeth) )
    hbeallfrac[ dset,ii] <- Nbu/sum(sasaburied(dset))
    hbealltot[dset,ii] <- Nbu
  }
}

#nbbutot[nbbutot==0] <- 1
#anbbutot[anbbutot==0] <- 1
#hbebutot[hbebutot==0] <- 1

par(mfrow=c(3,3))
par(mar=c(4,3,4,1))
par(mgp=c(2,1,0))

mx <- max(nballfrac,anballfrac)
plot(5:30,5:30,type='n',ylim=c(0,mx),xlab="nbthresh",ylab="buried uns frac",main="NB bu frac")
for(ii in 1:length(dsets)) lines(5:30,nballfrac[dsets[ii], 5:30],col=ii)
plot(anbthresh,anbthresh,type='n',ylim=c(0,mx),ylab="buried uns frac",main="ANB bu frac")
for(ii in 1:length(dsets)) lines(anbthresh,anballfrac[dsets[ii], ],col=ii)
plot(hbethresh,hbethresh,type='n',ylim=c(0,mx),ylab="buried uns frac",main="HBE bu frac")
for(ii in 1:length(dsets)) lines(hbethresh,hbeallfrac[dsets[ii], ],col=ii)

rnative <- nballfrac['decoybp',]/nballfrac['nativebp',]
rpdb    <- nballfrac['decoybp',]/nballfrac['pdbbp',]
mx <- 3#max(rnative,rpdb,na.rm=T)
plot(5:30,5:30,type='n',ylim=c(0,mx), xlab="nbthresh",ylab="comp. buried uns frac",
     main="NB Comp Frac")
lines(5:30,rnative[5:30],col=2)
lines(5:30,rpdb[5:30],col=3)

rnative <- anballfrac['decoybp',]/anballfrac['nativebp',]
rpdb    <- anballfrac['decoybp',]/anballfrac['pdbbp',]
#mx <- max(rnative,rpdb,na.rm=T)
plot(anbthresh,anbthresh,type='n',ylim=c(0,mx),ylab="comp. buried uns frac",
     main="ANB Comp Frac")
lines(anbthresh,rnative,col=2)
lines(anbthresh,rpdb,   col=3)

rnative <- hbeallfrac['decoybp',]/hbeallfrac['nativebp',]
rpdb    <- hbeallfrac['decoybp',]/hbeallfrac['pdbbp',]
#mx <- max(rnative,rpdb,na.rm=T)
plot(hbethresh,hbethresh,type='n',ylim=c(0,mx),ylab="comp. buried uns frac",
     main="HBE Comp Frac")
lines(hbethresh,rnative,col=2)
lines(hbethresh,rpdb,   col=3)

mx <- max(nballtot[,])
plot(5:30,5:30,type='n',ylim=c(10,mx), log='y',xlab="buried uns count",
     main="NB Samp. Count")
for(ii in 1:length(dsets)) lines(5:30,nballtot[dsets[ii], 5:30],col=ii)

mx <- max(anballtot[,])
plot(anbthresh,anbthresh,type='n',ylim=c(10,mx), log='y',xlab="buried uns count",
     main="ANB Samp. Count")
for(ii in 1:length(dsets)) lines(anbthresh,anballtot[dsets[ii], ],col=ii)

mx <- max(hbealltot[,])
plot(hbethresh,rep(1000,20),type='n',ylim=c(10,mx), log='y',xlab="buried uns count",
     main="HBE Samp. Count")
for(ii in 1:length(dsets)) lines(hbethresh,hbealltot[dsets[ii], ],col=ii)


dev.print(file="../fig/unsfrac_decoy_native_pdb.ps",horizontal=F,device=postscript)


###########################################
# unsfrac vs anb burial by atomtype
##########################################

plotbu <- function(bufrac,butot,atomtype,thresh,mx1=NA,mx2=NA,mx3=NA) {
  if(is.na(mx1)) mx1 <- max(bufrac)
  plot(thresh,thresh,type='n',ylim=c(0,mx1),ylab="buried uns frac",
       main=paste("unsfrac type",atomtype))
  for(jj in 1:length(dsets))
    lines(thresh,bufrac[dsets[jj], ],col=jj)
  
  rnative <- bufrac['decoybp',]/bufrac['nativebp',]
  rpdb    <- bufrac['decoybp',]/bufrac['pdbbp',]
  if(is.na(mx2)) mx2 <- max(rnative,rpdb,na.rm=T)
  plot(thresh,thresh,type='n',ylim=c(0,mx2),ylab="comp. buried uns frac",
       main=paste("comp frac",atomtype))
  lines(thresh,rnative,col=2)
  lines(thresh,rpdb,   col=3)
  
  if(is.na(mx3)) mx3 <- max(butot)
  plot(thresh,thresh,type='n',ylim=c(10,mx3), log='y',xlab="buried uns count",
       main=paste("count type",atomtype))
  for(jj in 1:length(dsets))
    lines(thresh,butot[dsets[jj], ],col=jj)
}


dsets <- c('decoybp','nativebp','pdbbp')
types <- sort(unique(decoybp$atomtype))

anbthresh <- quantile(decoybp$anb10,0:19/20)

anbbufrac <- array(0,c(length(types),length(dsets),20))
anbbutot  <- array(0,c(length(types),length(dsets),20))
dimnames(anbbufrac) <- list(paste("type",types),dsets,paste('anb',names(anbthresh)))
dimnames(anbbutot)  <- list(paste("type",types),dsets,paste('anb',names(anbthresh)))

sasaburied <- function(dset)       rep(T,length(get(dset)$nb))
nbburied   <- function(dset,nbth)  get(dset)$nb    >= nbth
anbburied  <- function(dset,anbth) get(dset)$anb10 >= anbth
fhbeuns    <- function(dset,hbeth) get(dset)$hbe   >= hbeth

for(dset in dsets) {
  hbeth <- -0.01
  for(ii in 1:20) {
    anbth <- anbthresh[ii]
    hbeth <- -0.01
    Nbu <- tapply(anbburied(dset,anbth)&fhbeuns(dset,hbeth),get(dset)$atomtype,sum)
    Nb  <- tapply(anbburied(dset,anbth),get(dset)$atomtype,sum)
    anbbufrac[,dset,ii] <- Nbu/Nb
    anbbutot[, dset,ii] <- Nbu
  }
}

par(mfcol=c(3,7))
par(mar=c(4,3,4,1))
par(mgp=c(2,1,0))


for(ii in 1:length(types))
  plotbu(anbbufrac[ii,,],anbbutot[ii,,],types[ii],anbthresh)
dev.print(file="../fig/unsfrac_vs_anb10_by_atomtype.ps",device=postscript)

mx1 <- max(anbbufrac)
mx2 <- 5
mx3 <- max(anbbutot)
for(ii in 1:length(types))
  plotbu(anbbufrac[ii,,],anbbutot[ii,,],types[ii],anbthresh,mx1,mx2,mx3)
dev.print(file="../fig/unsfrac_vs_anb10_by_atomtype_same_axes.ps",device=postscript)


###########################################
# unsfrac vs anb burial by atomtype for various pdbres
##########################################

plotbupdbres <- function(bufrac,butot,atomtype,thresh,mx1=NA,mx2=NA,mx3=NA) {
  color <- c(1,2,"#00FF00","#00DD00","#00BB00","#009900","#007700","#005500")

  if(is.na(mx1)) mx1 <- max(bufrac)
  plot(thresh,thresh,type='n',ylim=c(0,mx1),ylab="buried uns frac",
       main=paste("unsfrac type",atomtype))
  for(jj in 1:length(dsets))
    lines(thresh,bufrac[dsets[jj], ],col=color[jj])
  
  rnative <- bufrac['decoybp',]/bufrac['nativebp',]
  rpdb    <- bufrac['decoybp',]/bufrac['pdbbp',]
  rpdb18  <- bufrac['decoybp',]/bufrac['pdbbp18',]
  rpdb16  <- bufrac['decoybp',]/bufrac['pdbbp16',]
  rpdb14  <- bufrac['decoybp',]/bufrac['pdbbp14',]
  rpdb12  <- bufrac['decoybp',]/bufrac['pdbbp12',]
  rpdb10  <- bufrac['decoybp',]/bufrac['pdbbp10',]
  if(is.na(mx2)) mx2 <- max(rnative,rpdb,na.rm=T)
  plot(thresh,thresh,type='n',ylim=c(0,mx2),ylab="comp. buried uns frac",
       main=paste("comp frac",atomtype))
  lines(thresh,rnative,col=2)
  lines(thresh,rpdb,     col=color[3])
  lines(thresh,rpdb18,   col=color[4])
  lines(thresh,rpdb16,   col=color[5])
  lines(thresh,rpdb14,   col=color[6])
  lines(thresh,rpdb12,   col=color[7])
  lines(thresh,rpdb10,   col=color[8])
  
  if(is.na(mx3)) mx3 <- max(butot)
  plot(thresh,thresh,type='n',ylim=c(10,mx3), log='y',xlab="buried uns count",
       main=paste("count type",atomtype))
  for(jj in 1:length(dsets))
    lines(thresh,butot[dsets[jj], ],col=color[jj])
}

pdbbp10 <- pdbbp[pdbresolution<=1.0,]
pdbbp12 <- pdbbp[pdbresolution<=1.2,]
pdbbp14 <- pdbbp[pdbresolution<=1.4,]
pdbbp16 <- pdbbp[pdbresolution<=1.6,]
pdbbp18 <- pdbbp[pdbresolution<=1.8,]

dsets <- c('decoybp','nativebp','pdbbp','pdbbp18',
           'pdbbp16','pdbbp14','pdbbp12','pdbbp10')
types <- sort(unique(decoybp$atomtype))

anbthresh <- quantile(decoybp$anb10,0:19/20)

anbbufrac <- array(0,c(length(types),length(dsets),20))
anbbutot  <- array(0,c(length(types),length(dsets),20))
dimnames(anbbufrac) <- list(paste("type",types),dsets,paste('anb',names(anbthresh)))
dimnames(anbbutot)  <- list(paste("type",types),dsets,paste('anb',names(anbthresh)))

sasaburied <- function(dset)       rep(T,length(get(dset)$nb))
nbburied   <- function(dset,nbth)  get(dset)$nb    >= nbth
anbburied  <- function(dset,anbth) get(dset)$anb10 >= anbth
fhbeuns    <- function(dset,hbeth) get(dset)$hbe   >= -0.01

for(dset in dsets) {
  print(dset)
  hbeth <- -0.01
  for(ii in 1:20) {
    anbth <- anbthresh[ii]
    hbeth <- -0.01
    Nbu <- tapply(anbburied(dset,anbth)&fhbeuns(dset,hbeth),get(dset)$atomtype,sum)
    Nb  <- tapply(anbburied(dset,anbth),get(dset)$atomtype,sum)
    anbbufrac[,dset,ii] <- Nbu/Nb
    anbbutot[, dset,ii] <- Nbu
  }
}

par(mfcol=c(3,7))
par(mar=c(4,3,4,1))
par(mgp=c(2,1,0))


for(ii in 1:length(types))
  plotbupdbres(anbbufrac[ii,,],anbbutot[ii,,],types[ii],anbthresh)
dev.print(file="../fig/unsfrac_vs_anb10_various_pdbres_by_atomtype.ps",device=postscript)

mx1 <- max(anbbufrac)
mx2 <- 5
mx3 <- max(anbbutot)
for(ii in 1:length(types))
  plotbupdbres(anbbufrac[ii,,],anbbutot[ii,,],types[ii],anbthresh,mx1,mx2,mx3)
dev.print(file="../fig/unsfrac_vs_anb10_various_pdbres_by_atomtype_same_axes.ps",device=postscript)


##########################################3
# various hbond def
###########################################

dsets <- c('decoybp','nativebp','pdbbp')
types <- sort(unique(decoybp$atomtype))

hbdthresh <- 0:19/19*1.5+2.5

hbdallbufrac <- array(0,c(length(dsets),20))
hbdallbutot  <- array(0,c(length(dsets),20))
dimnames(hbdallbufrac) <- list(dsets,paste('hbD',hbdthresh))
dimnames(hbdallbutot)  <- list(dsets,paste('hbD',hbdthresh))
hbdbufrac <- array(0,c(length(types),length(dsets),20))
hbdbutot  <- array(0,c(length(types),length(dsets),20))
dimnames(hbdbufrac) <- list(paste("type",types),dsets,paste('hbD',hbdthresh))
dimnames(hbdbutot)  <- list(paste("type",types),dsets,paste('hbD',hbdthresh))

sasaburied <- function(dset)       rep(T,length(get(dset)$nb))
hbduns     <- function(dset,hbdth) get(dset)$dist >= hbdth

for(dset in dsets) {
  for(ii in 1:20) {
    hbdth <- hbdthresh[ii]
    Nbu <- tapply(sasaburied(dset)&hbduns(dset,hbdth),get(dset)$atomtype,sum)
    Nb  <- tapply(sasaburied(dset),get(dset)$atomtype,sum)
    hbdbufrac[,dset,ii] <- Nbu/Nb
    hbdbutot[, dset,ii] <- Nbu
    hbdallbufrac[dset,ii] <- sum(Nbu)/sum(Nb)
    hbdallbutot[dset,ii]  <- sum(Nbu)
  }
}

par(mfcol=c(3,8))
par(mar=c(4,3,4,1))
par(mgp=c(2,1,0))

plotbu(hbdallbufrac,hbdallbutot,'all',hbdthresh)
for(ii in 1:length(types))
  plotbu(hbdbufrac[ii,,],hbdbutot[ii,,],types[ii],hbdthresh)
dev.print(file="../fig/unsfrac_vs_hbond_dist_by_atomtype.ps",device=postscript)

mx1 <- max(hbdbufrac)
mx2 <- 3
mx3 <- max(hbdbutot)
plotbu(hbdallbufrac,hbdallbutot,'all',hbdthresh,mx1,mx2,mx3)
for(ii in 1:length(types))
  plotbu(hbdbufrac[ii,,],hbdbutot[ii,,],types[ii],hbdthresh,mx1,mx2,mx3)
dev.print(file="../fig/unsfrac_vs_hbond_dist_by_atomtype_same_axes.ps",device=postscript)


###########################################
## resolution effects
##########################################

tmp <- pdbinfo$resolution
names(tmp) <- pdbinfo$IDs
pdbresolution <- tmp[pdbbp$pdb]

resthresh <- 0:19/19+1
pdbresallfrac <- array(0,c(length(dsets),length(resthresh)))
pdbresalltot  <- array(0,c(length(dsets),length(resthresh)))
pdbresbufrac <- array(0,c(length(types),length(dsets),length(resthresh)))
pdbresbutot  <- array(0,c(length(types),length(dsets),length(resthresh)))
dimnames(pdbresallfrac) <- list(dsets,NULL)
dimnames(pdbresalltot)  <- list(dsets,NULL)
dimnames(pdbresbufrac)  <- list(types,dsets,NULL)
dimnames(pdbresbutot)   <- list(types,dsets,NULL)

Ndecoy  <- dim(decoybp)[1]
Nnative <- dim(nativebp)[1]
pdbresallfrac['decoybp',]  <- sum(decoybp$hbe >= -0.01)/Ndecoy
pdbresbufrac [,'decoybp',] <- tapply(decoybp$hbe >= -0.01,decoybp$atomtype,function(x) sum(x)/length(x))
pdbresallfrac[ 'nativebp',] <- sum(nativebp$hbe >= -0.01)/Nnative
pdbresbufrac [,'nativebp',] <- tapply(nativebp$hbe >= -0.01,nativebp$atomtype,function(x) sum(x)/length(x))
pdbresalltot[ 'decoybp',] <- sum(decoybp$hbe >= -0.01)
pdbresbutot [,'decoybp',] <- tapply(decoybp$hbe >= -0.01,decoybp$atomtype,sum)
pdbresalltot[ 'nativebp',] <- sum(nativebp$hbe >= -0.01)
pdbresbutot [,'nativebp',] <- tapply(nativebp$hbe >= -0.01,nativebp$atomtype,sum)
Iuns <- pdbbp$hbe >= -0.01
for(ii in 1:length(resthresh)) {
  resth <- resthresh[ii]
  Ires <- pdbresolution <= resth
  Nb  <- tapply(Ires,       pdbbp$atomtype,sum)
  Nbu <- tapply(Ires & Iuns,pdbbp$atomtype,sum)
  pdbresbufrac[,'pdbbp',ii] <- Nbu/Nb
  pdbresbutot [,'pdbbp',ii] <- Nbu
  pdbresallfrac['pdbbp',ii] <- sum(Nbu)/sum(Nb)
  pdbresalltot ['pdbbp',ii] <- sum(Nbu)
}

par(mfcol=c(3,8))
par(mar=c(4,3,4,1))
par(mgp=c(2,1,0))

plotbu(pdbresallfrac,pdbresalltot,'all',resthresh)
for(ii in 1:length(types))
  plotbu(pdbresbufrac[ii,,],pdbresbutot[ii,,],types[ii],resthresh)
dev.print(file="../fig/unsfrac_vs_pdbres_by_atomtype.ps",device=postscript)

mx1 <- max(pdbresbufrac)
mx2 <- 3
mx3 <- max(pdbresbutot)


###########################################
#
############################################

d1 <- density(decoybp$nb,1)
d2 <- density(nativebp$nb,1)
d3 <- density(pdbbp$nb,1)
plot(d1)
lines(d2)
lines(d3)


types <- sort(unique(pdbbp$atomtype))
densities5 <- list()
densities75 <- list()
densities10 <- list()
for(ii in 1:length(types)) {
  densities5[[ii]] <- density(pdbbp$anb5[pdbbp$atomtype==types[ii]],bw=1)
  densities75[[ii]] <- density(pdbbp$anb75[pdbbp$atomtype==types[ii]],bw=3)
  densities10[[ii]] <- density(pdbbp$anb10[pdbbp$atomtype==types[ii]],bw=10)
}
mx5  <- max(as.numeric(lapply(densities5, function(x) max(x$y))))
mx75 <- max(as.numeric(lapply(densities75,function(x) max(x$y))))
mx10 <- max(as.numeric(lapply(densities10,function(x) max(x$y))))

par(mfrow=c(3,1))

plot(densities5[[1]],type='n',xlim=c(min(pdbbp$anb5),max(pdbbp$anb5)),ylim=c(0,mx5),
     main="Density of atom nb 5.0 given sasa14==0")
for(ii in 1:length(densities5))
  lines(densities5[[ii]],col=ii)
legend(70,0.07,paste('type',types),col=1:7,lwd=2)

plot(densities75[[1]],type='n',xlim=c(min(pdbbp$anb75),max(pdbbp$anb75)),ylim=c(0,mx75),
     main="Density of atom nb 7.5 given sasa14==0")
for(ii in 1:length(densities75))
  lines(densities75[[ii]],col=ii)

plot(densities10[[1]],type='n',xlim=c(min(pdbbp$anb10),max(pdbbp$anb10)),ylim=c(0,mx10),
     main="Density of atom nb 10.0 given sasa14==0")
for(ii in 1:length(densities10))
  lines(densities10[[ii]],col=ii)


dev.print(device=postscript,file='../fig/burial_atom_nb_by_atomtypetype.ps',horizontal=F)

####################################################
# hbond def play
####################################################

Irand <- runif(dim(pdbbp)[1]) < 10000/dim(pdbbp)[1]
par(mfrow=c(3,2))
I <- Irand
plot(pdbbp$cosa[I],pdbbp$cosh[I],pch='.',xlab="H-A-B cosine",ylab="D-H-A cosine",
     main="Hbond angle parameters")
I <- pdbbp$hbe < -0.01 & Irand
plot(pdbbp$cosa[I],pdbbp$cosh[I],pch='.',xlab="H-A-B cosine",ylab="D-H-A cosine",
     main="Hbond angle parameters for HBE <= -0.01")

I <- Irand
plot(pdbbp$dist[I],pdbbp$cosh[I],pch='.',xlab="H-A dist",ylab="D-H-A cosine",
     main="Hbond dist vs D-H-A angle",xlim=c(1,4))
I <- pdbbp$hbe < -0.01 & Irand
plot(pdbbp$dist[I],pdbbp$cosh[I],pch='.',xlab="H-A dist",ylab="D-H-A cosine",
     main="Hbond dist vs D-H-A angle for HBE <= -0.01",xlim=c(1,4))

I <- Irand
plot(pdbbp$dist[I],pdbbp$cosa[I],pch='.',xlab="H-A dist",ylab="H-A-B cosine",
     main="Hbond dist vs H-A-B angle",xlim=c(1,4))
I <- pdbbp$hbe < -0.01 & Irand
plot(pdbbp$dist[I],pdbbp$cosa[I],pch='.',xlab="H-A dist",ylab="H-A-B cosine",
     main="Hbond dist vs H-A-B angle for HBE <= -0.01",xlim=c(1,4))

dev.print(device=postscript,file='../fig/hbond_angle_dist_scatter.ps',horizontal=F)


####################################################3
# discrimination tests
####################################################
library(gtools,lib.loc='~/software/Rlib')
library(gdata, lib.loc='~/software/Rlib')
library(gplots,lib.loc='~/software/Rlib')
library(ROCR,  lib.loc='~/software/Rlib')

# simple counts

anbthresh <- quantile(decoybp$anb10,0:19/19)
fanbbur <- function(x,anbth) x$anb10 >= anbth
fhbeuns <- function(x,hbeth) x$hbe >= hbeth

fbur <- function(x) fanbbur(x,anbthresh[10])
funs <- function(x) fhbeuns(x, -0.01)


decoypdbbur    <- tapply(fbur(decoybp),              decoybp$pdb,sum)
decoypdbburuns <- tapply(fbur(decoybp)&funs(decoybp),decoybp$pdb,sum)
decoypdbbur <- pmax(1,decoypdbbur)
