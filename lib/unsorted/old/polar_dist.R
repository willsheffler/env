polar.dist.names <- c("pdb","res1","aa1","atom1","at1","atype1",
                            "res2","aa2","atom2","at2","atype2","dist")

pd <- as.data.frame(read.table('~/pdb_energies/pdb_polar_dist_lt5.data'))[,-1:-2]
names(pd) <- polar.dist.names
save(pd,file="~/data/rdata/pd.pdb_polar_dist_lt5.rdata")

dd <- as.data.frame(read.table('decoy_energy_polar/ds2_polar_dist_lt5.data'))[,-1:-2]
names(dd) <- polar.dist.names
save(dd,file="~/data/rdata/dd.ds2_polar_dist_lt5.rdata")

load(file="~/data/rdata/pd.pdb_polar_dist_lt5.rdata")
load(file="~/data/rdata/dd.ds2_polar_dist_lt5.rdata")

pd.22.22.den <- density(pd$dist[pd$at1==22&pd$at2==22],0.050)
dd.22.22.den <- density(dd$dist[dd$at1==22&dd$at2==22],0.050)
plot (pd.22.22.den,col=2)
lines(dd.22.22.den,col=1)

atypes <- unique(pd$at1)
pdden <- list()
ddden <- list()
pn    <- list()
dn    <- list()
for(t1 in atypes ) {
  for(t2 in atypes ) {
    f <- paste(t1,t2,sep='.')
    print(f)
    pn[[f]] <- sum(pd$at1==t1 & pd$at2==t2 & pd$dist < 3)
    dn[[f]] <- sum(dd$at1==t1 & dd$at2==t2 & dd$dist < 3)
    if(pn[[f]] > 0 & dn[[f]] > 0) {
      pdden[[f]]  <- density(pd$dist[ pd$at1==t1 & pd$at2==t2 ], 0.050 )
      ddden[[f]]  <- density(dd$dist[ dd$at1==t1 & dd$at2==t2 ], 0.050 )
      #pdden[[f]]$y <- pdden[[f]]$y / pdden[[f]]$x**2
      #ddden[[f]]$y <- ddden[[f]]$y / ddden[[f]]$x**2
    }
  }
}

XMAX <- 5

par(mfrow=c(6,6))
par(cex=0.5)
par(mar=c(1,1,1,1))
par(mgp=c(1,0,0))
par(tck=0)
for(ii in 1:36) {
  n <- names(pn)[ii]
  if(pn[[n]] > 0 & dn[[n]] > 0) {
    ym <- max(max(pdden[[n]]$y[pdden[[n]]$x<XMAX]),max(ddden[[n]]$y[ddden[[n]]$x<XMAX]))
    plot (pdden[[n]],col=2,main=paste(n,pn[[n]],"/",dn[[n]]),ylim=c(0,ym),xlim=c(1.5,XMAX))
    lines(ddden[[n]],col=1)
  } else {
    plot(1,1)
  }
}

dev.copy2eps(file="~/project/features/polar_dist_lt5_scale3.eps")


polar.dist.names <- c("pdb","res1","aa1","atom1","at1","atype1",
                            "res2","aa2","atom2","at2","atype2","dist")

read.big.data <- function(fname1,fname2,fname3,reader=read.table,
                          takecols=T,colnames=NA,
                          facfunc=NA,facN=99999999) {
  fnames <- paste(fname1,fname2,fname3,sep='')
  for(file in fnames) {
    tmpdata <- reader(file)[,takecols]
    names(tmpdata) <- colnames
    if(!is.na(facfunc)) {
      fac <- facfunc(tmpdata)
      for(f in levels(fac)) {
        II <- fac==f
        if(sum(II) <= facN )
          d <- tmpdata[II,]
        else
          d <- tmpdata[which(II)[runif(sum(II))<facN/sum(II)],]
        if(f==levels(fac)[1])
          d <- x
        else
          d <- rbind(d,x)
      }
    }
    if(file==fnames[1])
      data <- d
    else
      data <- rbind(data,d)
  }
  return(data)
}



select10k <- function(x,N=10000) {
  if(dim(x)[1] < N)
    return(x)
  else
    return(x[ runif(dim(x)[1]) < N/dim(x)[1] ],)
}

II         <- pd$at1 > pd$at2
tmp        <- pd$at1[II]
pd$at1[II] <- pd$at2[II]
pd$at2[II] <- tmp
pairtypes  <- as.factor(paste(as.character(pd$at1),as.character(pd$at2),sep='/'))
pd <- tapply(pd,pairtypes,select10k)
