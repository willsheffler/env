load('pdbpolar.rdata')

aanames1 <- c("A","R","N","D","C","Q","E","G","H","I",
              "L","K","M","F","P","S","T","W","Y","V")
aanames3 <- c("ALA","ARG","ASN","ASP","CYS","GLN","GLU","GLY","HIS","ILE",
              "LEU","LYS","MET","PHE","PRO","SER","THR","TRP","TYR","VAL")


resname2rosettaid <- 1:20
names(resname2rosettaid) <- c("ALA","CYS","ASP","GLU","PHE","GLY","HIS",
                              "ILE","LYS","LEU","MET","ASN","PRO","GLN",
                              "ARG","SER","THR","VAL","TRP","TYR")

Npolar <- dim(pdbpolar)[1]
Npoints <- 20000
cut <- runif(Npolar) < Npoints / Npolar
isH <- pdbpolar$atomtype==22 | pdbpolar$atomtype==25
isHcut <- isH[cut]

###############################################
### burial measure densities

par(mfrow=c(4,2))
par(las=2)

fracs <- c()
fracs[1] <- sum(pdbpolar$sasa14==0)/Npolar
fracs[2] <- sum(pdbpolar$sasa14[isH]==0) / sum(isH)
fracs[3] <- sum(pdbpolar$sasa14[!isH]==0) / sum(!isH)
fracs[4] <- sum(pdbpolar$sasa10==0)/Npolar
fracs[5] <- sum(pdbpolar$sasa10[isH]==0) / sum(isH)
fracs[6] <- sum(pdbpolar$sasa10[!isH]==0) / sum(!isH)
names(fracs) <- c("frac buried sasa14",
               "frac H buried sasa14",
               "frac A buried sasa14",
               "frac buried sasa10",
               "frac H buried sasa10",
               "frac A buried sasa10"     )
barplot(fracs, col=c(1,2,3,1,2,3))


d <- density(pdbpolar$sasa14)
do <- density(pdbpolar$sasa14[!isH])
dh <- density(pdbpolar$sasa14[ isH])
plot(d,main="SASA14 histogram",xlim=c(0,10))
plot(do,main="SASA14 histogram (red=H)",xlim=c(0,10),
     ylim=c(0,max(do$y,dh$y)), col=1)
lines(dh,col=2)

d <- density(pdbpolar$lk)
do <- density(pdbpolar$lk[!isH])
dh <- density(pdbpolar$lk[ isH])
plot(d,main="LK histogram")
plot(do,main="LK histogram (red=H)", ylim=c(0,max(do$y,dh$y)), col=1)
lines(dh,col=2)

d <- density(pdbpolar$lk2)
do <- density(pdbpolar$lk2[!isH])
dh <- density(pdbpolar$lk2[ isH])
plot(d,main="LK2 histogram")
plot(do,main="LK2 histogram (red=H)", ylim=c(0,max(do$y,dh$y)), col=1)
lines(dh,col=2)

d <- density(pdbpolar$lk4)
do <- density(pdbpolar$lk4[!isH])
dh <- density(pdbpolar$lk4[ isH])
plot(d,main="LK4 histogram")
plot(do,main="LK4 histogram (red=H)", ylim=c(0,max(do$y,dh$y)), col=1)
lines(dh,col=2)

dev.print(device=postscript,file='../fig/burial_measures_densities.ps',horizontal=F)

################################################
### burial measure comparisons

par(mfrow=c(3,2))
plot(pdbpolar$sasa14[cut]+runif(sum(cut))*0.4 ,
     pdbpolar$lk[cut] ,  pch='.', col=isHcut+1 ,
     main="LK vs sasa14 (red=hydrogen)",
     xlab="sasa14",  ylab="LK")
do <- density(pdbpolar$lk[pdbpolar$sasa14==0 & !isH])
dh <- density(pdbpolar$lk[pdbpolar$sasa14==0 & isH])
plot(do, col=1, ylim=c(0,max(do$y,dh$y)),
     main="LK density for sasa14==0 (red=hydrogen)",xlab="LK burial")
lines(dh,col=2)

plot(pdbpolar$sasa14[cut]+runif(sum(cut))*0.4 ,
     pdbpolar$lk2[cut] ,
     pch='.', col=isHcut+1 ,   main="LK2 vs sasa 14 (red=hydrogen)",
     xlab="sasa14",     ylab="LK2")
do <- density(pdbpolar$lk2[pdbpolar$sasa14==0 & !isH])
dh <- density(pdbpolar$lk2[pdbpolar$sasa14==0 & isH])
plot(do, col=1, ylim=c(0,max(do$y,dh$y)),
     main="LK2 density for sasa14==0 (red=hydrogen)",xlab="LK2 burial")
lines(dh,col=2)

plot(pdbpolar$sasa14[cut]+runif(sum(cut))*0.4 ,
     pdbpolar$lk4[cut] ,
     pch='.',  col=isHcut+1 ,    main="LK4 vs sasa 14 (red=hydrogen)",
     xlab="sasa14",    ylab="LK4")
do <- density(pdbpolar$lk4[pdbpolar$sasa14==0 & !isH])
dh <- density(pdbpolar$lk4[pdbpolar$sasa14==0 & isH])
plot(do, col=1, ylim=c(0,max(do$y,dh$y)),
     main="LK4 density for sasa14==0 (red=hydrogen)",xlab="LK4 burial")
lines(dh,col=2)

dev.print(device=postscript,file='../fig/burial_measures_vs_sasa14.ps',horizontal=F)


par(mfrow=c(2,2))
plot(pdbpolar$sasa14[cut]+runif(sum(cut))*0.4 ,
     pdbpolar$sasa10[cut]+runif(sum(cut))*0.3 ,
     pch='.',main="sasa10 vs sasa 14", col=isHcut+1 ,
     xlab="sasa14",ylab="sasa10")
plot(pdbpolar$lk[cut] ,
     pdbpolar$lk2[cut] ,
     pch='.',main="LK vs LK2", col=isHcut+1 ,
     xlab="LK",ylab="LK2")
plot(pdbpolar$lk[cut] ,
     pdbpolar$lk4[cut] ,
     pch='.',main="LK vs LK2", col=isHcut+1 ,
     xlab="LK",ylab="LK4")
plot(pdbpolar$lk2[cut] ,
     pdbpolar$lk4[cut] ,
     pch='.',main="LK2 vs LK4", col=isHcut+1 ,
     xlab="LK2",ylab="LK4")
dev.print(device=postscript,file='../fig/burial_measures_scatter1.ps',horizontal=F)



plot(pdbpolar$nb[cut]+runif(sum(cut))*1.0 ,
     pdbpolar$sasa14[cut]+runif(sum(cut))*0.4 ,
     pch='.',
     #ain="sasa10 vs sasa 14", col=isHcut+1 ,
     #lab="sasa14",ylab="sasa10"
     )

plot(pdbpolar$nb[cut]+runif(sum(cut))*1.0 ,
     pdbpolar$lk[cut]+runif(sum(cut)) ,
     pch='.',
     #ain="sasa10 vs sasa 14", col=isHcut+1 ,
     #lab="sasa14",ylab="sasa10"
     )

plot(pdbpolar$nb[cut]+runif(sum(cut))*1.0 ,
     pdbpolar$lk2[cut]+runif(sum(cut)) ,
     pch='.',
     #ain="sasa10 vs sasa 14", col=isHcut+1 ,
     #lab="sasa14",ylab="sasa10"
     )

plot(pdbpolar$nb[cut]+runif(sum(cut))*1.0 ,
     pdbpolar$lk4[cut]+runif(sum(cut)) ,
     pch='.',
     #ain="sasa10 vs sasa 14", col=isHcut+1 ,
     #lab="sasa14",ylab="sasa10"
     )
