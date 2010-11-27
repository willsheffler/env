#data <- read.table('allE_from_pdb.data',row.names=NULL,header=T)
#ljdata <- data[,1:6]
#save(ljdata,file='/users/sheffler/project/packing_index/rdata/ljdata.rdata')
#load('/users/sheffler/project/packing_index/rdata/ljdata.rdata')
#sasadata <- read.table('sasa_14_and_5.data',header=T)
#packdata <- cbind(ljdata,sasadata[,4:5])
#save(packdata,file='/users/sheffler/project/packing_index/rdata/pdbpackdata.rdata')
#packdata <- read.table('/users/sheffler/project/packing_index/packing_stats',header=T)
#save(packdata,file='/users/sheffler/project/packing_index/rdata/pdbpackdata.rdata')

load('/users/sheffler/project/packing_index/rdata/pdbpackdata.rdata')

RES <- c("ALA","CYS","ASP","GLU","PHE","GLY","HIS","ILE","LYS","LEU",
         "MET","ASN","PRO","GLN","ARG","SER","THR","VAL","TRP","TYR")




#number of occurrences binned on nb
gm <- 0
for(aa in RES) {
  m <- max(hist(pdbdata$nb[pdbdata$aa==aa],breaks=5:30,plot=F)$counts)
  gm <- max(gm,m)
}
par(mfrow=c(4,5))
for(aa in RES) {
  idx <- pdbdata$aa==aa
  m <- max(hist(pdbdata$nb[idx],breaks=5:30,plot=F)$counts)
  hist(pdbdata$nb[idx],breaks=5:30,
       main=paste(aa," #samp max ",m,sep=''))
}


# number of occurances binned on sasafrac
# many residues are buried
par(mfrow=c(4,5))
for(aa in RES) {
  idx <- pdbdata$aa==aa
  b <- seq(from=0,to=1,by=0.05)
  m <- max(hist(pdbdata$sasafrac[idx],breaks=b,plot=F)$counts)
  hist(pdbdata$sasafrac[idx],breaks=b,
       main=paste(aa," #samp max ",m,sep=''))
}


# histograms of LJ scores for residue aa
#for(aa in 1:20) {
  aa <- "MET"
  lb <- floor(min(pdbdata$tlj[pdbdata$aa==aa]))
  ub <- ceiling(max(pdbdata$tlj[pdbdata$aa==aa]))
  histrange = seq(from=lb,to=ub,by=0.5)

  par(mar=c(3,2,2,1))
  par(mfrow=c(5,6))
  #  par(mfg=c(1,5,6,5))
  for(nb in 5:30) {
    idx <- pdbdata$aa==aa & pdbdata$nb==nb
    num <- sum( idx )
    hist( pdbdata$tlj[idx] , breaks=histrange,
         main=paste(aa," nb",nb,' ',num,sep=''))
  }
  #fname = paste("/users/sheffler/project/packing_index/fig/tlj_hist_res_",RES[aa],".eps",sep='')
  #dev.copy2eps(file=fname)
#}



# histograms of LJ scores for all residues with nb neighbors
#for(nb in 5:30) {
  nb <- 20
  lb <- floor(min(pdbdata$tlj[pdbdata$nb==nb]))
  ub <- ceiling(max(pdbdata$tlj[pdbdata$nb==nb]))
  histrange = seq(from=lb,to=ub,by=0.5)

  par(mar=c(3,2,2,1))
  par(mfrow=c(4,5))
  #  par(mfg=c(1,5,6,5))
  for(aa in RES) {
    idx <- pdbdata$aa==aa & pdbdata$nb==nb
    num <- sum( idx )
    hist( pdbdata$tlj[idx] , breaks=histrange,
         main=paste(aa," nb",nb,' ',num,sep=''))
  }
  #fname = paste("/users/sheffler/project/packing_index/fig/tlj_hist_nb_",nb,".eps",sep='')
  #dev.copy2eps(file=fname)
#}



# scatter plot of LJ energy vs. # neighbors for residue aa
par(mfrow=c(4,5))
for(aa in RES) {
#  aa = "LYS"
#  par(mfrow=c(1,1))
  idx <- pdbdata$aa==aa
  c <- cor(pdbdata$nb[idx],pdbdata$tlj[idx])
  plot( pdbdata$nb[idx] +runif(sum(idx))*1.0-0.5 , pdbdata$tlj[idx] , pch='.',
       main=paste(aa,"LJ vs. #nb cor",format(c,digits=2)))
}


# sasapack summary showing sasa5-sasa14 vs sasa14 for all residues
par(mfrow=c(4,5))
par(mar=c(3,2,2,1))
sasapack <- pdbdata$sasa5 - pdbdata$sasa14
for(aa in RES) {
  #aa <- "ALA"
  #par(mfrow=c(1,1))
  idx <- pdbdata$aa==aa
  plot(pdbdata$sasa14[idx],sasapack[idx],pch='.',
       xlim=c(0,60),ylim=c(0,60),
       main=paste(aa," sasapack ",sum(idx),sep=''))
}




# sasapack histograms per residue
aa <- "LEU"
par(mfrow=c(4,5))
sasapack <- pdbdata$sasa5 - pdbdata$sasa14
idx1 <- pdbdata$aa==aa
sasa14 <- pdbdata$sasa14[idx1]
for(bin in 1:20) {
  idx2 <- sasa14 > bin*3-3 & sasa14 <= bin*3
  tmp <- sasapack[idx1][idx2]
  tmp[tmp>60] <- 60
  hist( tmp , breaks=seq(from=0,to=60,by=3),
       main=paste(aa," bin",3*bin-3,'-',3*bin," ",sum(idx2),sep=''))
}


sasapacknsamp <- array(0,c(20,20))
sasapackmean <- array(0,c(20,20))
sasapacksd   <- array(0,c(20,20))
# output sasapack statistics
for(ii in 1:20) {
  aa <- RES[ii]
  idx1 <- pdbdata$aa==aa
  sasa14 <- pdbdata$sasa14[idx1]
  for(bin in 1:20) {
    idx2 <- sasa14 > bin*3-3 & sasa14 <= bin*3    
    sasapacknsamp[ii,bin] <- sum(idx2)
    n <- sum(idx2)
    if(sum(idx2 > 0)) {
      m <- mean(pdbdata$sasa5[idx1][idx2])
      s <- sd(pdbdata$sasa5[idx1][idx2])
    } else {
      m <- NA
      s <- NA
    }
    sasapacknsamp[ii,bin] <- n
    sasapackmean[ii,bin]  <- m
    sasapacksd[ii,bin]    <- s
    print(paste(aa,bin,n,m,s))
  }
}


# sasapack histograms by bin
lb <- 6
ub <- 9
par(mfrow=c(4,5))
sasapack <- pdbdata$sasa5 - pdbdata$sasa14
idx1 <- pdbdata$sasa14 > lb & pdbdata$sasa14 <= ub
aatmp <- pdbdata$aa[idx1]
for(aa in RES) {
  idx2 <- aatmp==aa
  tmp <- sasapack[idx1][idx2]
  tmp[tmp>60] <- 60
  hist( tmp , breaks=seq(from=0,to=60,by=3),
       main=paste(aa," bin",lb,'-',ub," ",sum(idx2),sep=''))
}


# sasafrac vs neighbors
n <- length(pdbdata$nb)
plot(pdbdata$nb+runif(n),pdbdata$sasafrac+runif(n)/100,pch='.')




# scatter plot of LJ energy vs. sasafrac for residue aa
par(mfrow=c(5,4))
par(mar=c(3,2,2,1))
for(aa in RES) {
  #par(mfrow=c(1,1))
  #aa = "LEU"
  idx <- pdbdata$aa==aa
  c <- cor(pdbdata$sasafrac[idx],pdbdata$tlj[idx])
  plot( pdbdata$sasafrac[idx] +runif(sum(idx))*0.01 , pdbdata$tlj[idx] , pch='.',
       main=paste(aa," LJ vs. SASA cor",format(c,digits=2)),xlim=c(0,.8))
}


# scatter plot of LJ energy vs. sasafrac for residue aa for better resolution structs
par(mfrow=c(4,5))
par(mar=c(3,2,2,1))
for(aa in RES) {
  #par(mfrow=c(1,1))
  #aa = "LEU"
  resolution <- pdbinfo$resolution[pdbdata$iter]
  idx <- pdbdata$aa==aa & resolution < 1.5
  c <- cor(pdbdata$sasafrac[idx],pdbdata$tlj[idx])
  plot( pdbdata$sasafrac[idx] +runif(sum(idx))*0.01 , pdbdata$tlj[idx] , pch='.',
       main=paste(aa," LJ vs. SASA cor",format(c,digits=2)),xlim=c(0,.02))
}


# scatter plot of ATR LJ energy vs. sasafrac for residue aa
par(mfrow=c(4,5))
par(mar=c(3,2,2,1))
for(aa in RES) {
  #par(mfrow=c(1,1))
  #aa = "LEU"
  idx <- pdbdata$aa==aa 
  c <- cor(pdbdata$sasafrac[idx],pdbdata$atr[idx])
  plot( pdbdata$sasafrac[idx] +runif(sum(idx))*0.01 , pdbdata$atr[idx] , pch='.',
       main=paste(aa," LJ vs. SASA cor",format(c,digits=2)),xlim=c(0,0.02))
  
}



# scatter plot of LJ energy vs. sasa10
# sasa10 gives more resolution near 0
par(mfrow=c(5,4))
par(mar=c(3,2,2,1))
for(aa in RES) {
  #par(mfrow=c(1,1))
  #aa = "LEU"
  idx <- pdbdata$aa==aa
  c <- cor(pdbdata$sasa10[idx],pdbdata$tlj[idx])
  l <- quantile(pdbdata$sasa10[idx],probs=c(0.90))
  plot( pdbdata$sasa10[idx] , pdbdata$tlj[idx] , pch='.',
       main=paste(aa," LJ vs. SASA cor",format(c,digits=2)) ,xlim=c(0,l))
}





# compute correlation between TLJ and nb/sasa
nbsasacor <- array(0,c(20,2))
rownames(nbsasacor) <- RES
colnames(nbsasacor) <- c("nbcor","sasacor")
for(ii in 1:20) {
  aa <- RES[ii]
  idx <- pdbdata$aa==aa
  nbsasacor[ii,1] <- cor(pdbdata$nb[idx],pdbdata$tlj[idx])
  nbsasacor[ii,2] <- cor(pdbdata$sasafrac[idx],pdbdata$tlj[idx])
}


par(mfrow=c(4,5))
par(mar=c(3,2,2,1))
#break up sasa==0 bin distribution plots
for (aa in RES) {
  idx1 <- pdbdata$sasafrac==0 & pdbdata$aa==aa
  nbth <- mean(pdbdata$nb[idx1])
  idx2 <- pdbdata$nb[idx1] < nbth
  plot(idx2 + runif(length(idx2))*0.9 , pdbdata$tlj[idx1], pch='.',
       main=paste(aa,"nbth",format(nbth,digits=1)))
}

#break up sasa==0 bin density plots
for (aa in RES) {
  idx1 <- pdbdata$sasafrac==0 & pdbdata$aa==aa
  nbth <- mean(pdbdata$nb[idx1])
  idx2 <- pdbdata$nb[idx1] < nbth
  d1 <- density(pdbdata$tlj[idx1][idx2])
  d2 <- density(pdbdata$tlj[idx1][!idx2])
  plot(d1, main=paste(aa,
# LJ vs. sasa05 for sasa14==0
par(mfrow=c(4,5))
par(mar=c(3,2,2,1))
idx1 <- pdbdata$sasa14==0
for(aa in RES) {
  idx2 <- pdbdata$aa[idx1]==aa
  a <- pdbdata$sasa5[idx1][idx2]
  b <- pdbdata$atr[idx1][idx2]
  plot( a , b , pch='.', main=paste(aa,format(cor(a,b),digits=2)))
}
"nbth",format(nbth,digits=3)))
  lines(d2,col=2)
}



#find which buried lysines have low LJ scores
w <- which( pdbdata$aa=="LEU" & pdbdata$sasafrac==0 & pdbdata$tlj > -1.5 )
pdbs <- as.character(pdbinfo$IDs[pdbdata$iter[w]])
aa  <- as.character(pdbdata$aa[w])
res <- as.character(pdbdata$res[w])
rep <- as.character(pdbdata$rep[w])
print(cbind(pdbs,aa,res,rep))

#find which buried lysines have low LJ ATR scores
w <- which( pdbdata$aa=="LEU" & pdbdata$sasafrac==0 & pdbdata$atr > -1.5 )
pdbs <- as.character(pdbinfo$IDs[pdbdata$iter[w]])
aa  <- as.character(pdbdata$aa[w])
res <- as.character(pdbdata$res[w])
rep <- as.character(pdbdata$rep[w])
print(cbind(pdbs,aa,res,rep))

#find which buried lysines have high LJ scores
w <- which( pdbdata$aa=="LEU" & pdbdata$sasafrac==0 & pdbdata$tlj < -6.5 )
pdbs <- as.character(pdbinfo$IDs[pdbdata$iter[w]])
aa  <- as.character(pdbdata$aa[w])
res <- as.character(pdbdata$res[w])
rep <- as.character(pdbdata$rep[w])
print(cbind(pdbs,aa,res,rep))



# LJ vs. sasa05 for sasafrac==0
par(mfrow=c(4,5))
par(mar=c(3,2,2,1))
idx1 <- pdbdata$sasafrac==0
#idx1 <- pdbdata$sasafrac > 0.45 & pdbdata$sasafrac < 0.50
for(aa in RES) {
  idx2 <- pdbdata$aa[idx1]==aa
  a <- pdbdata$sasa5[idx1][idx2]
  b <- pdbdata$atr[idx1][idx2]
  plot( a , b , pch='.', main=paste(aa,format(cor(a,b),digits=2)))
}


# LJ vs. sasa05 for various sasafrac
aa <- 'ALA'
idx1 <- pdbdata$aa==aa
par(mfrow=c(4,5))
par(mar=c(3,2,2,1))
for(bin in 1:20) {
  idx2 <- pdbdata$sasafrac[idx1] > 0.05*bin-0.05 & pdbdata$sasafrac[idx1] < 0.05*bin
  a <- pdbdata$sasa5[idx1][idx2]
  b <- pdbdata$atr[idx1][idx2]
  plot( a , b , pch='.',
       main=paste(aa,format(cor(a,b),digits=2), ' bin ',0.05*bin-0.05,'-',0.05*bin))
}



# LJ vs. sasa05 for sasafrac==0 with decoys in red
par(mfrow=c(5,4))
par(mar=c(3,2,2,1))
idx1 <- pdbdata$sasafrac==0
#idx1 <- pdbdata$sasafrac > 0.45 & pdbdata$sasafrac < 0.50
ik1 <- kiradata$sasafrac==0
for(aa in RES) {
  idx2 <- pdbdata$aa[idx1]==aa
  a <- pdbdata$sasa5[idx1][idx2]
  b <- pdbdata$atr[idx1][idx2]
  plot( a , b , pch='.', main=paste(aa,format(cor(a,b),digits=2)))

  ik2 <- kiradata$aa[ik1]==aa
  a <- kiradata$sasa5[ik1][ik2]
  b <- kiradata$atr[ik1][ik2]
  points( a , b , pch=1, col=2 )

}



# plot density of LJ for pdb vs. decoys
par(mfrow=c(4,5))
par(mar=c(3,2,2,1))
for(aa in RES) {
  idxpdb <- pdbdata$aa==aa
  dljpdb <- density( pdbdata$tlj[idxpdb] )
  idxkira <- kiradata$aa==aa
  dljkira <- density( kiradata$tlj[idxkira] )
  m <- max(max(dljpdb$y),max(dljkira$y))
  plot(dljpdb,ylim=c(0,m),
       main=paste(aa," samp ",sum(idxkira)))
  lines(dljkira,col=2)
}


# plot density of sasa5 for pdb vs. decoys
par(mfrow=c(4,5))
par(mar=c(3,2,2,1))
for(aa in RES) {
  idxpdb <- pdbdata$aa==aa & pdbdata$sasafrac==0
  dljpdb <- density( pdbdata$sasa5[idxpdb] )
  idxkira <- kiradata$aa==aa & kiradata$sasafrac==0
  dljkira <- density( kiradata$sasa5[idxkira] )
  ym <- max(max(dljpdb$y),max(dljkira$y))
  plot(dljpdb,ylim=c(0,ym),xlim=c(0,100),
       main=paste(aa," samp ",sum(idxkira)))
  lines(dljkira,col=2)
}



# sasa5-sasa14 vs sasafrac for all residues for decoys vs. pdb
par(mfrow=c(5,4))
par(mar=c(3,2,2,1))
sasapack <- pdbdata$sasa5 - pdbdata$sasa14
sasapackkira <- kiradata$sasa5 - kiradata$sasa14
for(aa in RES) {
  #aa <- "ALA"
  #par(mfrow=c(1,1))
  idx <- pdbdata$aa==aa
  idxkira <- kiradata$aa==aa
  plot(pdbdata$sasafrac[idx],sasapack[idx],pch='.',
       xlim=c(0,1),ylim=c(0,60),
       main=paste(aa," sasapack ",sum(idxkira),sep=''))

  points(kiradata$sasafrac[idxkira],sasapackkira[idxkira],pch='.',col=2 )
}


# TLJ vs sasafrac for all residues for decoys vs. pdb
par(mfrow=c(5,4))
par(mar=c(3,2,2,1))
for(aa in RES) {
  #aa <- "ALA"
  #par(mfrow=c(1,1))
  idx <- pdbdata$aa==aa
  idxkira <- kiradata$aa==aa
  plot(pdbdata$sasafrac[idx],pdbdata$tlj[idx],pch='.',
       xlim=c(0,1),
       main=paste(aa," TLJ ",sum(idxkira),sep=''))
  points(kiradata$sasafrac[idxkira],kiradata$tlj[idxkira],pch='.',col=2 )
}


# density of TLJ binned on sasafrac for a given residue
par(mfrow=c(5,4))
par(mar=c(3,2,2,1))
aa <- "ALA"
i1p <- pdbdata$aa==aa
i1k <- kiradata$aa==aa
for(bin in 1:20) {
  i2p <-  bin*0.05-0.05 <=  pdbdata$sasafrac[i1p] &  pdbdata$sasafrac[i1p] < bin*0.05
  i2k <-  bin*0.05-0.05 <= kiradata$sasafrac[i1k] & kiradata$sasafrac[i1k] < bin*0.05
  dp <- density( pdbdata$tlj[i1p][i2p])
  dk <- density(kiradata$tlj[i1k][i2k])
  ym <- max(max(dp$y),max(dk$y))
  plot(dp,ylim=c(0,ym),
       main=paste(aa,' ',sum(i2k)," bin ",bin*0.05-0.05,'-',bin*0.05,sep=''))
  lines(dk,col=2)
}


# density of TLJ binned on sasafrac for a given bin
par(mfrow=c(5,4))
par(mar=c(3,2,2,1))
bin <- 1
i1p <-  bin*0.05-0.05 <=  pdbdata$sasafrac &  pdbdata$sasafrac < bin*0.05
i1k <-  bin*0.05-0.05 <= kiradata$sasafrac & kiradata$sasafrac < bin*0.05
for(aa in RES) {
  i2p <- pdbdata$aa[i1p]==aa
  i2k <- kiradata$aa[i1k]==aa
  dp <- density( pdbdata$tlj[i1p][i2p])
  dk <- density(kiradata$tlj[i1k][i2k])
  ym <- max(max(dp$y),max(dk$y))
  plot(dp,ylim=c(0,ym),
       main=paste(aa,' ',sum(i2k)," bin ",bin*0.05-0.05,'-',bin*0.05,sep=''))
  lines(dk,col=2)
}


# density of sasapack binned on sasafrac for a given residue
par(mfrow=c(5,4))
par(mar=c(3,2,2,1))
aa <- "ALA"
i1p <- pdbdata$aa==aa
i1k <- kiradata$aa==aa
for(bin in 1:20) {
  i2p <-  bin*0.05-0.05 <=  pdbdata$sasafrac[i1p] &  pdbdata$sasafrac[i1p] < bin*0.05
  i2k <-  bin*0.05-0.05 <= kiradata$sasafrac[i1k] & kiradata$sasafrac[i1k] < bin*0.05
  dp <- density( pdbdata$sasa5[i1p][i2p] -  pdbdata$sasa14[i1p][i2p])
  dk <- density(kiradata$sasa5[i1k][i2k] - kiradata$sasa14[i1k][i2k])
  ym <- max(max(dp$y),max(dk$y))
  plot(dp,ylim=c(0,ym),
       main=paste(aa,' ',sum(i2k)," bin ",bin*0.05-0.05,'-',bin*0.05,sep=''))
  lines(dk,col=2)
}


# density of sasapack binned on sasafrac for a given bin
par(mfrow=c(5,4))
par(mar=c(3,2,2,1))
bin <- 1
i1p <-  bin*0.05-0.05 <=  pdbdata$sasafrac &  pdbdata$sasafrac < bin*0.05
i1k <-  bin*0.05-0.05 <= kiradata$sasafrac & kiradata$sasafrac < bin*0.05
for(aa in RES) {
  i2p <- pdbdata$aa[i1p]==aa
  i2k <- kiradata$aa[i1k]==aa
  dp <- density( pdbdata$sasa5[i1p][i2p] -  pdbdata$sasa14[i1p][i2p])
  dk <- density(kiradata$sasa5[i1k][i2k] - kiradata$sasa14[i1k][i2k])
  ym <- max(max(dp$y),max(dk$y))
  plot(dp,ylim=c(0,ym),
       main=paste(aa,' ',sum(i2k)," bin ",bin*0.05-0.05,'-',bin*0.05,sep=''))
  lines(dk,col=2)
}


# compute stats for LJ sasa05 joint distribution
msasa <- array(0,c(20,11))
matr  <- array(0,c(20,11))
ssasa <- array(0,c(20,11))
satr  <- array(0,c(20,11))
rho   <- array(0,c(20,11))
for(r in 1:20){
  idx1 <- pdbdata$aa==RES[r]
  idx2 <- pdbdata$sasafrac[idx1] == 0
  msasa[r,1] <- mean( pdbdata$sasa5[idx1][idx2] )
  ssasa[r,1] <- sd(   pdbdata$sasa5[idx1][idx2] )
  matr[r,1]  <- mean( pdbdata$atr[idx1][idx2] )
  satr[r,1]  <- sd(   pdbdata$atr[idx1][idx2] )
  rho[r,1]   <- cor(  pdbdata$sasa5[idx1][idx2], pdbdata$atr[idx1][idx2] )
  for(b in 1:10){
    print(c(r,b))
    idx2 <- pdbdata$sasafrac[idx1] > 0.1*b-0.1 & pdbdata$sasafrac[idx1] <= 0.1*b
    msasa[r,b+1] <- mean( pdbdata$sasa5[idx1][idx2] )
    ssasa[r,b+1] <- sd(   pdbdata$sasa5[idx1][idx2] )
    matr[r,b+1]  <- mean( pdbdata$atr[idx1][idx2] )
    satr[r,b+1]  <- sd(   pdbdata$atr[idx1][idx2] )
    rho[r,b+1]   <- cor(  pdbdata$sasa5[idx1][idx2], pdbdata$atr[idx1][idx2] )
  }
}


par(mfrow=c(5,4))
par(mar=c(3,2,2,1))
for(r in 1:20){
  plot( matr[r,] , main=paste(RES[r],"ATR mean") )
}
dev.print(device=postscript,horizontal=F,file="fig/pdb_packing_stats_atr_mean.ps")

par(mfrow=c(5,4))
par(mar=c(3,2,2,1))
for(r in 1:20){
  plot( satr[r,] , main=paste(RES[r],"ATR sd") )
}
dev.print(device=postscript,horizontal=F,file="fig/pdb_packing_stats_atr_sd.ps")

par(mfrow=c(5,4))
par(mar=c(3,2,2,1))
for(r in 1:20){
  plot( msasa[r,] , main=paste(RES[r],"sasa5 mean") )
}
dev.print(device=postscript,horizontal=F,file="fig/pdb_packing_stats_sasa5_mean.ps")

par(mfrow=c(5,4))
par(mar=c(3,2,2,1))
for(r in 1:20){
  plot( ssasa[r,] , main=paste(RES[r],"sasa5 sd") )
}
dev.print(device=postscript,horizontal=F,file="fig/pdb_packing_stats_sasa5_sd.ps")

par(mfrow=c(5,4))
par(mar=c(3,2,2,1))
for(r in 1:20){
  plot( rho[r,] , main=paste(RES[r],"ATR sasa5 cor") )
}
dev.print(device=postscript,horizontal=F,file="fig/pdb_packing_stats_cor_atr_sasa5.ps")

binlb <- array(0,c(10))
binub <- array(0,c(10))
binlb[ 1] <- -0.1; binub[ 1] <- 0.00
binlb[ 2] <- 0.00; binub[ 2] <- 0.02
binlb[ 3] <- 0.02; binub[ 3] <- 0.05
binlb[ 4] <- 0.05; binub[ 4] <- 0.09
binlb[ 5] <- 0.09; binub[ 5] <- 0.15
binlb[ 6] <- 0.15; binub[ 6] <- 0.22
binlb[ 7] <- 0.22; binub[ 7] <- 0.30
binlb[ 8] <- 0.30; binub[ 8] <- 0.45
binlb[ 9] <- 0.45; binub[ 9] <- 0.65
binlb[10] <- 0.65; binub[10] <- 1.00

sasaprob <- array(0,c(20,10,75))
sasacount <- array(0,c(20,10))
# compute counts for sasapack
for(r in 1:20) {
  idx1 <- pdbdata$aa==RES[r]
  for(b in 1:10) {
    print(c(r,b))
    idx2 <-  binlb[b] < pdbdata$sasafrac[idx1] & pdbdata$sasafrac[idx1] <= binub[b]
    tmp <- pdbdata$sasa5[idx1][idx2]-pdbdata$sasa14[idx1][idx2]
    sasacount[r,b] <- length(tmp)
    for(ii in 1:75){
      sasaprob[r,b,ii] <- sum(tmp>ii)/sasacount[r,b]
    }
  }
}
save(sasaprob,file='rdata/sasaprob.rdata')
save(sasacount,file='rdata/sasacount.rdata')

par(mfrow=c(5,4))
par(mar=c(3,2,2,1))
colors <- rainbow(10)
for(r in 1:20) {
  plot(sasaprob[r,1,],main=paste(RES[r],"sasaprob"),type='l', col=colors[1])
  for(bin in 2:10){
    lines( sasaprob[r,bin,] , col=colors[bin] )
  }
}
dev.copy2eps(file="fig/sasaprob_probabilities.eps")

par(mfrow=c(5,4))
par(mar=c(3,2,2,1))
for(r in 1:20){
  m <- max(sasacount[r,])
  plot(sasacount[r,], ylim=c(0,m),
       main=paste(RES[r],"mx",m,"mn",min(sasacount[r,])))
}
dev.copy2eps(file="fig/sasaprob_counts.eps")

# write sasaprob data file
write(20,file="sasa_prob_cdf.txt")
for(r in 1:20) {
  write(paste(r,RES[r],75,10),file="sasa_prob_cdf.txt",append=T)
  write(sasaprob[r,,],file="sasa_prob_cdf.txt",ncolumns=10,append=T)
}


# sasabin  <- ceiling(pdbdata$sasafrac*10)+1
# sasadiff <- pdbdata$sasa5 - pdbdata$sasa14
calc.sasaprob <- function(ii) {
  sasabin  <- ceiling(pdbdata$sasafrac[ii]*10)+1
  aa <- which(RES==pdbdata$aa[ii])
  sasadiff <- floor(pdbdata$sasa5[ii] - pdbdata$sasa14[ii])
  sasadiff <- max(1,min(75,sasadiff),sasadiff)
  return(sasaprob[aa,sasabin,sasadiff])
}

tmp <- 1:207
for(ii in 1:207) {
  tmp[ii] <- calc.sasaprob(ii)
}
