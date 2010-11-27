dr <- list()
Nprots <- 0
for(set in c("bqian_caspbench","kira_hom",'sraman_nmr','phil_homolog','kira_bench')) {
  print(paste('loading dr.',set,'.ds3.region.rdata',sep=''))
  load(paste('dr.',set,'.ds3.region.rdata',sep=''))
  regionscores <- regionscores[ regionscores$cat %in% c('gsgr','gsbr') ,]
  p <- 0.7;
  regionscores$both7  <- p*(1-regionscores$bup7 ) + (1-p)*regionscores$spk7  
  regionscores$both10 <- p*(1-regionscores$bup10) + (1-p)*regionscores$spk10 
  Nprots = Nprots + length(unique(regionscores$prot))
  dr[[set]] <- regionscores
  rm(regionscores)
}
#load(paste('dr.',set,'.ds3.region.rdata',sep=''))

screendisplay <- T
screendisplay <- F

set <- 'bqian_caspbench'
for(set in c("bqian_caspbench","kira_hom",'sraman_nmr','phil_homolog','kira_bench')) {

  for(prot in unique(dr[[set]]$prot)) {
    
    nprot <- length(unique(dr[[set]]$prot))
    if(screendisplay){
      x11()
    } else {
      png(filename=paste('fig/region_',set,'_',prot,'.png',sep=''),
          height=6*200,width=5*180)
    }
    I <- dr[[set]]$prot==prot

    par(mfcol=c(6,5))
    par(mar=c(1,1,2,1))
    par(mgp=c(1,0,0))
    par(tck=0)
    par(cex=0.75)
    colors=rainbow(30)[1:25]
    colors= c(colors,rep(colors[25],20))
    defaultxlim <- c(0,5)
    defaultylim=NULL
    print(paste('plotting ','fig/region_',set,'_',prot,sep=''))
    
    for(score in c('atr','res','spk','bup','both')) {
      score7  <- dr[[set]][[paste(score, '7',sep='')]]
      score10 <- dr[[set]][[paste(score,'10',sep='')]]
      for(rms in c('rmsabs','rmsrel','rms50')) {
        xlim <- defaultxlim; if(rms=="rmsabs") xlim <- c(0,15)
        for(reg in c("7","10")){
          rmsreg <- rms; if(rms %in% c('rmsabs','rmsrel')) rmsreg  <- paste(rms,reg,sep='')
          plot(dr[[set]][[rmsreg]][I],  get(paste('score',reg,sep=''))[I] ,
               pch='.', col=colors[dr[[set]][[paste('nb',reg,sep='')]][I]],
               xlim=xlim,ylim=defaultylim,
               main=paste(prot,reg,rmsreg,score))
        }
      }
    }
    if(!screendisplay)
      print(paste('writing file ','fig/region_',set,'_',prot,'.png',sep=''))
    dev.off()
  }
}

dev.print(file='fig/region_bqian_caspbench_spk_bup.ps',device=postscript,horizontal=T)
dev.print(file='fig/region_bqian_caspbench_spk_bup.ps',device=postscript,horizontal=T)
dev.print(file='fig/region_bqian_caspbench_spk_bup.ps',device=postscript,horizontal=T)

    
finddim <- function(n,height,width) {
  h <- ceiling(n/sqrt(n)*height/width)
  v <- ceiling(n/h)
  return(c(h,v))
}

x11()
png(filename="fig/region_all_prots_bup_spk_vs_rms50_rms100_3positions.png",width=900,height=1200)

par(mfrow=finddim(4*61,1,1))#c(6,3))
par(mar=c(0.2,0.2,1.0,0.2))
par(mgp=c(2000,2000,2000))
par(tck=0)
par(cex=0.5)
par(oma=c(0,0,4,0))

colors7 = rainbow(16)[1:13]
colors7 = c(colors7,rep(colors7[13],20))
colors10= rainbow(30)[1:25]
colors10= c(colors10,rep(colors10[25],20))
defaultxlim <- c(0,5)
defaultylim <- c(0,0.5)


print(paste('plotting ','fig/region_',set,'_',prot,sep=''))

set <- 'bqian_caspbench'
for(set in c("bqian_caspbench","kira_hom",'sraman_nmr','phil_homolog','kira_bench')) {
  for(prot in unique(dr[[set]]$prot)) {
    print(paste("plotting",set,prot))
    
    #I <- dr[[set]]$prot==prot
    #plot( dr[[set]]$rms100[I], dr[[set]]$both7[I], pch='.',col='black',
    #     xlim=defaultxlim,ylim=defaultylim,main=prot )
    #points(dr[[set]]$rms50[I], dr[[set]]$both7[I], pch='.', col='grey')
    
    I <- dr[[set]]$prot==prot
    Nres <- max(dr[[set]]$res[I])
    color <- rainbow(1.2*Nres)[1:Nres]
    color <- c('red','blue','green')
    names(color) <- c('H',"L","S")
    I <- I & dr[[set]]$res %in% c(10,floor(Nres/2),Nres-10)
    color[10] = 'red'
    color[floor(Nres/2)] <- 'black'
    color[Nres-10] <- 'blue'
    
    tmpc <- c(rep('black',sum(I)))#,rep('red',sum(I)))
    tmps <- c(dr[[set]]$ss[I])#,dr[[set]]$res[I])
    tmpr <- c(dr[[set]]$res[I])#,dr[[set]]$res[I])
    tmpx <- c(dr[[set]]$rms100[I])#,dr[[set]]$rms50[I])
    tmpy <- c(dr[[set]]$both10[I])#,dr[[set]]$both10[I])
    o <- order(runif(sum(I)))#*2))

    
    plot( tmpx[o] , tmpy[o] ,  pch='+', col=color[tmpr[o]],#
         xlim=defaultxlim,ylim=defaultylim,xlab=NULL, ylab=NULL,
         main=paste(substr(set,1,10),prot) )

    Ir <- I &  dr[[set]]$res == 10
    plot( dr[[set]]$rms100[Ir] , dr[[set]]$both10[Ir] ,  pch='+', col='red',
         xlab=NULL, ylab=NULL,
         main=paste(substr(set,1,10),prot) )

    Ir <- I &  dr[[set]]$res == floor(Nres/2)
    plot( dr[[set]]$rms100[Ir] , dr[[set]]$both10[Ir] ,  pch='+', col='black',
         xlab=NULL, ylab=NULL,
         main=paste(substr(set,1,10),prot) )

    Ir <- I &  dr[[set]]$res == Nres-10
    plot( dr[[set]]$rms100[Ir] , dr[[set]]$both10[Ir] ,  pch='+', col='blue',#
         xlab=NULL, ylab=NULL,
         main=paste(substr(set,1,10),prot) )

    #I <- dr[[set]]$prot==prot & dr[[set]]$cat=='gsgr'
    #plot( dr[[set]]$rms50[I], dr[[set]]$both10[I], pch='.',col='green',
    #     xlim=defaultxlim,ylim=defaultylim,main=prot )
    #I <- dr[[set]]$prot==prot & dr[[set]]$cat=='gsbr'
    #points(dr[[set]]$rms50[I], dr[[set]]$both10[I], pch='.', col='red')
  }
  
}

mtext("0.7*BUP + 0.3*SPK score vs. RMS100 helix=red sheet=green loop=blue; Xaxis rms 0-5, Yaxis score 0-0.5",
      outer=T)

dev.off()


marksselements <- function(ss) {
  e <- c()
  g <- 1
  for(ii in 1:(length(ss)-1)) {
    e[ii] <- g
    if(ss[ii]!=ss[ii+1]) 
      g <- g+1
  }
  e[length(ss)] <- g

  count <- 0
  for(ii in 1:e[length(e)]){
    if(sum(e==ii)<5) {
      e[e==ii] <- 0      
    } else {
      count <- count+1
      e[e==ii] <- count
    }
  }
  return( e)
}

for(set in names(dr)) {
  dr[[set]]$element <- rep(0,length(dr[[set]]$ss))
  for(prot in unique(dr[[set]]$prot)) {
    print(paste(set,prot))
    Iprot <- dr[[set]]$prot==prot
    for(pdb in unique(dr[[set]]$pdb[Iprot])) {
      print(pdb)
      Ipdb <- dr[[set]]$pdb[Iprot]==pdb
      print(sum(Ipdb))
      ss <- dr[[set]]$ss[Iprot][Ipdb]
      print(ss)
      dr[[set]]$element[Iprot][Ipdb] <- marksselements(ss)
    }
  }

}

