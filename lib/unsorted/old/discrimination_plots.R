load('~/project/features/discrimination/rdata/ds.ds3.rtuf.decoyscores.rdata')
source('~/scripts/buried_polar/buried_polar_util.R')

setprots <- levels(ds$setprot)
rocss <- array(0,c(length(setprots),4))
colnames(rocss) <- c('ss','rocscore','rocsasa','rocbpe')

#test <- function(x) {

#par(mfrow=c(4,5))
par(mfrow=c(1,1))
par(mgp=c(2,1,0))
par(mar=c(3,3,4,1))
par(cex=1.0)

#for(ii in (4*x-3):(4*x)){#length(setprots)) {
for(ii in 1:length(setprots)) {

  cat  <-unique(as.character(ds$category[ds$setprot==setprots[ii]]))
  if(any(ds$setprot==setprots[ii]) && 'gsbr'%in%cat&& 'gsbr'%in%cat&& 'relax_native'%in%cat){

    #print(setprots[ii])
    
    set  <- strsplit(setprots[ii],' ')[[1]][1]
    prot <- strsplit(setprots[ii],' ')[[1]][2]
    print(c(set,prot))
    
    all  <- ds [ ds$setprot==setprots[ii] , ]
    gsbr <- all[ all$category=='gsbr' ,]
    gsgr <- all[ all$category=='gsgr' ,]
    nat  <- all[ all$category=='relax_native' ,]
    H <- gsbr$H[1]
    E <- gsbr$E[1]
    L <- gsbr$L[1]
    
    xlim <- c(0,7)#max(gsbr$rms))
    
    ylim <- c(min(all$score),min(max(all$score),0))
    plot(0:1,0:1,xlim=xlim,,ylim=ylim,type='n', xlab="RMSD", ylab="SCORE",
         main=paste(setprots[ii],"score\n",'H',round(H,2),'E',round(E,2),'L',round(L,2)))
    points(gsbr$rms,gsbr$score,col=2,pch=16)
    points(gsgr$rms,gsgr$score,col=3,pch=16)
    points(nat$rms,nat$score,col=1,pch=16)
    dev.copy2eps(file=paste("~/project/research_reports/fig/decoy_examples/",set,"_",prot,"_score.eps",sep=""))
    
    #ylim <- c(min(all$sasaprob),max(all$sasaprob))
    ##ylim <- c(0,1)
    #plot(0:1,0:1,xlim=xlim,,ylim=ylim,type='n',
    #     main=paste(setprots[ii],"sasaprob score\n",'H',round(H,2),'E',round(E,2),'L',round(L,2)))
    #points(gsbr$rms,gsbr$sasaprob,col=2,pch=16)
    #points(gsgr$rms,gsgr$sasaprob,col=3,pch=16)
    #points(nat$rms, nat$sasaprob, col=1,pch=16)
    
    #ylim <- c(min(all$bpe),max(all$bpe))
    ##ylim <- c(0,1)
    #plot(0:1,0:1,xlim=xlim,,ylim=ylim,type='n',
    #     main=paste(setprots[ii],"buried polar score\n",'H',round(H,2),'E',round(E,2),'L',round(L,2)))
    #points(gsbr$rms,gsbr$bpe,col=2)
    #points(gsgr$rms,gsgr$bpe,col=3)
    #points(nat$rms, nat$bpe, col=1)
    
    ylim <- c(min(-all$sasaprob-3*all$bpe),max(-all$sasaprob-3*all$bpe))
    #ylim <- c(0,1)
    plot(0:1,0:1,xlim=xlim,,ylim=ylim,type='n', xlab="RMSD", ylab="BUP_SPK SCORE",
         main=paste(setprots[ii],"buried polar score\n",'H',round(H,2),'E',round(E,2),'L',round(L,2)))
    points(gsbr$rms,-(gsbr$sasaprob+3*gsbr$bpe),col=2,pch=16)
    points(gsgr$rms,-(gsgr$sasaprob+3*gsgr$bpe),col=3,pch=16)
    points(nat$rms, -(nat$sasaprob+3*nat$bpe), col=1,pch=16)
    dev.copy2eps(file=paste("~/project/research_reports/fig/decoy_examples/",set,"_",prot,"_spkbup.eps",sep=""))
    
#    ylim <- c(min(all$bpenaive),max(all$bpenaive))
#    #ylim <- c(0,1)
#    plot(0:1,0:1,xlim=xlim,,ylim=ylim,type='n',
#         main=paste(setprots[ii],"bpenaive score\n",'H',round(H,2),'E',round(E,2),'L',round(L,2)))
#    points(gsbr$rms,gsbr$bpenaive,col=2)
#    points(gsgr$rms,gsgr$bpenaive,col=3)
#    points(nat$rms, nat$bpenaive, col=1)

    #posvals <- list(score=-gsgr$score,sasaprob=gsgr$sasaprob,bpe=gsgr$bpe,bpenaive=gsgr$bpenaive)
    #negvals <- list(score=-gsbr$score,sasaprob=gsbr$sasaprob,bpe=gsbr$bpe,bpenaive=gsbr$bpenaive)
    #roc <- calc.roc.categories(posvals,negvals)
    #plot.roc(roc,colors=c('black','blue','red','pink'))

    #rocss[ii,'ss' ] <- H-E
    #rocss[ii,'rocscore'] <- roc$score$score
    #rocss[ii,'rocsasa']  <- roc$sasaprob$score
    #rocss[ii,'rocbpe']   <- roc$bpe$score
    
  } else {
    plot(0:1,0:1,type='n',main=paste(setprots[ii],"score"))
    plot(0:1,0:1,type='n',main=paste(setprots[ii],"sasaprob score"))
    plot(0:1,0:1,type='n',main=paste(setprots[ii],"buried polar score"))
    plot(0:1,0:1,type='n',main=paste(setprots[ii],"naive buried polar score"))
  }
  #if(ii%%4==0) {
  #  dev.print(device=postscript,horizontal=T,
  #            file=paste('fig/score_v_rms_',ii/4,'.ps',sep=''))
  #}
}
#}

if(0) {
plot(rocss[,1],rocss[,2],main="SS (high=helix) vs. Score Roc")
lines(c(-10,10),c(0.5,0.5),col=2)
plot(rocss[,1],rocss[,3],main="SS (high=helix) vs. Sasaprob Roc")
lines(c(-10,10),c(0.5,0.5),col=2)
plot(rocss[,1],rocss[,4],main="SS (high=helix) vs. Buried Polar Roc")
lines(c(-10,10),c(0.5,0.5),col=2)

dev.print(device=postscript,horizontal=T,
          file=paste('fig/score_v_rms_',ceiling(ii/4),'.ps',sep=''))

}


col <- 1 + as.numeric(ds$setprot)%%7
plot(ds$rms,ds$score,pch='.',xlim=c(0,7),ylim=c(-3,0),col=col)
     
