#################################################################3
#        with SS and more decoys
#################################################################


##### start analysis

# decided not to use these because I don't quite get what ROCR is doing
#library(gtools,lib.loc='~/software/Rlib')
#library(gdata, lib.loc='~/software/Rlib')
#library(gplot,lib.loc='~/software/Rlib')
#library(ROCR,  lib.loc='~/software/Rlib')



calc.joint.bins <- function(bp,fieldlist,namesubs,fsep=".",nsep="") {
  if(missing(namesubs))
    sub <- list()
  bins <- list()
  bins[['nobins']] <- as.factor(rep('nobins',dim(bp)[1]))
  for(ii in 1:length(fieldlist)) {
    fields <- fieldlist[[ii]]
    name <- paste(fieldlist[[ii]],collapse=nsep)
    for(sub in namesubs)
      name <- sub(sub[1],sub[2],name)
    print(paste(c("calculating joint factor on fields:",fields),collapse=" "))
    jointfield <- as.character(bp[[fields[1]]])
    if(length(fields)>1) {
      for(jj in 2:length(fields)) {
        field <- fields[jj]
        jointfield <- paste(jointfield,as.character(bp[[field]]),sep=fsep)
      }
    }
    bins[[name]] <- as.factor(jointfield)
  }
  return(bins)  
}

calc.frac.unsatisfied <- function(unsatisfied,bins) {
  uf <- list()
  for(bin in names(bins)) {
    ct <- tapply(unsatisfied, bins[[bin]], length)
    uc <- tapply(unsatisfied, bins[[bin]], sum)
    uf[[bin]] <- uc/ct
  }
  return(uf)
}

calc.pred <- function(bins,predictor,satisfied=F,rm.satisfied=F) {
  pred <- list()
  for(bin in names(bins)) {
    tmp <- bins[[bin]]
    levels(tmp) <- predictor[[bin]]
    tmp <- as.numeric(as.character(tmp))
    if(rm.satisfied)
      tmp <- tmp[!satisfied]
    else if(any(satisfied))
      tmp[satisfied] <- 1
    pred[[bin]] <- tmp
  }
  return(pred)
}


# function to plot roc curves
# labels must be true and false
calc.roc <- function(posval,negval,n=10,highgood=T) {
  pred  <- c(posval,negval)
  label <- c(rep(T,length(posval)),rep(F,length(negval)))
  tpr <- 0:(n-1)/(n-1)
  fpr <- 0:(n-1)/(n-1)
  q <- quantile(pred,0:(n-1)/(n-1))  
  for(ii in 2:(n-1)) {
    x <- sum(pred<=q[ii] & !label)/sum( !label )
    y <- sum(pred<=q[ii] &  label)/sum(  label )
    fpr[ii] <- x
    tpr[ii] <- y
  }
  if(highgood) {
    tpr <- rev(1-tpr)
    fpr <- rev(1-fpr)
  }
  score <- 0
  for(ii in 2:n) {
    #print(paste(ii,score, mean(c(tpr[ii-1],tpr[ii])),abs(fpr[ii]-fpr[ii-1])))
    score <- score + mean(c(tpr[ii-1],tpr[ii]))*abs(fpr[ii]-fpr[ii-1])
  }
  return(list(y=tpr,x=fpr,score=score))
}

calc.roc.categories <- function(posvals,negvals,n=10,highgood=T) {
  if(any(names(posvals)!=names(negvals)))
    print("ERROR calc.roc.categories names for pos and neg don't match")
  roc <- list()  
  for(bin in names(posvals)) {
    roc[[bin]] <- calc.roc(posvals[[bin]],negvals[[bin]],n,highgood)
  }
  return(roc)
}



plot.roc <- function(roc,title="Some ROC Curves",colors=NULL,file=NULL,legx=0.50,legy=0.35,xlim=c(0,1),ylim=c(0,1)) {
  rocscores <- format(round(as.numeric(lapply(roc,function(x) x$score)),3))
  if(is.null(colors))
    colors = rainbow(length(roc))
  plot (c(0,1),c(0,1) , col=1 , type='n',main=title,,xlim=xlim,ylim=ylim,
        xlab="Fraction False Positives",ylab="Fraction True Positives")
  for(ii in 1:length(roc)) {
    lines(roc[[ii]], col=colors[ii],lwd=2)
    #lines(roc[[ii]]$x,roc[[ii]]$y**(1/3), col=colors[ii],lwd=2)
  }
  legend(legx,legy,c(paste(rocscores,'on:',names(roc))), col=c(colors),
         lwd=c(rep(3,length(roc))))
  if(!is.null(file))
    dev.print(file=file,device=postscript,horizontal=F)       
}



mydensity <- function(x) {
  if(length(x)==1) return(list(x=c(x,x),y=c(0,1)))
  else             return(density(x))
}

gt  <- function(a,b) return( a >  b )
gte <- function(a,b) return( a >= b )
lt  <- function(a,b) return( a <  b )
lte <- function(a,b) return( a <= b )

padzeros <- function(x,n=2) {
  x <- as.character(x)
  y <- x
  for(ii in 1:length(x))
    y[ii] <- paste( c(rep("0",(n-nchar(x[ii]))) , x[ii]) , sep="" , collapse='' )
  return(y)
}

calcbins <- function(data,breaks){
  bin <- rep(NA,length(data))
  for(ii in 1:(length(breaks)-1)) {
    lcmp <- lte
    rcmp <- lt
    if( ii == 1) {
      lcmp <- lte
    } else if( breaks[ii] == breaks[ii+1] ) {
      lcmp <- lte
    } else if( breaks[ii-1] == breaks[ii] ) {
      lcmp <- lt
    }
    if( ii == length(breaks)-1 ) {
      rcmp <- lte
    } else if(breaks[ii]==breaks[ii+1]) {
      rcmp <- lte
    } else if(breaks[ii+1]==breaks[ii+2]) {
      rcmp <- lt
    }
    #lbl <- paste(format(breaks[ii],  digits=3),'-',
    #             format(breaks[ii+1],digits=3),sep='')
    lbl <- padzeros(ii,2)
    #lbl <- padzeros(ii,nchar(as.character(length(breaks)-1)))
    bin[ lcmp(breaks[ii],data) & rcmp(data,breaks[ii+1]) ] <- lbl
  }
  bin <- as.factor(bin)
  return(bin)
}
#tmp <- calcbins(pdbe$sasafrac,c(0,0,0.5,1,1))

calcbinsbyfactor <- function(data,fac,breaklist) {
  bin <- c()
  for(l in levels(fac)){
    b <- calcbins(data[fac==l],breaklist[[l]])
    bin[fac==l] <- as.character(b)
  }
  bin <- as.factor(bin)
  return(bin)
}

# TODO: write joint bins function
calcjointbins <- function(data,breaks) {
  
}

#calcsasabins <- function(sasafrac) {
#  return(calcbins(sasafrac,c(0,0,0.02,0.05,0.09,0.15,0.22,0.30,0.45,0.60,1.0)))
#}


quantilebreaks <- function(data,nbins) {
  # TODO!!! make this take big modes as special case rather than just 0!
  N <- length(data)
  if(sum(data==0)/N > 1/nbins) {
    return(c(0,0,quantile(data[data!=0],1:(nbins-1)/(nbins-1))))
  } else {
    return(quantile(data,0:nbins/nbins))
  }
}

quantilebreaksbyfactor <- function(data,fac,nbins) {
  breaks <- as.list(1:length(levels(fac)))
  for(ii in 1:length(levels(fac))){
    breaks[[ii]] <- quantilebreaks(data[fac==levels(fac)[ii]],nbins)
    breaks[[ii]][1] <- min(data)
    breaks[[ii]][length(breaks[[ii]])] <- max(data)
  }
  names(breaks) <- levels(fac)
  return(breaks)
}

#jointquantilebreaks <- function(data,nbins) {
#  n <- dim(data)[,2]
#  breaks <- as.list()
#  breaks1 <- quantilebreaks(data,nbins1)
#  breaks2 <- quantilebreaks(data,nbins2)
#  return(list(breaks1,breaks2))
#}

# TODO: write joint breaks function
# should specify desired # per bin and some kind of weight to prefer more
# bins for one "factor" vs. another
# do this based on quantiles of whitened data!

# X <- cbind(pdbe$nb+runif(dim(pdbe)[1])-0.5,pdbe$sasafrac)
# pc <- prcomp(X)
# C <- t(t(X)-pc$center)
# S <- t(t(C)/pc$sdev)
# R <- S %*% pc$rotation

# nb <- pdbe$nb + runif(length(pdbe$nb))
# sf <- pdbe$sasafrac
# nb <- ( nb - mean(nb) ) / sd( nb )
# sf <- ( sf - mean(sf) ) / sd( sf )

# plot( nb , sf , pch='.')

list2array <- function(l) {
  if(is.numeric(l))
    return(l)
  ii <- 1
  while(is.null(l[[ii]]))
    ii <- ii+1        
  d <- c(dim(l),length(l[[ii]]))
  n <- append(dimnames(l),list(names(l[[ii]])))
  a <- array(NA,d)
  dimnames(a) <- n
  if(length(dim(l))==2) {
    for(ii in 1:d[1])
      for(jj in 1:d[2])
        if(!is.null(l[[ii,jj]]))
          a[ii,jj,] <- l[[ii,jj]]
    return(a)
  } else if( length(dim(l))==3 ) {
    for(ii in 1:d[1])
      for(jj in 1:d[2])
        for(kk in 1:d[3])
          if(!is.null(l[[ii,jj,kk]]))
            a[ii,jj,kk,] <- l[[ii,jj,kk]]
    return(a)
  } else {
    print(paste('list dim',length(dim(l)),'not handled yet'))
  }
}
