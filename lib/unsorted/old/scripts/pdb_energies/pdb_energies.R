#pdbe <- read.table('pdb_energies.data',header=F)[,-1]
#enames <- c('atre','repe','sole','probe','dune',
#            'intrae','hbe','paire','rese','tlje')
#names(pdbe) <- c('prot','pdb','res','aa','ss','nb',
#                 'sasafrac','sasa14','sasa5',
#                 enames,'sasapack','sasaprob'
#                 )
#sasadiff <- pdbe$sasa5 - pdbe$sasa14
#pdbe <- cbind( pdbe , sasadiff )
#save(pdbe,file='pdbe.rdata')

load('pdbe.rdata')

enames <- c('atre','repe','sole','probe','dune',
            'intrae','hbe','paire','rese','tlje','sasadiff')

resnames <- c("ALA","CYS","ASP","GLU","PHE","GLY","HIS","ILE","LYS","LEU",
              "MET","ASN","PRO","GLN","ARG","SER","THR","VAL","TRP","TYR")


#pdbe$sasafrac[pdbe$sasafrac > 1] <- 1

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
  bin <- c()
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

#sfbin <- calcbinsquantile(pdbe$sasafrac,10)
#nbbin <- calcbinsquantile(pdbe$nb+runif(length(pdbe$nb))-0.5,10)
#print(levels(sfbin))
#print(levels(nbbin))

#sf <- pdbe$sasafrac
#nb <- pdbe$nb
#fuzzysf <- sf - runif(length(sf))/1000
#fuzzysf[sf==0] <- 0

#sfqbreaks <- quantilebreaksbyfactor(fuzzysf,pdbe$aa,10)
#sfbin     <- calcbinsbyfactor(      fuzzysf,pdbe$aa,sfqbreaks)
#fuzzynb   <- nb+runif(length(nb))
#nbqbreaks <- quantilebreaksbyfactor(fuzzynb,pdbe$aa,10)
#nbbin     <- calcbinsbyfactor(      fuzzynb,pdbe$aa,nbqbreaks)

calcestats <- function(data,enames,bins) {
  estats <- array(0,c(length(enames),
                      length(levels(data$aa)),
                      length(levels(data$ss)),
                      length(bins),
                      16,
                      103))
  dimnames(estats) <- list(enames,
                           levels(data$aa),
                           levels(data$ss),
                           names(bins),
                           paste('bin',padzeros(1:16)),
                           c('count','sdev',paste('q',0:100/100)))
  names(dimnames(estats)) <- c('etype','res','ss','bins','bin','state')
  quant <- function(x) quantile(x,0:100/100)
  for(e in enames) {
    for(b in names(bins))  {
      print(paste(e,b))
      estats[e,,,b,,'count'] <- list2array(tapply(data[,e],list(data$aa,data$ss,bins[[b]]),length))
      estats[e,,,b,,'sdev' ] <- list2array(tapply(data[,e],list(data$aa,data$ss,bins[[b]]),sd))
      estats[e,,,b,,3:103  ] <- list2array(tapply(data[,e],list(data$aa,data$ss,bins[[b]]),quant))
    }
  }
  estats[,,,,,'count'] [is.na(estats[,,,,,'count'])] <- 0

  return(estats)
}

nsfbin <- c(16,8,4,2,1)
nnbbin <- c(1,2,4,8,16)
pdbebins <- as.list(1:length(nsfbin))
names(pdbebins) <- c('nb01sf16','nb02sf08','nb04sf04','nb08sf02','nb16sf01')
for(ii in 1:length(nsfbin)) {
  sfbreaks <- quantilebreaksbyfactor(fuzzysf,pdbe$aa,nsfbin[ii])
  nbbreaks <- quantilebreaksbyfactor(fuzzynb,pdbe$aa,nnbbin[ii])
  sfbin    <- calcbinsbyfactor(      fuzzysf,pdbe$aa,sfbreaks)
  nbbin    <- calcbinsbyfactor(      fuzzynb,pdbe$aa,nbbreaks)
  pdbebins[[ii]] <- as.factor(paste('nb',nbbin,'sf',sfbin,sep=''))
}
enames <- enames
data <- pdbe[1:100,]
bins <- lapply(pdbebins,function(x) x[1:100])

#tmp <- calcestats(data,enames,bins)

pdbestats <- calcestats(pdbe,enames,pdbebins)
save(pdbestats,file='pdbestats.rdata')


# todo look at histograms of counts in bins
# figure out why the counts suck so much
# compute bins for each SS/AA combo!!
apply(pdbestats[,,,,,'count'],4,density)
plot(c(),c(),type='n',xlim=
for(ii in 1:5)


  
par(mfrow=c(4,5))
for(ii in 1:length(RES)) {  
  xmin
  plot(
}
