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

setwd('~/project/features/energy_v_pdb')
#load('pdbe.rdata')
load('~/data/rdata/project_features_energy_v_pdb_pe.pdb.energies.rdata')
pdbe <- pe

enames <- c('sasa14')  #('atre','repe','sole','probe','dune',
              #'intrae','hbe','paire','rese','tlje','sasapack')

resnames <- c("ALA","CYS","ASP","GLU","PHE","GLY","HIS","ILE","LYS","LEU",
              "MET","ASN","PRO","GLN","ARG","SER","THR","VAL","TRP","TYR")
res2num  <- 1:20
names(res2num) <- resnames

source('~/scripts/R-utils.R')

#########################################3
#   oct 4 05
#####################################3

pdbe$sasafrac <- pmin(pdbe$sasafrac,1)
pdbe$nb       <- pmax(pdbe$nb,6)
r <- runif(dim(pdbe)[1])-0.5

aa    <- as.character(pdbe$aa)
ss    <- as.character(pdbe$ss)
aass  <- paste(aa,'ss',ss,sep='')
aaf   <- as.factor(aa)
ssf   <- as.factor(ss)
aassf <- as.factor(aass)

brk <- list()
brk$nb      <- c(5:30+0.5)
#brk$sf      <- quantilebreaks(pdbe$sasafrac,25)
brk$aasf    <- quantilebreaksbyfactor(pdbe$sasafrac,aaf,  25)
brk$sssf    <- quantilebreaksbyfactor(pdbe$sasafrac,ssf,  25)
brk$aasssf  <- quantilebreaksbyfactor(pdbe$sasafrac,aassf,25)

#  AAAssSsf01nb01

bin <- list()
#bin$none   <- as.factor(rep('...ss.sf..nb..',length(aass)))
#bin$aa     <- as.factor(paste(aa,'ss.sf..nb..',sep=''))
#bin$ss     <- as.factor(paste('...ss',ss,'sf..nb..',sep=''))
#bin$aass   <- as.factor(paste(aass,'sf..nb..',sep=''))
#bin$nb     <- as.factor(paste('...ss.sf..nb', as.character(calcbins(pdbe$nb,brk$nb)),sep=''))
#bin$aanb   <- as.factor(paste(aa,'ss.sf..nb', as.character(calcbins(pdbe$nb,brk$nb)),sep=''))
bin$aassnb <- as.factor(paste(aass,'sf..nb',  as.character(calcbins(pdbe$nb,brk$nb)),sep=''))
#bin$sf     <- as.factor(paste('...ss.sf',     as.character(calcbins(pdbe$sasafrac,brk$sf)),'nb..',sep=''))
#bin$aasf   <- as.factor(paste(aa,'ss.sf',as.character(calcbinsbyfactor(pdbe$sasafrac,aaf,brk$aasf)),'nb..',sep=''))
#bin$sssf   <- as.factor(paste('...ss',ss,'sf',as.character(calcbinsbyfactor(pdbe$sasafrac,ssf,brk$sssf)),'nb..',sep=''))
#bin$aasssf <- as.factor(paste(aass,'sf',as.character(calcbinsbyfactor(pdbe$sasafrac,aassf,brk$aasssf)),'nb..',sep=''))

#brk$aasf5    <- quantilebreaksbyfactor(pdbe$sasafrac,aaf,  5)
#brk$aasssf5  <- quantilebreaksbyfactor(pdbe$sasafrac,aassf,5)
#bin.aasf5    <- calcbinsbyfactor(pdbe$sasafrac,aaf,  brk$aasf5  )
#bin.aasssf5  <- calcbinsbyfactor(pdbe$sasafrac,aassf,brk$aasssf5)
#aasf         <- as.factor(paste(aa,  'sf',as.character(bin.aasf5  ),sep=''))
#aasssf       <- as.factor(paste(aass,'sf',as.character(bin.aasssf5),sep=''))
#brk$aanb5    <- quantilebreaksbyfactor(pdbe$nb+r,aasf,  5)
#brk$aassnb5  <- quantilebreaksbyfactor(pdbe$nb+r,aasssf,5)
#bin.aanb5    <- calcbinsbyfactor(pdbe$nb,aasf,  brk$aanb5  )
#bin.aassnb5  <- calcbinsbyfactor(pdbe$nb,aasssf,brk$aassnb5)
#bin$aasfnb   <- as.factor(paste(aa,'ss.sf',as.character(bin.aasf5),  'nb',as.character(bin.aanb5),  sep=''))
#bin$aasssfnb <- as.factor(paste(aass, 'sf',as.character(bin.aasssf5),'nb',as.character(bin.aassnb5),sep=''))
#rm(bin.aasf5,bin.aasssf5,bin.aanb5,bin.aassnb5)

sd <- no <- d.sd <- d.no <- list()
for(etype in enames) {
  #etype <- 'atre'
  sd[[etype]] <- no[[etype]] <- list()
  for(b in names(bin)) {
    sd[[etype]][[b]] <- tapply(pdbe[[etype]],bin[[b]],sd)
    no[[etype]][[b]] <- tapply(pdbe[[etype]],bin[[b]],length)
  }
  
  d.sd[[etype]] <- d.no[[etype]] <- list()
  for(b in names(bin)) {
    d.sd[[etype]][[b]] <- mydensity(sd[[etype]][[b]])
    d.no[[etype]][[b]] <- mydensity(no[[etype]][[b]]/mean(no[[etype]][[b]]))
  }
}

par(mfrow=c(4,3))
par(mar=c(3,3,3,1))
for(etype in enames) {
  xlim <- c(0,max(as.numeric(lapply(d.sd[[etype]],function(x) max(x$x)))))
  ylim <- c(0,max(as.numeric(lapply(d.sd[[etype]],function(x) max(x$y)))))
  plot(0,0,type='n',xlim=xlim,ylim=ylim,main=paste(etype,"std. dev."))
  col=list(none='black',aa='green',nb='blue',sf='red',sfnb='purple')
  for(b in names(bin)) {
    lines(d.sd[[etype]][[b]],col=col[[b]])
  }
}
xlim <- c(0,4)#max(as.numeric(lapply(d.no[[etype]],function(x) max(x$x)))))
ylim <- c(0,max(as.numeric(lapply(d.no[[etype]],function(x) max(x$y)))))
plot(0,0,type='n',xlim=xlim,ylim=ylim,main=paste("counts"))
col=list(none='black',aa='green',nb='blue',sf='red',sfnb='purple')
for(b in names(bin)) {
  lines(d.no[[etype]][[b]],col=col[[b]])
}
dev.print(file='fig/051905_pdb_energy_binning_std_dev_and_counts.ps',device=postscript,horizontal=F)

getaa     <- function(f) return(substr(f,1,3))
getss     <- function(f) return(substr(f,6,6))
getlowsf  <- function(f,offset=0){
  if(substr(f, 9,10)=="..")
    return("..")
  if(substr(f,13,14)=="..")
    return(brk$sf[as.numeric(substr(f, 9,10))+offset])
  return(brk$aasf5[[substr(f,1,3)]][as.numeric(substr(f, 9,10))+offset])
}
gethighsf <- function(f){
  getlowsf(f,1)
}
getlownb  <- function(f,offset=0){
  if(substr(f,13,14)=="..")
    return("..")
  if(substr(f, 9,10)=="..")
    return(brk$nb[as.numeric(substr(f,13,14))+offset])
  return(brk$aanb5[[paste(substr(f,1,3),substring(f,7,10),sep='')]][as.numeric(substr(f,13,14))+offset])
}
gethighnb <- function(f){
  getlownb(f,1)
}
  
q <- list()
for(etype in enames) {
  d   <- as.data.frame(round(list2array(tapply(pdbe[[etype]],bin$aassnb,function(x) quantile(x,0:100/100))),3))
  nsamp <- tapply(pdbe[[etype]],bin$aassnb,length)
  n   <- rownames(d)
  aa  <- getaa(n)  
  aai <- res2num[aa]
  ss  <- getss(n)
  #lowsf  <- round(as.numeric(lapply(n,getlowsf )),3)
  #highsf <- round(as.numeric(lapply(n,gethighsf)),3)
  lowsf <- highsf <- rep(0,length(ss))
  lownb  <- round(as.numeric(lapply(n,getlownb )),3)
  highnb <- round(as.numeric(lapply(n,gethighnb)),3)
  equantiles <- cbind(1:length(n),aa,aai,ss,lowsf,highsf,lownb,highnb,nsamp,d)
  colnames(equantiles) <- c('index','aa','aai','ss','lowsf','highsf','lownb','highnb','nsamp',
                            paste('q',padzeros(0:100,3),sep=''))
  write.data.frame(equantiles,file=paste('energy_quantile__',etype,'__aa_ss_nb.data',sep=''),pretty=T,
                   meta=paste("META Nbins",dim(d)[1],'Nquantiles',dim(d)[2],
                              'Nssbin',length(unique(ss)),'Nsfbin 5','Nnbbin 5'))
}


#####  write raw distributions

for(binningmode in c("none","aa","ss","aass","nb","aanb","aassnb","sf",
                     "aasf","sssf","aasssf","aasfnb","aasssfnb")) {
  print(binningmode)
  setwd(paste("/users/sheffler/pdb_energies/binned_raw_data/",binningmode,sep=""))
  for(etype in enames) {
    print(paste(" ",etype))
    nsamp <- tapply(pdbe[[etype]],bin[[binningmode]],length)
    d   <- tapply(pdbe[[etype]],bin[[binningmode]],sort)
    n   <- rownames(d)
    aa  <- getaa(n)  
    aai <- res2num[aa]; if(is.na(aai)[1]) aai <- ".."
    ss  <- getss(n)
    lowsf  <- as.character(lapply(n,getlowsf ))
    highsf <- as.character(lapply(n,gethighsf))
    lownb  <- as.character(lapply(n,getlownb ))
    highnb <- as.character(lapply(n,gethighnb))
    metadata <- as.data.frame(list(etype,aa,aai,ss,lowsf,highsf,lownb,highnb,nsamp))
    rownames(metadata) <- rownames(d)
    colnames(metadata) <- c('energy_type','aa','aai','ss','lowsf','highsf','lownb','highnb','nsamp')
    #print(metadata[dim(metadata)[1],])
    for( b in names(d)) {
      #print(paste("  ",b),sep='')
      l <- length(d[[b]])    
      fname <- paste("raw_data_energy_",etype,"_res",b,"_nsamp_",l,".data",sep="")
      write.data.frame(metadata[b,],file=fname,pretty=T)
      write("sorted raw data values:",file=fname,append=T)
      write(d[[b]],file=fname,append=T,ncolumns=1)    
    }
  }
}


#####  write raw distributions INCLUDES SASA

for(binningmode in c("none","aa","ss","aass","nb","aanb","aassnb")) {
  print(binningmode)
  setwd(paste("/users/sheffler/pdb_energies/binned_raw_data_sasa/",binningmode,sep=""))
  for(etype in c("sasa5","sasa14")) {
    print(paste(" ",etype))
    nsamp <- tapply(pdbe[[etype]],bin[[binningmode]],length)
    d   <- tapply(pdbe[[etype]],bin[[binningmode]],sort)
    n   <- rownames(d)
    aa  <- getaa(n)  
    aai <- res2num[aa]; if(is.na(aai)[1]) aai <- ".."
    ss  <- getss(n)
    lowsf  <- as.character(lapply(n,getlowsf ))
    highsf <- as.character(lapply(n,gethighsf))
    lownb  <- as.character(lapply(n,getlownb ))
    highnb <- as.character(lapply(n,gethighnb))
    metadata <- as.data.frame(list(etype,aa,aai,ss,lowsf,highsf,lownb,highnb,nsamp))
    rownames(metadata) <- rownames(d)
    colnames(metadata) <- c('energy_type','aa','aai','ss','lowsf','highsf','lownb','highnb','nsamp')
    #print(metadata[dim(metadata)[1],])
    for( b in names(d)) {
      #print(paste("  ",b),sep='')
      l <- length(d[[b]])    
      fname <- paste("raw_data_energy_",etype,"_res",b,"_nsamp_",l,".data",sep="")
      write.data.frame(metadata[b,],file=fname,pretty=T)
      write("sorted raw data values:",file=fname,append=T)
      write(d[[b]],file=fname,append=T,ncolumns=1)    
    }
  }
}



#############################################


                                        #sfbin <- calcbinsquantile(pdbe$sasafrac,10)
#nbbin <- calcbinsquantile(pdbe$nb+runif(length(pdbe$nb))-0.5,10)
#print(levels(sfbin))
#print(levels(nbbin))

#sf <- pdbe$sasafrac
#nb <- pdbe$nb
#fuzzysf <- sf - runif(length(sf))/1000
#fuzzysf[sf==0] <- 0

#sfqbreaks <- quantilebreaksbyfactor(fuzzysf,aass,10)
#sfbin     <- calcbinsbyfactor(      fuzzysf,aass,sfqbreaks)
#fuzzynb   <- nb+runif(length(nb))
#nbqbreaks <- quantilebreaksbyfactor(fuzzynb,aass,10)
#nbbin     <- calcbinsbyfactor(      fuzzynb,aass,nbqbreaks)

#calcestats <- function(data,enames,bins) {
#  estats <- array(0,c(length(enames),
#                      length(levels(data$aa)),
#                      length(levels(data$ss)),
#                      length(bins),
#                      16,
#                      103))
#  dimnames(estats) <- list(enames,
#                           levels(data$aa),
#                           levels(data$ss),
#                           names(bins),
#                           paste('bin',padzeros(1:16)),
#                           c('count','sdev',paste('q',0:100/100)))
#  names(dimnames(estats)) <- c('etype','res','ss','bins','bin','state')
#  quant <- function(x) quantile(x,0:100/100)
#  for(e in enames) {
#    for(b in names(bins))  {
#      print(paste(e,b))
#      estats[e,,,b,,'count'] <- list2array(tapply(data[,e],list(data$aa,data$ss,bins[[b]]),length))
#      estats[e,,,b,,'sdev' ] <- list2array(tapply(data[,e],list(data$aa,data$ss,bins[[b]]),sd))
#      estats[e,,,b,,3:103  ] <- list2array(tapply(data[,e],list(data$aa,data$ss,bins[[b]]),quant))
#    }
#  }
#  estats[,,,,,'count'] [is.na(estats[,,,,,'count'])] <- 0

#  return(estats)
#}

#nsfbin <- c(16,8,4,2,1)
#nnbbin <- c(1,2,4,8,16)
#pdbebins <- as.list(1:length(nsfbin))

#names(pdbebins) <- c('nb01sf16','nb02sf08','nb04sf04','nb08sf02','nb16sf01')
#for(ii in 1:length(nsfbin)) {
#  sfbreaks <- quantilebreaksbyfactor(fuzzysf,pdbe$aa,nsfbin[ii])
#  nbbreaks <- quantilebreaksbyfactor(fuzzynb,pdbe$aa,nnbbin[ii])
#  sfbin    <- calcbinsbyfactor(      fuzzysf,pdbe$aa,sfbreaks)
#  nbbin    <- calcbinsbyfactor(      fuzzynb,pdbe$aa,nbbreaks)
#  pdbebins[[ii]] <- as.factor(paste('nb',nbbin,'sf',sfbin,sep=''))
#}
#enames <- enames
#data <- pdbe[1:100,]
#bins <- lapply(pdbebins,function(x) x[1:100])

##tmp <- calcestats(data,enames,bins)

#pdbestats <- calcestats(pdbe,enames,pdbebins)
#save(pdbestats,file='pdbestats.rdata')


## todo look at histograms of counts in bins
## figure out why the counts suck so much
## compute bins for each SS/AA combo!!

 
#pdbe$sasafrac <- pmin(1,pdbe$sasafrac)
#pdbe$repe   <- pmin(2,pdbe$repe)
#pdbe$tlje   <- pdbe$atre+pdbe$repe
#pdbe$intrae <- pmin(2,pdbe$intrae)
#pdbe$rese   <- pmin(99,pdbe$rese)
