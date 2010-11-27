
# include some basic utilities from a separate file
source('~/scripts/R-utils.R') #CHANGEME!

# some dataset must be loaded or created
# energy.data is an R 'data.frame' with colums for
# each energy types and some metadata
# type names(energy.data) to see the names of the columes
#setwd('~/project/features/energy_v_pdb')
load("~/data/rdata/project_features_energy_v_pdb_pe.pdb.energies.rdata")
load("~/data/rdata/project_features_energy_v_pdb_de.decoyset2.energy.rdata")

train.data <- pe
#test.data  <- pe
test.data  <- de

##########################################
##    FIX ME! BRK on PDB, BIN on DECOY!!
##########################################


# impose cutoffs for burial measures
train.data$sasafrac <- pmin(train.data$sasafrac,1)
train.data$nb       <- pmax(train.data$nb,6)
test.data$sasafrac  <- pmin(test.data$sasafrac,1)
test.data$nb        <- pmax(test.data$nb,6)
# random vector to spread out discrete "nb" values
r <- runif(dim(train.data)[1])-0.5


# (training) discrets bins these are R "factors"
aa    <- as.character(train.data$aa)
ss    <- as.character(train.data$ss)
aass  <- paste(aa,'ss',ss,sep='')
aaf   <- as.factor(aa)
ssf   <- as.factor(ss)
aassf <- as.factor(aass)

# these are break vectors for the scalar
# "brk" is a list with names
# breaks functions are in R-utils.R
# quantilebreaks(data,Nbins): returns Nbins+1 breaks for vector data
brk <- list()
brk$nb      <- c(5:30+0.5)
brk$sf      <- quantilebreaks(train.data$sasafrac,25)
brk$aasf    <- quantilebreaksbyfactor(train.data$sasafrac,aaf,  25)
brk$sssf    <- quantilebreaksbyfactor(train.data$sasafrac,ssf,  25)
brk$aasssf  <- quantilebreaksbyfactor(train.data$sasafrac,aassf,25)


# more brk's are computed here
# these are 5 bins on sasa frac (to be further divided in 5 bins each by nb)
# these breaks/bins are used only as intermediates to do binning on nb & sf
brk$aasf5    <- quantilebreaksbyfactor(train.data$sasafrac,aaf,  5)
brk$aasssf5  <- quantilebreaksbyfactor(train.data$sasafrac,aassf,5)
bin.aasf5    <- calcbinsbyfactor(train.data$sasafrac,aaf,  brk$aasf5  )
bin.aasssf5  <- calcbinsbyfactor(train.data$sasafrac,aassf,brk$aasssf5)
aasf5        <- as.factor(paste(aa,  'sf',as.character(bin.aasf5  ),sep=''))
aasssf5      <- as.factor(paste(aass,'sf',as.character(bin.aasssf5),sep=''))

# these breaks/bins are for nb & sf
brk$aanb5    <- quantilebreaksbyfactor(train.data$nb+r,aasf5,  5)
brk$aassnb5  <- quantilebreaksbyfactor(train.data$nb+r,aasssf5,5)
bin.aanb5    <- calcbinsbyfactor(train.data$nb,aasf5,  brk$aanb5  )
bin.aassnb5  <- calcbinsbyfactor(train.data$nb,aasssf5,brk$aassnb5)



# (testing) discrets bins these are R "factors"
aa    <- as.character(test.data$aa)
ss    <- as.character(test.data$ss)
aass  <- paste(aa,'ss',ss,sep='')
aaf   <- as.factor(aa)
ssf   <- as.factor(ss)
aassf <- as.factor(aass)



# these are the bins
# "bin" is a named list
# bin factor levels have names like this:
# RESssXsfYnbZ, RES= 3 letter res code; X = H,S, or L; Y = sf bin (01-25 or 01-05), X=nb bin (01-25 or 01-05)
# actual bins differ for each aa,ss. actual breaks for bin can be found in the "brk" list
# calcbins is in R-util.R
bin <- list()
bin$none     <- as.factor(rep('...ss.sf..nb..',length(aass)))
bin$aa       <- as.factor(paste(aa,'ss.sf..nb..',sep=''))
bin$ss       <- as.factor(paste('...ss',ss,'sf..nb..',sep=''))
bin$aass     <- as.factor(paste(aass,'sf..nb..',sep=''))
bin$nb       <- as.factor(paste('...ss.sf..nb', as.character(calcbins(test.data$nb,brk$nb)),sep=''))
bin$aanb     <- as.factor(paste(aa,'ss.sf..nb', as.character(calcbins(test.data$nb,brk$nb)),sep=''))
bin$aassnb   <- as.factor(paste(aass,'sf..nb',  as.character(calcbins(test.data$nb,brk$nb)),sep=''))
bin$sf       <- as.factor(paste('...ss.sf',     as.character(calcbins(test.data$sasafrac,brk$sf)),'nb..',sep=''))
bin$aasf     <- as.factor(paste(aa,'ss.sf',as.character(calcbinsbyfactor(test.data$sasafrac,aaf,brk$aasf)),'nb..',sep=''))
bin$sssf     <- as.factor(paste('...ss',ss,'sf',as.character(calcbinsbyfactor(test.data$sasafrac,ssf,brk$sssf)),'nb..',sep=''))
bin$aasssf   <- as.factor(paste(aass,'sf',as.character(calcbinsbyfactor(test.data$sasafrac,aassf,brk$aasssf)),'nb..',sep=''))
bin$aasfnb   <- as.factor(paste(aa,'ss.sf',as.character(bin.aasf5),  'nb',as.character(bin.aanb5),  sep=''))
bin$aasssfnb <- as.factor(paste(aass, 'sf',as.character(bin.aasssf5),'nb',as.character(bin.aassnb5),sep=''))



# little utility functions to parse out info from the bin factor names
getaa     <- function(f) return(substr(f,1,3))
getss     <- function(f) return(substr(f,6,6))


getlowsf  <- function(f,offset=0){
  if(substr(f, 9,10)=="..")
    return("..")
  if(substr(f,13,14)=="..")
    return(brk$sf[as.numeric(substr(f, 9,10))+offset])
  return(brk$aasf5[[substr(f,1,3)]][as.numeric(substr(f, 9,10))+offset])
}
getlowsf  <- function(f,offset=0){
  if(substr(f,09,10)=="..")
    return("..")
  if(substr(f,13,14)=="..")
    return(brk$sf[as.numeric(substr(f,09,10))+offset])
  if(substr(f,6,6)=='.') {
    return(brk$aasf5[[paste(substr(f,1,3),sep='')]]
           [as.numeric(substr(f,13,14))+offset])
  } else {
    return(brk$aasssf5[[paste(substr(f,1,6),sep='')]]
           [as.numeric(substr(f,13,14))+offset])
  }
}

gethighsf <- function(f){
  getlowsf(f,1)
}


#brk$nb      <- c(5:30+0.5)
#brk$aanb5    <- quantilebreaksbyfactor(train.data$nb+r,aasf5,  5)
#brk$aassnb5  <- quantilebreaksbyfactor(train.data$nb+r,aasssf5,5)


getlownb  <- function(f,offset=0){
  if(substr(f,13,14)=="..")
    return("..")
  if(substr(f, 9,10)=="..")
    return(brk$nb[as.numeric(substr(f,13,14))+offset])
  if(substr(f,6,6)=='.') {
    return(brk$aanb5[[paste(substr(f,1,3),substring(f,7,10),sep='')]]
           [as.numeric(substr(f,13,14))+offset])
  } else {
    return(brk$aassnb5[[paste(substr(f,1,10),sep='')]]
           [as.numeric(substr(f,13,14))+offset])
  }
}
gethighnb <- function(f){
  getlownb(f,1)
}


#enames <- c('atre','repe','sole','probe','dune',
#            'intrae','hbe','paire','rese','tlje','sasapack',
#            'sasa14','sasa5','atrrepsol')

test.data$atrrepsol <- test.data$atre + test.data$repe + test.data$sole
enames <- c('atrrepsol') # example new energy

#####  write raw distributions

# loop over different sets of bins
for(binningmode in c("none","aa","ss","aass","nb","aanb","aassnb","sf",
                     "aasf","sssf","aasssf","aasfnb","aasssfnb")) {
  print(binningmode)
 
  # make and set WD
  wd <- paste("/users/sheffler/data/ENERGY_BINNED/ds2_binned_raw_data_atrrepsol/",binningmode,sep="")
  system(paste("mkdir -p",wd))
  setwd(wd)
  rm(wd)
  
  # loop over energy types (must be col names of train/test.data)
  for(etype in enames) {

    # don't bin sasa on anything involving sasafrac
    # (sasafrac = sasa14/max(sasa14) per sesidue)
    if( ! etype %in% c('sasa5','sasa14') |
       binningmode %in% c("none","aa","ss","aass","nb","aanb","aassnb") ) {
      
      print(paste(" ",etype))
      
      # this is where the actual binning is done
      d   <- tapply(test.data[[etype]],bin[[binningmode]],sort)
      
      # this is all "metadata"
      nsamp <- tapply(test.data[[etype]],bin[[binningmode]],length)
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
      colnames(metadata) <- c('energy_type','aa','aai','ss',
                              'lowsf','highsf','lownb','highnb','nsamp')
                                        #print(metadata[dim(metadata)[1],])
      # loop over the different bins
      for( b in names(d)) {
                                        #print(paste("  ",b),sep='')
        # write out metadata
        l <- length(d[[b]])    
        fname <- paste("raw_data_energy_",etype,"_res",b,"_nsamp_",l,".data",sep="")
        write.data.frame(metadata[b,],file=fname,pretty=T)
        write("sorted raw data values:",file=fname,append=T)

        # this is where the actual data is written
        write(d[[b]],file=fname,append=T,ncolumns=1)    
      }

    #}
      
    }
  }
}
  
