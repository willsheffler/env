 
read.fasc.file <- function(fname){
  #print(paste('reading fasc:',fname))
  minimal.fasc.names <- c('filename','score','rms')
  head <- read.table(fname,header=T,nrows=1,comment.char='')
  if(names(head)[1]!="filename" | any( ! minimal.fasc.names %in% names(head)) ) {
    print(paste('WARNING: bad header, ignoring file: ',fname))
    return(NULL)
  }
  t <- read.table(fname,header=T,comment.char="",
                  blank.lines.skip=T,fill=T)[,minimal.fasc.names]
  if(any(is.na(t[,1]))) {
    print(paste("WARNING: rows NA, ignoring file",fname))
    return(NULL)
  }
  
  falsehead <- which(is.na(sapply(t$rms,as.numeric)))
  if( length(falsehead) > 0) {
    t <- t[-falsehead,]
    for( n in minimal.fasc.names[-c(1,length(minimal.fasc.names))] ) {
      t[,n] <- as.numeric(t[,n])
    }
  }

  t <- t[ t$score <= 0 , ]
  
  tmp <- which(apply(t,1,function(x) sum(is.na(x)))>0)
  if( length(tmp) > 0) {
    print(paste("WARNING:NA lines in score file",fname))
    nalines <- c()
    for(w in tmp) {
      nalines <- c(nalines,w)
      if(w > 1) nalines <- c(nalines,w-1)
      if(w < dim(t)[1]) nalines <- c(nalines,w+1)
    }
    print(paste("     removing line ",nalines))
    t <- t[-nalines,]
  }
  
  Ibad <-  t$rms < 0
  if( any(Ibad) ) {
    print(paste("WARNING: some lines bad in score file",fname))
    print(paste("     lines: ",which(Ibad)))
    t <- t[ Ibad , ]
  }
  for(ii in 1:length(minimal.fasc.names)) {
    I <- is.na(t[,ii])
    if(any(I))
      t <- t[!I,]
  }
  #print(paste(fname,sum(is.na(t))))
  names(t) <- minimal.fasc.names
  #names(t)[names(t)=="mx."] <- "mx#"
  return(t)
}

#for(file in files) { read.fasc.file(file) }

joint.fasc.file <- function(scorefiles) {
  ii <- 0 ; scores <- NULL
  while(is.null(scores)) {
    ii <- ii+1
    scores <- read.fasc.file(scorefiles[ii])
  }
  for(file in scorefiles){ #[-1]) { # testing!
    print(paste("reading fasc:",file))
    tmp <- read.fasc.file(file)
    if(!is.null(tmp)) {
      print(paste('adding fasc: ',file))
      filename <- as.character(tmp$filename)
      pdb <- as.character(lapply(strsplit(gsub('//','/',filename),'/'),function(x) x[length(x)]))
      #path <- (strsplit(gsub('//','/',file),'/'))[[1]]
      #path <- paste(path[-length(path)],collapse='/')
      path <- substr(file,1,nchar(file)-5)
      filename <- paste(path,pdb,sep='/')
      tmp$filename <- filename
      scores <- rbind(scores,tmp)
    } 
  }
  scores$filename <- as.factor(sub('[.].*','.pdb',as.character(scores$filename)))  
  first <- tapply(1:length(scores$filename),scores$filename,function(x) x[1])
  scores <- scores[first,]
  if(length(first) != length(scores$filename))
    print(paste("removing",length(scores$filename)-length(first),
                "duplicate 'filename's from score files"))
  return(scores)
}

scorefiles <- as.character(read.table('~/data/decoys/ds2.score.files.whip')$V1)
ds2.score <- joint.fasc.file(scorefiles)
#rownames(ds2.score) <- as.character(ds2.score$filename)
#ds2.score <- ds2.score[,-1]
#save(ds2.score,file="ds2.score.dataset2.score.rms.rdata")
 
