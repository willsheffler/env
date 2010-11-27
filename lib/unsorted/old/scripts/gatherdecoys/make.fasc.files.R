
read.fasc.file <- function(fname){
  #print(paste('reading fasc:',fname))
  minimal.fasc.names <- c("filename","score","env","pair","vdw","hs","ss",
                          "sheet","cb","rsigma","hb_srbb","hb_lrbb","rg",
                          "co","rama","bk_tot","fa_atr","fa_rep","fa_sol",
                          "hbsc","fa_dun","fa_pair","fa_prob","gsolt",
                          "rms","rms_min","srms","mxrms","mxlgE",
                          "description")
  head <- read.table(fname,header=T,nrows=1,comment.char='')
  if(names(head)[1]!="filename" | any( ! minimal.fasc.names %in% names(head)) ) {
    print(paste('WARNING: bad header, ignoring file: ',fname))
    return(NULL)
  }
  t <- read.table(fname,header=T,comment.char="",
                  blank.lines.skip=T,fill=T)[,minimal.fasc.names]

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

#pdbdirs <- c(
#"/data/kira/homology_models/homolog_vall//1a4p/1a4p_fold_1dt7_clusterbc_cst_l3/",
#"/data/kira/homology_models/homolog_vall//1a4p/1a4p_fold_1dt7_cst_l3/",
#"/data/kira/homology_models/homolog_vall//1a4p/1a4p_fold_1jwd_clusterbc_cst_l3/",
#"/data/kira/homology_models/homolog_vall//1a4p/1a4p_fold_1jwd_cst_l3/",
#"/data/kira/homology_models/homolog_vall//1a4p/1a4p_fold_1jwd_idscore_cst_l3/"
#)

joint.pdb.list <- function(pdbdirs) {
  pdbs <- as.data.frame(list(name=c(),loc=c()))
  for(ii in 1:length(pdbdirs)) {
    print(paste("reading pdbs from",pdbdirs[ii]))
    tmp  <- dir(pdbdirs[ii],pattern='.pdb')
    name <- sub('[.].*','.pdb',tmp)
    loc  <- rep(pdbdirs[ii],length(name))
    du   <- strsplit(system(paste('du -a',pdbdirs[ii],"| grep 'pdb$' "),T),'\t')
    size <- as.numeric(lapply(du,function(x) as.numeric(x[[1]]) ))
    file <- as.character(lapply(du,function(x) x[[2]]))
    file <- substring( file,nchar(pdbdirs[[ii]])+1 )
    if(length(grep('/',file)))
       file <- file[-grep('/',file)]
    if(length(file)!=length(name)) {
      print("WARNING: len of pdb names and du files differ")
      print(paste("     ",pdbdirs[ii]))
    }
    tmp <- name %in% file
    if(min(tmp)==0) {
      print(paste("WARNING: some pdb names not in du names"))
      print( name [!tmp ])
      name <- name[tmp]
      loc  <- loc [tmp]
    }
    names(size) <- file    
    tosmall <- size[name] < 4
    #tosmall[ !name %in% file ] <- 
    pdbs <- rbind(pdbs,as.data.frame(list(name=name[!tosmall],loc=loc[!tosmall])))
  }
  if(length(levels(pdbs$name))!=length(pdbs$name))
    print(paste("warning:",length(pdbs$name)-length(levels(pdbs$name)),
                "duplicate pdb names!!"))
  first <- tapply(1:length(pdbs$name),pdbs$name,function(x) x[1])
  pdbs  <- pdbs[first,]
  return(pdbs)
}

padleft  <- function(s,l){
  if(nchar(s)!=l)
    return(paste(c(rep(' ',l-nchar(s)),s),sep='',collapse=''))
  return(s)
}
padright <- function(s,l){
  if(nchar(s)!=l)
    return(paste(c(s,rep(' ',l-nchar(s))),sep='',collapse=''))
  return(s)
}

write.data.frame <- function(data,file,pretty=F) {
  print(paste("writting data frame",file))
  n <- colnames(data)
  top <- as.data.frame(as.list(n))
  colnames(top) <- n
  for(ii in 1:length(n))
    top[,ii] <- as.character(top[,ii])
  data <- rbind(top,data)
  if(pretty) {
    for(ii in 1:length(n)) {
      #print(paste("padding column ",ii," (slow; find better way!)"))
      #l <- max(nchar(data[,ii]))
      data[,ii] <- format(data[,ii],justify='left') #sapply(data[,ii],function(x) padright(x,l))
    }
  }
  write.table(data,file=file,quote=F,row.names=F,col.names=F)
}

pickatmost <- function(data,n=200) {
  if(dim(data)[1] <= n)
    return(data)
  return( data[order(runif(dim(data)[1]))[1:n] ,])
}



scorefiles <- scan('log/score_files.txt', list(""))[[1]]
pdbdirs    <- scan('log/decoy_dirs.txt',  list(""))[[1]]

if(0) {
  scorefiles <- scan('~/data/decoys/sraman_nmr/1hb6/log/score_files.txt', list(""))[[1]]
  pdbdirs    <- scan('~/data/decoys/sraman_nmr/1hb6/log/decoy_dirs.txt',  list(""))[[1]]
}

print(scorefiles)
print(pdbdirs)

#fasc <- make.fasc.files(scorefiles,pdbdirs)


print("makeing joint fasc")
fasc <- joint.fasc.file(scorefiles)
print("makeing joint pdb list")
pdbs <- joint.pdb.list(pdbdirs)
Nfasc <- dim(fasc)[1]
Npdb  <- dim(pdbs)[1]
print(paste('starting with',dim(pdbs)[1],'pdbs and',dim(fasc)[1],'scores'))
print("culling scores without pdbs")
fasc <- fasc[fasc$filename %in% pdbs$name,]
print("culling pdbs without scores")
pdbs <- pdbs[pdbs$name %in% fasc$filename,]
print(paste('ending with',dim(pdbs)[1],'pdbs and',dim(fasc)[1],'scores'))
if(length(pdbs$name)!=length(fasc$filename)) {
  print("ERROR: num pdbs not same as num scores")
  return
}
N <- dim(fasc)[1]
print(paste("culled",Npdb-N,'pdbs and',Nfasc-N,'scores'))
print("sorting pdbs and fascs")
pdbs <- pdbs[order(as.character(pdbs$name)),]
fasc <- fasc[order(as.character(fasc$filename)),]
fasc$filename <- paste(as.character(pdbs$loc),as.character(pdbs$name),sep='/')

if(sum(fasc$rms==0) > 10)
  print("WARNING: more than 10 decoys have rms == 0!!")

orms     <- order(fasc$rms)
oscore   <- order(fasc$score)
rmscut   <- quantile(fasc$rms,0.2)
scorecut <- quantile(fasc$score,0.3)
ogsbr    <- order(fasc$score+(fasc$rms<=rmscut)*99999)
ogsgr    <- order(fasc$rms+(fasc$score>=scorecut)*99999)

lowscore <- fasc[oscore[1:200],]
lowrms   <- fasc[orms[1:200],]
gsbr     <- fasc[ogsbr[1:200],]
gsgr     <- fasc[ogsgr[1:200],]

Iolap <- gsbr$filename %in% gsgr$filename
if(any(Iolap))
  print(paste("WARNING",sum(Iolap),"gsgr overlap with gsbr! "))
gsbr <- gsbr[ !Iolap  , ]
#gsgr <- gsgr[ !gsgr$filename %in% gsbr$filename  , ]

rmsbreaks <- 0:10
Nrms <- length(rmsbreaks)-1
rms <- as.list(1:Nrms)
for(ii in 1:Nrms) {
  orms <- order(fasc$score+(!rmsbreaks[ii] <= fasc$rms |
                            !fasc$rms      <= rmsbreaks[ii+1] |
                            !fasc$rms  <= rmsbreaks[ii+1] 
                            )*99999)  
  orms <- orms[ which(rmsbreaks[ii]    <= fasc$rms[orms] &
                      fasc$rms[orms]   <= rmsbreaks[ii+1] &
                      fasc$score[orms] <= scorecut 
                      ) ]
  orms <- orms[1:min(200,length(orms))]
  rms[[ii]] <- fasc[orms,]
}
names(rms) <- paste(0:9,1:10,sep='-')
numrms <- as.numeric(lapply(rms,function(x) length(x[,1])))
names(numrms) <- names(rms)
write.table(rmsbreaks,file='log/rms.breaks',col.names=F,row.names=F)

nat <- NULL
if(file.exists('relax_native.fasc'))
  nat <- read.fasc.file('relax_native.fasc')

print("writting out .fasc files")
write.data.frame(fasc,file='all.fasc',pretty=T)
for(n in c('lowscore','lowrms','gsbr','gsgr')){
  write.data.frame(get(n),file=paste(n,'.fasc',sep=''),pretty=T)
  if(!file.exists(n)) dir.create(n)
}
for(ii in 1:Nrms){
  write.data.frame(rms[[ii]],file=paste('rms',ii-1,'.fasc',sep=''))
  if(!file.exists(paste('rms',ii-1,sep='')))
    dir.create(paste('rms',ii-1,sep=''))
}

print("generating plots")
pch <- '+'
nat <- NULL

postscript(file="plot_gsbr_gsgr.ps",horizontal=F)
xlim <- c( 0.0, quantile(fasc$rms,0.90))
ylim <- c( min(fasc$score) , min(0,quantile(fasc$score,0.90)))
if(!is.null(nat))
  ylim <- c(min(ylim[1],nat$score),max(ylim[2],quantile(nat$score,.9)))
plot(  fasc$rms, fasc$score, xlim=xlim, ylim=ylim, pch=pch)
points(gsbr$rms, gsbr$score, col="red", pch=pch)
points(gsgr$rms, gsgr$score, col="blue", pch=pch)
if(!is.null(nat))
  points(nat$rms, nat$score, col="orange",pch=pch)
dev.off()

postscript(file="plot_lowscore_lowrms.ps",horizontal=F)
sinr <- lowscore$filename %in% lowrms  $filename
rins <- lowrms  $filename %in% lowscore$filename
plot(  fasc$rms, fasc$score, xlim=xlim, ylim=ylim, pch=pch)
points(lowrms  $rms[!rins], lowrms  $score[!rins], col="red",    pch=pch)
points(lowscore$rms[!sinr], lowscore$score[!sinr], col="blue",   pch=pch)
points(lowscore$rms[ sinr], lowscore$score[ sinr], col="violet", pch=pch)
if(!is.null(nat)) 
  points(nat$rms, nat$score, col="orange",pch=pch)
dev.off()

postscript(file="plot_rms_categories.ps",horizontal=F)
plot(  fasc$rms, fasc$score, xlim=xlim, ylim=ylim, pch=pch)
mn <- ceiling(min(fasc$rms)+0.000000001)
mx <- ceiling(xlim[2])
print(paste("mn",mn,"mx",mx))
mn <- max(mn,1)
mx <- min(mx,12)
colors <- c(rep('red',mn-1),rainbow(mx-mn+1),rep('violet',12-mx))
print(colors)
for(ii in 1:Nrms){
  if(length(rms[[ii]]$rms))
    points(rms[[ii]]$rms,rms[[ii]]$score, col=colors[ii],pch=pch)
}
if(!is.null(nat))
  points(nat$rms, nat$score, col="orange",pch=pch)
dev.off()
