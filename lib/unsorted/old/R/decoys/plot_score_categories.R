
plotcategories <- function(decoydir,dataset,prot,xlim=NULL,ylim=NULL,pch='+') {
  dir <- paste(decoydir,dataset,'/',prot,sep='')
  allfile  <- paste(dir,'/all.fasc',sep='')
  gsbrfile <- paste(dir,'/gsbr.fasc',sep='')
  gsgrfile <- paste(dir,'/gsgr.fasc',sep='')
  natfile  <- paste(dir,'/relax_native.fasc',sep='')
  nat  <- NULL
  if(file.exists(allfile)) {
    all  <- read.table(allfile, header=T,comment.char='')
    gsbr <- read.table(gsbrfile,header=T,comment.char='?')
    gsgr <- read.table(gsgrfile,header=T,comment.char='?')
    if( file.exists(natfile) )
      nat <- read.table(natfile,header=T,comment.char='?')
  } else {
    all  <- as.data.frame(list(rms=c(1),score=c(-1)))
    gsbr <- as.data.frame(list(rms=c(1),score=c(-1)))
    gsgr <- as.data.frame(list(rms=c(1),score=c(-1)))
  }
  all  <- all [  all$score <= 0 , ]
  gsbr <- gsbr[ gsbr$score <= 0 , ]
  gsgr <- gsgr[ gsgr$score <= 0 , ]
  I1 <- gsbr$filedataset %in% gsgr$filedataset
  I2 <- gsgr$filedataset %in% gsbr$filedataset
  if( any(I1) | any(I2) )
    print(paste("WARNING: gsgr and gsbr overlap!",dataset,prot))
  
  if(is.null(xlim))
    xlim <- c( 0.0, quantile(all$rms,0.90))
  if(is.null(ylim)) 
    ylim <- c( min(all$score) , max(min(quantile(all$score,.9),0)) )
  
  plot(all$rms,  all$score, xlim=xlim, ylim=ylim, pch='.',
       main=paste(dataset,'\n',prot,'nsamp:',dim(all)[1]))
  points(gsbr$rms, gsbr$score, col="red", pch=pch)
  points(gsgr$rms, gsgr$score, col="blue", pch=pch)
  if(!is.null(nat))
    points(nat$rms, nat$score, col="orange",pch=pch)
    
}

if(0)
  plotcategories('~/data/decoys/','bqian_caspbench','t224')

datasetplot <- function (decoydir,dataset,prots,figdir,xlim,ylim) {
  n <- length(prots)
  y <- floor(sqrt(n))
  x <- ceiling(n/y)
  par(mfrow=c(x,y))
  for( prot in prots )
    plotcategories(decoydir,dataset,prot)
  plotfile <- paste(figdir,'/categories_',dataset,'_relative.ps',sep='')
  dev.print(device=postscript,width=8.5,height=10,plotfile,horizontal=F)
  par(mfrow=c(x,y))
  for( prot in prots ) plotcategories(decoydir,dataset,prot,xlim,ylim)
  plotfile <- paste(figdir,'/categories_',dataset,'_absolute.ps',sep='')
  dev.print(device=postscript,width=8.5,height=10,plotfile,horizontal=F)
}


makeplots <- function() {
  
  decoydir <- "/users/sheffler/data/decoys/"
  figdir   <- "/users/sheffler/project/features/decoydata/fig"
  x11(width=8.5*2/3,height=11*2/3)
  par(mar=c(2,2,3,2))
  
  xlim  <- c(0,8)
  ylim  <- c(-320,0)  
  prots <- c('t196','t199','t205','t206','t223','t224','t232','t234','t249' )
  datasetplot(decoydir,'bqian_caspbench',prots,figdir,xlim,ylim)
  
  
  xlim  <- c(0,14)
  ylim  <- c(-450,0)
  prots <- c('1a32','1acp','1b72','1ig5','1r69','1tig','1vii',
             '1aa3','1ail','1bf4','1pgx','1tif','1utg','256b')
  datasetplot(decoydir,'kira_fullatombenchmark',prots,figdir,xlim,ylim)
  
  xlim  <- c(0,10)
  ylim  <- c(-350,0)  
  prots <- c('1fvk',  '1hb6',  '1hoe',  '1who',  '3il8',  '4rnt')
  datasetplot(decoydir,'sraman_nmr',prots,figdir,xlim,ylim)

  xlim  <- c(0,4)
  ylim  <- c(-1300,0)
  prots <- c('1a4p','1agi','1ar0','1b07','1bmb','1cfy','1d2n','1erv','1pva',
             '1acf','1ahn','1aw2','1be9','1btn','1crb','1dol','1ig5','1qav')
  datasetplot(decoydir,'kira_hom',prots,figdir,xlim,ylim)

  dev.off()

}


## check my scoring...
#dir <- '/users/sheffler/data/decoys/bqian/t196'
#setwd(dir)
#test <- read.table('testrn.fasc',header=T,comment.char='?')
#actual <- read.table('actualrn.fasc',header=T,comment.char='?')
#actual <- actual[-201,]

