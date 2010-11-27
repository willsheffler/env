files <- commandArgs()[-1]
files <- files[ -which( substr(files,1,2)=="--" ) ]

c("FILE","SCORE" ,"VDW" ,"ENV" ,"PAIR" ,
  "SSPAIR" ,"HS" ,"RSIGMA" ,"RAMACHANDRAN",
  "FA_ATR", "FA_REP", "FA_SOL", "FA_PAIR",
  "FA_DUN", "FA_PROB", "HB_SRBB", "HB_LRBB",
  "HB_SC", "FA_INTRA", "GSOLT", "SASA", "SASA_PACK",
  "RMSD", "LOOP_TIME", "RELAX_TIME", "RMSD0", "ROUND",
  "description" )

#files <- c("phil_r11di2.score","phil_r11o2f.score","phil_r11tig.score","phil_r11af7.score","phil_r11dtj.score","phil_r11ogw.score","phil_r12reb.score","phil_r11b72.score","phil_r11mky.score","phil_r11r69.score","phil_r11csp.score","phil_r11mla.score","phil_r11shf.score","phil_r11dcj.score","phil_r11n0u.score","phil_r11tif.score")

names <- list()
scores <- list()
rms   <- list()

par(mfrow=c(4,4))
for(f in files) {

  s <- strsplit(f,'/')[[1]]
  #print(s)
  if(length(s)>1) {
    file <- s[[length(s)]]
    path <- paste(s[-length(s)],'/',sep='',collapse='')
  } else {
    file <- s[[1]]
    path <- "./"
  }
  #print(path)
  #print(file)
  #q(save="no")
  if(substr(file,nchar(file)-5,nchar(file))!=".score") {
    print(paste("ignoring non-score file",file))
  } else {
    d <- read.table(paste(path,file,sep=''),header=T,fill=T)
    colnames(d) <- tolower(colnames(d))
    for( n in tolower(c("SCORE" ,"VDW" ,"ENV" ,"PAIR" ,"SSPAIR" ,"HS" ,"RSIGMA" ,
                        "RAMACHANDRAN", "FA_ATR", "FA_REP", "FA_SOL", "FA_PAIR",
                        "FA_DUN", "FA_PROB", "HB_SRBB", "HB_LRBB", "HB_SC",
                        "FA_INTRA", "GSOLT", "SASA", "SASA_PACK", "RMSD",
                        "LOOP_TIME", "RELAX_TIME", "RMSD0", "ROUND","sheet",
                        "cb","rg","co","contact","maxsub","barcode","cst",
                        "chainbreak") ))  {
      if(n %in% names(d)) {
        #print(n)
        d[[n]] <- as.numeric(as.character(d[[n]]))        
        w <-  -which(is.na(d[[n]])) 
        #print(w)
        if(length(w))          
          d <- d[ -which(is.na(d[[n]])) , ]
      }
    }
    
    if("rms" %in% names(d))
      d$rmsd = d$rms
    else if("rmsd" %in% names(d))
      d$rms  = d$rmsd
    
    name <- substr(file,1,nchar(file)-6)
    print(name)
    assign(name,d)
    save(list=c(name),file=paste(path,name,'.rdata',sep=''))

  }

  m <- max(d$rms)
  postscript(file=paste(path,'/',name,"_score_v_rms.ps",sep=''),horizontal=T)
  plot(d$rms,d$score,pch='.',xlim=c(0,m),xlab="RMSD",ylab="Score",
       main=paste("Score vs. RMSD for",name))
  dev.off()

  names  <- append(names,name)
  scores <- append(scores, list(d$score))
  rms    <- append(rms,    list(d$rms))

}


postscript(file="summary_score_vs_rms.ps",horizontal=T)

w <- ceiling(sqrt(length(files)))
h <- ceiling(length(files)/w)
par(mfrow=c(h,w))
par(mar=c(2,2,3,1))
par(mgp=c(2,1,0))
m <- max(as.numeric(lapply(rms,max)))
#print(m)
#print(rms)
for(ii in 1:length(names)) {
  plot(rms[[ii]],scores[[ii]],pch='.',xlab="RMSD",ylab="Score",xlim=c(0,m),
       main=paste("Score vs. RMSD for",names[[ii]]))
}
dev.off()

q(save="no")
