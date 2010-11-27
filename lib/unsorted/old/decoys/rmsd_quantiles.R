source("/users/sheffler/scripts/R/write.data.frame.R")

rmsq <- function(path) {
  s <- read.table(paste(path,'/all.fasc',sep=''),header=T,comment.char="")
  return( quantile(s$rms,0:9/9) )
}

decoydir <- "/users/sheffler/data/decoys/"
q <- array(0,c(0,10))
n <- c()

prots <- c('t196','t199','t205','t206','t223','t224','t232','t234','t249' )
for ( p in prots ) {
  f <- paste(decoydir,'bqian_caspbench',p,sep='/')
  print(f)
  q <- rbind(q,rmsq(f))
  n <- append(n,p)
}

prots <- c('1acp','1b72','1ig5','1r69','1tig','1vii',
           '1aa3','1ail','1bf4','1pgx','1tif','1utg','256b')
for ( p in prots ) {
  f <- paste(decoydir,'kira_fullatombenchmark',p,sep='/')
  print(f)
  q <- rbind(q,rmsq(f))
  n <- append(n,p)
}

prots <- c('1fvk',  '1hb6',  '1hoe',  '1who',  '3il8',  '4rnt')
for ( p in prots ) {
  f <- paste(decoydir,'sraman_nmr',p,sep='/')
  print(f)
  q <- rbind(q,rmsq(f))
  n <- append(n,p)
}

prots <- c('1a4p','1agi','1ar0','1b07','1bmb','1cfy','1d2n','1erv','1pva',
           '1acf','1ahn','1aw2','1be9','1btn','1crb','1dol','1ig5','1qav')
for ( p in prots ) {
  f <- paste(decoydir,'kira_hom',p,sep='/')
  print(f)
  q <- rbind(q,rmsq(f))
  n <- append(n,p)
}

prots <- c('1af7','1csp','1di2','1mky','1n0u','1ogw','1shf','1tig',
           '1b72','1dcj','1dtj','1mla','1r69','1tif','2reb')#,'1o2f'
for ( p in prots ) {
  f <- paste(decoydir,'phil_homolog',p,sep='/')
  print(f)
  q <- rbind(q,rmsq(f))
  n <- append(n,p)
}

rownames(q) <- n
write.data.frame(q,file="sheffler_decoy_set_2_rms_quantiles.txt",pretty=T)
