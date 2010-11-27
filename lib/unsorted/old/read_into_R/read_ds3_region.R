decoydir <- "/users/sheffler/data/decoys"
decoys <- list()
decoys$bqian_caspbench <- c("t196","t199","t205","t206","t223","t224","t232","t234","t249")
decoys$phil_homolog <-    c("1af7","1csp","1di2","1mky","1n0u","1ogw","1shf","1tig",
                            "1b72","1dcj","1dtj","1mla","1r69","1tif")#,"2reb" )
decoys$kira_hom <-        c("1a4p","1agi","1ar0","1b07","1bmb","1cfy","1d2n","1erv","1pva",
                            "1acf","1ahn","1aw2","1be9","1btn","1crb","1dol","1ig5","1qav" )
decoys$sraman_nmr <-      c("1fvk","1hb6","1hoe","1i27","1who","2int","3il8","4rnt")
decoys$kira_bench <-      c("1acp","1b72","1ig5","1r69","1tig","1vii",
                            "1aa3","1ail","1bf4","1pgx","1tif","1utg")#,"256b","1a32")




calc.decoy.region.scores <- function(input.categories,sets) {
  if(missing(sets))
    sets <- names(decoys)
  for(set in sets) {
    regionscores <- data.frame()
    for(prot in decoys[[set]]) {
      if(missing(input.categories)) {
        rand <- dir(paste(decoydir,set,prot,'data',sep='/'),pattern='rand_.*_region.data$')
        rand <- as.character(lapply(rand,function(x)substring(x,1,nchar(x)-12)))
        if(length(rand)==0)
          categories <- c('gsbr','gsgr','relax_native')
        else
          categories <- c('gsbr','gsgr','relax_native',rand[1])
        if(set=='phil_homolog')
          categories <- c('all',categories)
      } else {
        categories <- input.categories
      }
      for(category in categories) {
        #print(paste(set,prot,category))
        fascfile <- paste(decoydir,'/',set,'/',prot,'/data/',category,'_score.fasc',sep='')
        regfile  <- paste(decoydir,'/',set,'/',prot,'/data/',category,'_region.data',sep='')
        if(!file.exists(regfile)){
          print(paste('no region.data for',set,prot,category))
        } else if(0) {
          #tmp <- read.table(regfile,nrows=1,header=F)
          #if(dim(tmp)[2]!=34)
          #  print(paste(dim(tmp)[2],regfile))
                                        #print(regfile
          reg <- read.table(regfile,nrows=-1,header=F)[,-1]
          #print(paste(dim(reg)[2],regfile))
          if(dim(reg)[2]!=33) {
            print(paste("bad dim on",set,prot,category,dim(reg)[2],collapse=' '))
          } else {
            print(paste(set,prot,category,dim(reg)[2],dim(reg)[1]))
            names(reg) <- c('prot','pdb','chain','res','aa','ss','nb',
                            'atr1','res1','spk1','bup1','rms1',
                            'rms25','rms50','rms75','rms100','rms125',
                            'r7', 'nb7', 'atr7', 'res7', 'spk7', 'bup7', 'rmsabs7','rmsrel7',
                            'r10','nb10','atr10','res10','spk10','bup10','rmsabs10','rmsrel10'
                            )
            #print(dim(reg))
            #print( reg)
            reg$prot <- prot
            reg$set <- set
            reg$cat <- category
            regionscores <- rbind(regionscores,reg)
          }
        }
      }
    }
    save(regionscores,file=paste('~/project/features/dr.',set,'.ds3.region.rdata',sep=''))
    rm(regionscores)
  }
}

calc.decoy.region.scores(c('gsgr','gsbr'),sets=c('sraman_nmr')) 
