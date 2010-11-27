decoy_bruied_polar.data  df_polar  log  relax_native_bruied_polar.data

decoybp  <- read.table('decoy_bruied_polar.data',header=F)[,-1]
nativebp <- read.table('relax_native_buried_polar.data',header=F)[,-1]
nbp <- c('prot','pdb','aa','res','atom',
         'hbe','hbet','dist','cosh','cosa',
         'atomtype','sasa10','sasa7',
         'lk1','lk2','lk4',
         'nb','sasafrac',
         'anb5','anb75','anb10','rms',)
names(decoybp) <- nbp
names(nativebp) <- nbp
save(decoybp,file='decoybp.rdata')
save(nativebp,file='nativebp.rdata')

setwd('/users/sheffler/pdb_buried_polar')
pdbbp <- read.table('pdb_buried_polar.data',header=F)[,-1]
names(pdbbp) <- nbp
save(pdbbp,file='pdbbp.rdata') 

######################

pdbburial <- read.table('pdb_burial.data',header=F)[,-1]
n <- c('prot','pdb','aa','res','atom',
       'hbe','hbet','dist','cosh','cosa',
       'atomtype','sasa14','sasa10','sasa7',
       'lk1','lk2','lk4','nb','sasafrac',
       'anb5','anb75','anb10',
       'lk1e','lk2e','lk4e'  ) 
names(pdbburial) <- n
pdbbp <- pdbburial[,n]

what <- list("","","",integer(0),"",
          double(0),double(0),double(0),double(0),double(0),
          integer(0),
          double(0),double(0),double(0),double(0),double(0),double(0),
          integer(0),double(0), integer(0), integer(0), integer(0),
          double(0),double(0),double(0) )
pdbburial <- scan('pdb_burial.data',what=what)


names(pdbburial) <- n
nbp <- c('pdb','aa','res','atom',
         'hbe','hbet','dist','cosh','cosa',
         'atomtype','nb','sasafrac',
         'anb5','anb75','anb10')
pdbbp <- pdbburial[pdbburial$sasa14==0,nbp]
