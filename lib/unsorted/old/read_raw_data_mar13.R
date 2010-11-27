atmp <- read.table("/users/sheffler/project/atom_distances/atype.data")
anames <- list()
for(ii in unique(as.numeric(atmp[,1])))
  anames[[ii]] <- as.character(atmp[atmp[,1]==ii,3])
atypes <- c()
for(ii in unique(as.numeric(atmp[,1])))
  atypes[ii] <- as.character(atmp[atmp[,1]==ii,2])[1]
aindex <- 1:25
names(aindex) <- atypes

pd <- array(0,c(10,10,25,25,1000))
dimnames(pd) <- list(c('anb000-050','anb051-100','anb101-150','anb151-200','anb201-250',
                       'anb251-300','anb301-350','anb351-400','anb401-450','anb451-500' ),
                     c('sep0','sep1','sep2','sep3-4','sep5-8','sep9-16',
                       'sep17-32','sep33-64','sep65-128','sep129-256'),
                     atypes, atypes, 1:1000/100-0.005  )
tmp1 <- read.table("/users/sheffler/w/pdb_energies/pdb_atype_dist.data0",header=F)
tmp2 <- read.table("/users/sheffler/w/pdb_energies/pdb_atype_dist.data1",header=F)
tmp2 <- read.table("/users/sheffler/w/pdb_energies/pdb_atype_dist.data2",header=F)
tmp2 <- read.table("/users/sheffler/w/pdb_energies/pdb_atype_dist.data3",header=F)
names(tmp1) <- c( 'anb','sep','at1','at2','dist','count' )
names(tmp2) <- c( 'anb','sep','at1','at2','dist','count' )

for(tmp in list(tmp1,tmp2,tmp3,tmp4)) {
  attach(tmp)
  idx1 <- ceiling(anb/60)
  idx2 <- ceiling(log2(sep+1))+1 # from 0!
  idx3 <- aindex[at1] 
  idx4 <- aindex[at2] 
  idx5 <- ceiling(dist*100) 
  d1 <- dim(pd)[2]*dim(pd)[3]*dim(pd)[4]*dim(pd)[5]
  d2 <-            dim(pd)[3]*dim(pd)[4]*dim(pd)[5]
  d3 <-                       dim(pd)[4]*dim(pd)[5]
  d4 <-                                  dim(pd)[5]
  d5 <-                                           1
  idx  <- (idx1-1)*d1 + (idx2-1)*d2 + (idx3-1)*d3 + (idx4-1)*d4 + (idx5-1)*d5 + 1
  pd[idx] <- count
  
  for(ii in 1:dim(tmp1)[1]) {
    if(ii%%1000==0)
      print(ii)
    pd[  ,
        ,
       aindex[at1[ii]] ,
        ,
        ]  <- dist[ii]
  }
}



