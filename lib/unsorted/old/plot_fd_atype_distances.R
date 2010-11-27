pd0 <- read.table("/work/sheffler/pdb_energies/pdb.list0.df_atype_distances")
pd1 <- read.table("/work/sheffler/pdb_energies/pdb.list1.df_atype_distances")
pd2 <- read.table("/work/sheffler/pdb_energies/pdb.list2.df_atype_distances")
pd3 <- read.table("/work/sheffler/pdb_energies/pdb.list3.df_atype_distances")
pd4 <- read.table("/work/sheffler/pdb_energies/pdb.list4.df_atype_distances")
pd5 <- read.table("/work/sheffler/pdb_energies/pdb.list5.df_atype_distances")
pd6 <- read.table("/work/sheffler/pdb_energies/pdb.list6.df_atype_distances")
pd7 <- read.table("/work/sheffler/pdb_energies/pdb.list7.df_atype_distances")
pd8 <- read.table("/work/sheffler/pdb_energies/pdb.list8.df_atype_distances")
pd9 <- read.table("/work/sheffler/pdb_energies/pdb.list9.df_atype_distances")
pd <- pd0
pd$V6 <- ( pd0$V6+pd1$V6+pd2$V6+pd3$V6+pd4$V6+
           pd5$V6+pd6$V6+pd7$V6+pd8$V6+pd9$V6) / 10
pd$V7 <- ( pd0$V7+pd1$V7+pd2$V7+pd3$V7+pd4$V7+
           pd5$V7+pd6$V7+pd7$V7+pd8$V7+pd9$V7 )

patk <- read.table('/work/sheffler/data/decoys/patrick_hb_fix/patrick_hb_fix.list.df_atype_distances',header=F)
jd <- read.table("/work/sheffler/data/decoys/baker_jim_his_fix/baker_jim_his_fix_cull.list.df_atype_distances",header=F)

names(patk) <- c("atypenum1","atypenum2","at1","at2","dist","relfreq","count")
names(pd)   <- c("atypenum1","atypenum2","at1","at2","dist","relfreq","count")
names(jd)   <- c("atypenum1","atypenum2","at1","at2","dist","relfreq","count")
jd <- jd[jd$atypenum1<=jd$atypenum2,]

p <- function(t1,t2) {
  idx <- patk$atypenum1==t1 & patk$atypenum2==t2 & patk$dist<5.0
  if(sum(idx)==0)
    return
  plot (patk$dist[idx],patk$count[idx]/sum(patk$count[idx]),
        type='l',xlim=c(1.5,5.0),
        main=paste(t1,t2))
  idx <-   pd$atypenum1==t1 &   pd$atypenum2==t2 &  pd$dist<5.0
  lines(  pd$dist[idx],  pd$count[idx]/sum(pd$count[idx]), col=2)
  idx <-   jd$atypenum1==t1 &   jd$atypenum2==t2 &  jd$dist<5.0
  lines(  jd$dist[idx],  jd$count[idx]/sum(jd$count[idx]), col=3)
}

par(mar=c(3,3,3,1))
par(mfrow=c(3,3))
for(ii in c(8,13,14,15,20,22,25)) {
  p(min(13,ii),max(13,ii))
}

dev.copy2eps(file="/users/sheffler/figures/new_jim_patrick_hb_fix_OH_group_distances.eps")
