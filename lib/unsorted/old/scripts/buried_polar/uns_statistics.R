######## load raw data ########

# setwd('/users/sheffler/project/packing_index/rdata')
# pdb.sasa0 <- read.table('../pdb_uns_list.txt',header=F)[,3:14]
# names(pdb.sasa0) <- c('pdb','aa','res','atom','hbond','atomname',
#                      'atomtype','sasa14','sasa10','sasa7','nb','sasafrac')
# save(pdb.sasa0,file='pdb.sasa0.rdata')
# kiranathbond <- read.table('../kira_decoys_native_hbond.txt',header=F)[,3:14]
# names(kiranathbond) <- c('pdb','aa','res','atom','hbond','atomname',
#                          'atomtype','sasa14','sasa10','sasa7','nb','sasafrac')
# save(kiranathbond,file='kiranathbond.rdata')



load('pdbdata.rdata')
load('pdbuns.rdata')

aanames1 <- c("A","R","N","D","C","Q","E","G","H","I",
              "L","K","M","F","P","S","T","W","Y","V")
aanames3 <- c("ALA","ARG","ASN","ASP","CYS","GLN","GLU","GLY","HIS","ILE",
              "LEU","LYS","MET","PHE","PRO","SER","THR","TRP","TYR","VAL")


resname2rosettaid <- 1:20
names(resname2rosettaid) <- c("ALA","CYS","ASP","GLU","PHE","GLY","HIS",
                              "ILE","LYS","LEU","MET","ASN","PRO","GLN",
                              "ARG","SER","THR","VAL","TRP","TYR")


# total residue counts
rescount <- tapply(pdbdata$aa,pdbdata$aa,function(x) length(x))
resmeannb    <- tapply(pdbdata$nb,pdbdata$aa,mean)
resmeansasafrac <- tapply(pdbdata$sasafrac,pdbdata$aa,mean)

par(mfrow=c(3,2))

# for each residue, what fraction have UNS BB O and BB H?
bbhfrac <- tapply(pdbuns$atomtype,pdbuns$aa,function(x) sum(x==25) ) / rescount
bbofrac <- tapply(pdbuns$atomtype,pdbuns$aa,function(x) sum(x==20) ) / rescount
par(las=2)
barplot( bbhfrac, main="backbone H UNS frac", ylim=c(0,0.16))
barplot( bbofrac, main="backbone O UNS frac", ylim=c(0,0.16))

par(las=0)
# is this just related to burial?
plot( bbhfrac , resmeannb , type='n',  main="BB H UNS frac vs. Mean NB" )
text( bbhfrac , resmeannb , levels(pdbuns$aa) )

plot( bbofrac , resmeannb , type='n',  main="BB O UNS frac vs. Mean NB" )
text( bbofrac , resmeannb , levels(pdbuns$aa) )

plot( bbhfrac , resmeansasafrac , type='n',  main="BB H UNS frac vs. Mean SASAfrac" )
text( bbhfrac , resmeansasafrac , levels(pdbuns$aa) )

plot( bbofrac , resmeansasafrac , type='n',  main="BB O UNS frac vs. Mean SASAfrac" )
text( bbofrac , resmeansasafrac , levels(pdbuns$aa) )

dev.print(device=postscript,file='../fig/uns_bb.ps',horizontal=F)


# what are the frequencies of UNS for all hbond don/acc types?

aa.atomname <- as.factor(paste(as.character(pdbuns$aa),as.character(pdbuns$atomname)))
aa.atomname.rescount <- rescount[substring(levels(aa.atomname),1,3)]
aa.atomname.unscount <- tapply(aa.atomname,aa.atomname,function(x) length(x))
aa.atomname.unsfrac <- aa.atomname.unscount / aa.atomname.rescount

par(mfrow=c(1,1))
par(las=2)
par(mar=c(4,6,4,2))
barplot( rev(aa.atomname.unsfrac) , main="UNS Frac For All Hbond Atom Names" ,
         horiz=T, xlab="Fraction Buried-UNS / Total Cases",
        xlim=c(0,0.3), cex.names=0.8 )
dev.print(device=postscript,file='../fig/uns_all_atomnames.ps',horizontal=F)

###################################################
# PDB burried hbond stats
###################################################

load('pdb.sasa0.rdata')
load('pdbdata.rdata')

# all residue stats (not just buried)
rescount <- tapply(pdbdata$aa,pdbdata$aa,function(x) length(x))
resmeannb    <- tapply(pdbdata$nb,pdbdata$aa,mean)
resmeansasafrac <- tapply(pdbdata$sasafrac,pdbdata$aa,mean)

###########################3
# overall stats
overall.stats <- function() {
  x <- sum(pdb.sasa0$hbond >= -0.01) / dim(pdb.sasa0)[1]
  print(paste("Fraction of buried polar atoms with uns hbonds",x))
  
  I <- pdb.sasa0$atomtype==20
  x <- sum(pdb.sasa0$hbond[I] >= -0.01) / sum(I)
  print(paste("Fraction of buried backbone C=O atoms with uns hbonds",x))
  x <- sum(pdb.sasa0$hbond[I] >= -0.01) / sum(rescount)
  print(paste("Fraction of all backbone C=O atoms with uns hbonds",x))
  
  I <- pdb.sasa0$atomtype==25
  x <- sum(pdb.sasa0$hbond[I] >= -0.01) / sum(I)
  print(paste("Fraction of buried backbone H-N atoms with uns hbonds",x))
  x <- sum(pdb.sasa0$hbond[I] >= -0.01) / sum(rescount)
  print(paste("Fraction of all backbone H-N atoms with uns hbonds",x))
}

###########################
# backbone uns frac
Ibbh <- pdb.sasa0$atomtype == 25
Ibbo <- pdb.sasa0$atomtype == 20
pdb.sasa0.bbh <- pdb.sasa0[Ibbh,]
pdb.sasa0.bbo <- pdb.sasa0[Ibbo,]

bbh.count    <- tapply(pdb.sasa0.bbh$hbond,pdb.sasa0.bbh$aa,function(x) length(x) )
bbh.unscount <- tapply(pdb.sasa0.bbh$hbond,pdb.sasa0.bbh$aa,function(x) sum(x >= -0.01) )
bbh.unsfrac  <- bbh.unscount / bbh.count
bbh.buryfrac <- bbh.count / rescount
bbh.meannb   <- tapply(pdb.sasa0.bbh$nb,pdb.sasa0.bbh$aa,mean)

bbo.count    <- tapply(pdb.sasa0.bbo$hbond,pdb.sasa0.bbo$aa,function(x) length(x) )
bbo.unscount <- tapply(pdb.sasa0.bbo$hbond,pdb.sasa0.bbo$aa,function(x) sum(x >= -0.01) )
bbo.unsfrac  <- bbo.unscount / bbo.count
bbo.buryfrac <- bbo.count / rescount
bbo.meannb   <- tapply(pdb.sasa0.bbo$nb,pdb.sasa0.bbo$aa,mean)

par(mfrow=c(3,2))
par(las=2)
barplot( bbh.unsfrac, main="backbone H UNS frac", ylim=c(0,0.23))
barplot( bbo.unsfrac, main="backbone O UNS frac", ylim=c(0,0.23))

par(las=0)
# is this just related to burial?
plot( bbh.unsfrac , bbh.meannb , type='n',  main="BB H UNS frac vs. Mean NB" )
text( bbh.unsfrac , bbh.meannb , names(bbh.unsfrac) )

plot( bbo.unsfrac , bbo.meannb , type='n',  main="BB O UNS frac vs. Mean NB" )
text( bbo.unsfrac , bbo.meannb , names(bbo.unsfrac) )

plot( bbh.unsfrac , bbh.buryfrac , type='n',  main="BB H UNS frac vs. buried frac" )
text( bbh.unsfrac , bbh.buryfrac , names(bbh.unsfrac) )

plot( bbo.unsfrac , bbo.buryfrac , type='n',  main="BB O UNS frac vs. buried frac" )
text( bbo.unsfrac , bbo.buryfrac , names(bbo.unsfrac) )

dev.print(device=postscript,file='../fig/uns_bb.ps',horizontal=F)



# all atom uns stats

aa.atomname <- as.factor(paste(as.character(pdb.sasa0$aa),
                               as.character(pdb.sasa0$atomname)))
aa.atomtype <- as.factor(paste(as.character(pdb.sasa0$aa),
                               as.character(pdb.sasa0$atomtype)))
aa.atomname.count <- tapply(aa.atomname,aa.atomname,function(x) length(x))
aa.atomtype.count <- tapply(aa.atomtype,aa.atomtype,function(x) length(x))
aa.atomname.unscount <- tapply(pdb.sasa0$hbond,aa.atomname,function(x) sum(x > -0.01))
aa.atomtype.unscount <- tapply(pdb.sasa0$hbond,aa.atomtype,function(x) sum(x > -0.01))
aa.atomname.unsfrac <- aa.atomname.unscount / aa.atomname.count
aa.atomtype.unsfrac <- aa.atomtype.unscount / aa.atomtype.count

par(mfrow=c(1,1))
par(las=2)
par(mar=c(4,6,4,2))
barplot( rev(aa.atomname.unsfrac) , main="UNS Frac For All Hbond Atom Names" ,
         horiz=T, xlab="Fraction Buried-UNS / Total Buried",
         cex.names=0.7)#, xlim=c(0,0.5) )
dev.print(device=postscript,file='../fig/unsfrac_all_atomnames.ps',horizontal=F)
barplot( rev(aa.atomtype.unsfrac) , main="UNS Frac For All Hbond Atom Types" ,
         horiz=T, xlab="Fraction Buried-UNS / Total Buried",
         cex.names=0.7)#, xlim=c(0,0.5) )
dev.print(device=postscript,file='../fig/unsfrac_all_atomtypes.ps',horizontal=F)

# stats on just side chain atoms
I <- nchar(levels(aa.atomname)) > 5
sc.atomnames <- levels(aa.atomname)[I]
tmp <- substring(levels(aa.atomtype),5,6)
I <- tmp != "20" & tmp != "25"
sc.atomtypes <- levels(aa.atomtype)[I]

par(mfrow=c(1,1))
par(las=2)
par(mar=c(4,6,4,2))
barplot( rev(aa.atomname.unsfrac[sc.atomnames]) ,
         main="UNS Frac For All Hbond Atom Names" ,
         horiz=T, xlab="Fraction Buried-UNS / Total Buried",
         cex.names=1)#, xlim=c(0,0.5) )
dev.print(device=postscript,file='../fig/unsfrac_sidechain_atomnames.ps',horizontal=F)
barplot( rev(aa.atomtype.unsfrac[sc.atomtypes]) ,
         main="UNS Frac For All Hbond Atom Types" ,
         horiz=T, xlab="Fraction Buried-UNS / Total Buried",
         cex.names=1)#, xlim=c(0,0.5) )
dev.print(device=postscript,file='../fig/unsfrac_sidechain_atomtypes.ps',horizontal=F)


#####################
# fraction of buried UNS each atom is responsible for

total.buried.uns <- sum(pdb.sasa0$hbond >= -0.01 )

aa.atomname <- as.factor(paste(as.character(pdb.sasa0$aa),
                               as.character(pdb.sasa0$atomname)))
aa.atomname.unscount <- tapply(pdb.sasa0$hbond,aa.atomname,function(x) sum(x >= -0.01))
aa.atomname.fractotal <- aa.atomname.unscount / total.buried.uns

I <- nchar(levels(aa.atomname)) > 5
sc.atomnames <- levels(aa.atomname)[I]
bbo.atomnames <- paste(aanames3,"O")
bbh.atomnames <- paste(aanames3,"H")


par(mfrow=c(3,1))
par(las=2)
par(mar=c(8,6,4,2))
barplot((aa.atomname.fractotal[bbo.atomnames]) ,
        main="BBO Fraction buried uns / total uns" ,
        xlab="",#BB O Fraction Buried-UNS / Total Buried UNS",
        cex.names=1)#, xlim=c(0,0.5) )
barplot((aa.atomname.fractotal[bbh.atomnames]) ,
        main="BBH Fraction buried uns / total uns" ,
        xlab="",#BB O Fraction Buried-UNS / Total Buried UNS",
        cex.names=1)#, xlim=c(0,0.5) )
barplot((aa.atomname.fractotal[sc.atomnames]) ,
        main="Sidechain Fraction buried uns / total uns" ,
        xlab="",#BB O Fraction Buried-UNS / Total Buried UNS",
        cex.names=1)#, xlim=c(0,0.5) )

dev.print(device=postscript,file='../fig/UNS_frac_total_uns.ps',horizontal=F)





# hbond energy histograms for backbone O
resnames <- names(rescount)

par(mfrow=c(5,4))
par(mar=c(3,4,3,1))
bbo.minE <- quantile(pdb.sasa0$hbond[pdb.sasa0$atomtype==20],0.001)
for( resname in resnames ) {
  atomname <- paste(resname,"O")
  I <- aa.atomname == atomname
  tmp <- pdb.sasa0$hbond[I]
  tmp[tmp<bbo.minE] <- bbo.minE
  hist(tmp,breaks=seq(from=bbo.minE-0.01,to=0.01,length=15),
       main=paste(atomname),xlab="HBondE")
}
dev.print(device=postscript,file='../fig/hbondE_bbo_hist.ps',horizontal=F)

par(mfrow=c(5,4))
par(mar=c(3,4,3,1))
bbh.minE <- quantile(pdb.sasa0$hbond[pdb.sasa0$atomtype==25],0.001)
for( resname in resnames ) {
  if(resname=="PRO") {
    plot(1,1,type='n',main="PRO H")
  } else {
    atomname <- paste(resname,"H")
    I <- aa.atomname == atomname
    tmp <- pdb.sasa0$hbond[I]
    tmp[tmp<bbh.minE] <- bbh.minE
    hist(tmp,#breaks=seq(from=bbh.minE-0.01,to=0.01,length=15),
         main=paste(atomname),xlab="HBondE")
  }
}
dev.print(device=postscript,file='../fig/hbondE_bbh_hist.ps',horizontal=F)

par(mfrow=c(6,5))
par(mar=c(3,4,3,1))
I <- pdb.sasa0$atomtype!=20 & pdb.sasa0$atomtype!=25
sc.minE <- quantile(pdb.sasa0$hbond[I],0.001)
for( atomname in sc.atomnames ) {
  I <- aa.atomname == atomname
  tmp <- pdb.sasa0$hbond[I]
  tmp[tmp<bbo.minE] <- bbo.minE
  hist(tmp,breaks=seq(from=bbo.minE-0.01,to=0.01,length=15),
       main=paste(atomname),xlab="HBondE")
}
dev.print(device=postscript,file='../fig/hbondE_sidechain_hist.ps',horizontal=F)

# hbond energy vs. residue burial

par(mfrow=c(5,4))
par(mar=c(3,4,3,1))
par(mgp=c(2,1,0))
for( resname in resnames ) {
  atomname <- paste(resname,"O")
  #I <- aa.atomname == atomname & runif(length(aa.atomname)) < 0.2
  I <- aa.atomname == atomname #& runif(length(aa.atomname)) < 0.2
  I <- I & runif(length(I)) < 4000/sum(I)
  print(paste(atomname,sum(I)))
  plot(pdb.sasa0$nb[I]+runif(sum(I)),pdb.sasa0$hbond[I],pch='.',ylim=c(-1.5,0),
       main=paste(atomname,"HB_E"),xlab="neighbors",ylab="hbond E")
}
dev.print(device=postscript,file='../fig/hbondE_bbo_energy_vs_burial.ps',horizontal=F)

par(mfrow=c(5,4))
par(mar=c(3,4,3,1))
par(mgp=c(2,1,0))
for( resname in resnames ) {
  if(resname=="PRO") {
    plot(1,1,type='n',main="PRO H HB_E",xlab="neighbors",ylab="hbond E")
  }else{
    atomname <- paste(resname,"H")
    #I <- aa.atomname == atomname & runif(length(aa.atomname)) < 0.2
    I <- aa.atomname == atomname #& runif(length(aa.atomname)) < 0.2
    I <- I & runif(length(I)) < 4000/sum(I)
    plot(pdb.sasa0$nb[I]+runif(sum(I)),pdb.sasa0$hbond[I],pch='.',ylim=c(-1.5,0),
         main=paste(atomname,"HB_E"),xlab="neighbors",ylab="hbond E")
  }
}
dev.print(device=postscript,file='../fig/hbondE_bbh_energy_vs_burial.ps',horizontal=F)


par(mfrow=c(6,5))
par(mar=c(3,4,3,1))
par(mgp=c(2,1,0))
#I <- pdb.sasa0$atomtype!=20 & pdb.sasa0$atomtype!=25
sc.minE <- quantile(pdb.sasa0$hbond[I],0.001)
for( atomname in sc.atomnames ) {
  I <- aa.atomname == atomname #& runif(length(aa.atomname)) < 0.2
  I <- I & runif(length(I)) < 2000/sum(I)
  plot(pdb.sasa0$nb[I]+runif(sum(I)),pdb.sasa0$hbond[I],pch='.',#ylim=c(-1.5,0),
       main=paste(atomname,"HB_E"),xlab="neighbors",ylab="hbond E")

}
dev.print(device=postscript,file='../fig/hbondE_sidechain_energy_vs_burial.ps',horizontal=F)

#######################################
# buried counts, uns counts, uns frac for vs #nb

par(mfrow=c(5,4))
par(mar=c(3,4,3,1))
par(mgp=c(2,1,0))
for( resname in resnames ) {
  atomname <- paste(resname,"O")
  I <- aa.atomname == atomname
  nb.count <- tapply(pdb.sasa0$hbond[I],pdb.sasa0$nb[I],function(x) length(x))
  nb.unscount <- tapply(pdb.sasa0$hbond[I],pdb.sasa0$nb[I],function(x) sum(x >= -0.01))
  nb.frac <- nb.unscount / nb.count
  nb.frac[nb.count < 50] <- 0
#  plot(as.numeric(names(nb.frac)),nb.frac,
#       main=paste(atomname,"UNS v NB"),xlab='neighbors',ylab='UNS frac')
  plot(as.numeric(names(nb.frac)),nb.count/sum(nb.count),type='l',
       main=paste(atomname,"# v NB"),xlab='neighbors',ylab='Counts')
  lines(as.numeric(names(nb.frac)),nb.unscount/sum(nb.unscount),col=2)
}
#dev.print(device=postscript,file='../fig/hbondE_bbo_unsfrac_vs_burial.ps',horizontal=F)
dev.print(device=postscript,file='../fig/UNS_bbo_counts_vs_burial.ps',horizontal=F)

par(mfrow=c(5,4))
par(mar=c(3,4,3,1))
par(mgp=c(2,1,0))
for( resname in resnames ) {
  if(resname=="PRO") {
    plot(1,1,type='n',main="PRO H UNS v NB",xlab="neighbors",ylab="hbond E")
  }else{
    atomname <- paste(resname,"H")
    I <- aa.atomname == atomname
    nb.count <- tapply(pdb.sasa0$hbond[I],pdb.sasa0$nb[I],function(x) length(x))
    nb.unscount <- tapply(pdb.sasa0$hbond[I],pdb.sasa0$nb[I],function(x) sum(x >= -0.01))
    nb.frac <- nb.unscount / nb.count
    nb.frac[nb.count < 50] <- 0
    #plot(as.numeric(names(nb.frac)),nb.frac,
    #     main=paste(atomname,"UNS v NB"),xlab='neighbors',ylab='UNS frac')
    plot(as.numeric(names(nb.frac)),nb.count/sum(nb.count),type='l',
         main=paste(atomname,"# v NB"),xlab='neighbors',ylab='Counts')
    lines(as.numeric(names(nb.frac)),nb.unscount/sum(nb.unscount),col=2)
  }
}
#dev.print(device=postscript,file='../fig/hbondE_bbh_unsfrac_vs_burial.ps',horizontal=F)
dev.print(device=postscript,file='../fig/UNS_bbh_counts_vs_burial.ps',horizontal=F)

par(mfrow=c(6,5))
par(mar=c(3,4,3,1))
par(mgp=c(2,1,0))
for( atomname in sc.names ) {
  #atomname <- paste(resname,"O")
  I <- aa.atomname == atomname
  nb.count <- tapply(pdb.sasa0$hbond[I],pdb.sasa0$nb[I],function(x) length(x))
  nb.unscount <- tapply(pdb.sasa0$hbond[I],pdb.sasa0$nb[I],function(x) sum(x >= -0.01))
  nb.frac <- nb.unscount / nb.count
  nb.frac[nb.count < 50] <- 0
#  plot(as.numeric(names(nb.frac)),nb.frac,
#       main=paste(atomname,"UNS v NB"),xlab='neighbors',ylab='UNS frac')
  plot(as.numeric(names(nb.frac)),nb.count/sum(nb.count),type='l',
       main=paste(atomname,"# v NB"),xlab='neighbors',ylab='Counts')
  lines(as.numeric(names(nb.frac)),nb.unscount/sum(nb.unscount),col=2)
}
#dev.print(device=postscript,file='../fig/UNS_sidechain_unsfrac_vs_burial.ps',horizontal=F)
dev.print(device=postscript,file='../fig/UNS_sidechain_counts_vs_burial.ps',horizontal=F)

# overall counts for all atoms vs. NB

par(mfrow=c(2,1))
par(mar=c(3,4,3,1))
par(mgp=c(2,1,0))
nb.count <- tapply(pdb.sasa0$hbond,pdb.sasa0$nb,function(x) length(x))
nb.unscount <- tapply(pdb.sasa0$hbond,pdb.sasa0$nb,function(x) sum(x >= -0.01))
nb.frac <- nb.unscount / nb.count
nb.frac[nb.count < 50] <- 0

plot(as.numeric(names(nb.frac)),nb.frac,type='l',
     main="All atoms UNS & burried / burried vs. NB",
     xlab='neighbors',ylab='UNS frac')
plot(as.numeric(names(nb.frac)),nb.count/sum(nb.count),type='l',
     main="All atoms counts v NB (buried=black, burieduns=red)",
     xlab='neighbors',ylab='Counts')
lines(as.numeric(names(nb.frac)),nb.unscount/sum(nb.unscount),col=2)
dev.print(device=postscript,file='../fig/UNS_total_unsfrac_counts.ps',horizontal=F)


###### # of cases vs. burial


##############################################
# stats for natives and decoys
##############################################

load('kirahbond.rdata')


# all residue stats (not just buried)
rescount <- tapply(pdbdata$aa,pdbdata$aa,function(x) length(x))
resmeannb    <- tapply(pdbdata$nb,pdbdata$aa,mean)
resmeansasafrac <- tapply(pdbdata$sasafrac,pdbdata$aa,mean)

# backbone uns stats
Ibbh.kira <- kirahbond$atomtype == 25
Ibbo.kira <- kirahbond$atomtype == 20
kirahbond.bbh <- kirahbond[Ibbh.kira,]
kirahbond.bbo <- kirahbond[Ibbo.kira,]

bbh.count.kira    <- tapply(kirahbond.bbh$hbond,kirahbond.bbh$aa,function(x) length(x) )
bbh.unscount.kira <- tapply(kirahbond.bbh$hbond,kirahbond.bbh$aa,function(x) sum(x >= -0.01) )
bbh.unsfrac.kira  <- bbh.unscount.kira / bbh.count.kira
bbh.meannb.kira   <- tapply(kirahbond.bbh$nb,kirahbond.bbh$aa,mean)

bbo.count.kira    <- tapply(kirahbond.bbo$hbond,kirahbond.bbo$aa,function(x) length(x) )
bbo.unscount.kira <- tapply(kirahbond.bbo$hbond,kirahbond.bbo$aa,function(x) sum(x >= -0.01) )
bbo.unsfrac.kira  <- bbo.unscount.kira / bbo.count.kira
bbo.meannb.kira   <- tapply(kirahbond.bbo$nb,kirahbond.bbo$aa,mean)

par(mfrow=c(2,1))
par(las=2)
barplot(rbind(bbh.unsfrac,bbh.unsfrac.kira),beside=T,
        main="backbone H UNS frac",,col=c(1,2),space=c(0.2,1) )
barplot(rbind(bbo.unsfrac,bbo.unsfrac.kira),beside=T,
        main="backbone O UNS frac",col=c(1,2),space=c(0.2,1) )

dev.print(device=postscript,file='../fig/uns_bb_decoys.ps',horizontal=F)




aa.atomname.kira <- as.factor(paste(as.character(kirahbond$aa),
                               as.character(kirahbond$atomname)))
aa.atomtype.kira <- as.factor(paste(as.character(kirahbond$aa),
                               as.character(kirahbond$atomtype)))
aa.atomname.count.kira <- tapply(aa.atomname.kira,aa.atomname.kira,function(x) length(x))
aa.atomtype.count.kira <- tapply(aa.atomtype.kira,aa.atomtype.kira,function(x) length(x))
aa.atomname.unscount.kira <- tapply(kirahbond$hbond,
                                    aa.atomname.kira,function(x) sum(x >= -0.01))
aa.atomtype.unscount.kira <- tapply(kirahbond$hbond,
                                    aa.atomtype.kira,function(x) sum(x >= -0.01))
aa.atomname.unsfrac.kira <- aa.atomname.unscount.kira / aa.atomname.count.kira
aa.atomtype.unsfrac.kira <- aa.atomtype.unscount.kira / aa.atomtype.count.kira

par(mfrow=c(1,1))
par(las=2)
par(mar=c(4,6,4,2))
barplot(rbind(aa.atomname.unsfrac.kira[rev(sc.names)],aa.atomname.unsfrac[rev(sc.names)]),
        beside=T,horiz=T,space=c(0,1),col=c(2,1),
        main="UNS frac for Side Chain Atoms black=native red=kira")
dev.print(device=postscript,file='../fig/uns_sc_decoys.ps',horizontal=F)



sc.names.by.type <- c("ARG 1HH1","ARG 1HH2","ARG 2HH1","ARG 2HH2","ARG HE",#NA,
                      "ASN 1HD2","ASN 2HD2",#NA,
                      "CYS HG",#NA,
                      "GLN 1HE2","GLN 2HE2",#NA,
                      "HIS HD1", "HIS HE2",#NA,
                      "LYS 1HZ",  "LYS 2HZ",  "LYS 3HZ",#NA,
                      "SER HG",#NA,
                      "THR HG1",#NA,
                      "TRP HE1",#NA,
                      "TYR HH", #NA, # 22
                      NA,#NA,
                      "SER OG", # NA,
                      "THR OG1",# NA,
                      "TYR OH", # 13
                      NA,#NA,
                      "ASP OD1",  "ASP OD2",#NA,
                      "GLU OE1",  "GLU OE2", # 14
                      NA,#NA,
                      "ASN OD1",#NA,
                      "GLN OE1", # 15
                      NA, #NA,
                      "HIS ND1",  "HIS NE2",  ) # 8

par(mfrow=c(1,1))
par(las=2)
par(mar=c(4,6,4,2))
barplot(rbind(aa.atomname.unsfrac.kira[rev(sc.names.by.type)],
              aa.atomname.unsfrac[rev(sc.names.by.type)]),
        beside=T,horiz=T,space=c(0,1),col=c(2,1),
        main="UNS frac for Side Chain Atoms black=native red=kira")
dev.print(device=postscript,file='../fig/uns_sc_decoys_by_type.ps',horizontal=F)


##################################################
##### with kira natives


# backbone uns stats
load('kiranathbond.rdata')

Ibbh.kiranat <- kiranathbond$atomtype == 25
Ibbo.kiranat <- kiranathbond$atomtype == 20
kiranathbond.bbh <- kiranathbond[Ibbh.kiranat,]
kiranathbond.bbo <- kiranathbond[Ibbo.kiranat,]

bbh.count.kiranat    <- tapply(kiranathbond.bbh$hbond,kiranathbond.bbh$aa,
                               function(x) length(x) )
bbh.unscount.kiranat <- tapply(kiranathbond.bbh$hbond,kiranathbond.bbh$aa,
                               function(x) sum(x > -0.01) )
bbh.unsfrac.kiranat  <- bbh.unscount.kiranat / bbh.count.kiranat
bbh.meannb.kiranat   <- tapply(kiranathbond.bbh$nb,kiranathbond.bbh$aa,mean)

bbo.count.kiranat    <- tapply(kiranathbond.bbo$hbond,kiranathbond.bbo$aa,
                               function(x) length(x) )
bbo.unscount.kiranat <- tapply(kiranathbond.bbo$hbond,kiranathbond.bbo$aa,
                               function(x) sum(x > -0.01) )
bbo.unsfrac.kiranat  <- bbo.unscount.kiranat / bbo.count.kiranat
bbo.meannb.kiranat   <- tapply(kiranathbond.bbo$nb,kiranathbond.bbo$aa,mean)

par(mfrow=c(2,1))
par(las=2)
barplot(rbind(bbh.unsfrac,bbh.unsfrac.kira,bbh.unsfrac.kiranat),beside=T,
        main="backbone H UNS frac",,col=c(1,2,3),space=c(0.2,1) )
barplot(rbind(bbo.unsfrac,bbo.unsfrac.kira,bbo.unsfrac.kiranat),beside=T,
        main="backbone O UNS frac",col=c(1,2,3),space=c(0.2,1) )

dev.print(device=postscript,file='../fig/uns_bb_decoys_with_native.ps',horizontal=F)


aa.atomname.kiranat <- as.factor(paste(as.character(kiranathbond$aa),
                               as.character(kiranathbond$atomname)))
aa.atomtype.kiranat <- as.factor(paste(as.character(kiranathbond$aa),
                               as.character(kiranathbond$atomtype)))
aa.atomname.count.kiranat <- tapply(aa.atomname.kiranat,aa.atomname.kiranat,
                                    function(x) length(x))
aa.atomtype.count.kiranat <- tapply(aa.atomtype.kiranat,aa.atomtype.kiranat,
                                    function(x) length(x))
aa.atomname.unscount.kiranat <- tapply(kiranathbond$hbond,
                                    aa.atomname.kiranat,function(x) sum(x >= -0.01))
aa.atomtype.unscount.kiranat <- tapply(kiranathbond$hbond,
                                    aa.atomtype.kiranat,function(x) sum(x >= -0.01))
aa.atomname.unsfrac.kiranat <- aa.atomname.unscount.kiranat / aa.atomname.count.kiranat
aa.atomtype.unsfrac.kiranat <- aa.atomtype.unscount.kiranat / aa.atomtype.count.kiranat

par(mfrow=c(1,1))
par(las=2)
par(mar=c(4,6,4,2))
barplot(rbind(aa.atomname.unsfrac.kiranat[rev(sc.names)],
              aa.atomname.unsfrac.kira[rev(sc.names)],
              aa.atomname.unsfrac[rev(sc.names)]),
        beside=T,horiz=T,space=c(0,1),col=c(3,2,1),
        main="UNS frac for Side Chain Atoms black=native red=kira")
dev.print(device=postscript,file='../fig/uns_sc_decoys_with_native.ps',horizontal=F)


# buried UNS independence

# why are OH groups likely to be buried uns?
# could be because hard to satisy both groups at same time?
count.sat.oh.groups <- function() {

  # both OH buried
  I <- pdb.sasa0$aa=="TYR" & 
       (pdb.sasa0$atomname=="OH" | pdb.sasa0$atomname=="HH") # and buried
  pdbres <- as.factor(paste(as.character(pdb.sasa0$pdb[I]),pdb.sasa0$res[I]))
  tyr.ohhh.buried <- length(pdbres) - length(levels(pdbres))
  print(tyr.ohhh.buried)
  # both OH buried & OH SAT
  I <- pdb.sasa0$aa=="TYR" & 
       ((pdb.sasa0$atomname=="OH"&pdb.sasa0$hbond < -0.01) | pdb.sasa0$atomname=="HH")
  pdbres <- as.factor(paste(as.character(pdb.sasa0$pdb[I]),pdb.sasa0$res[I]))
  tyr.oh.sat <- length(pdbres) - length(levels(pdbres))
  print(tyr.oh.sat)
  # both OH buried & OH SAT
  I <- pdb.sasa0$aa=="TYR" & 
       ((pdb.sasa0$atomname=="HH"&pdb.sasa0$hbond < -0.01) | pdb.sasa0$atomname=="OH")
  pdbres <- as.factor(paste(as.character(pdb.sasa0$pdb[I]),pdb.sasa0$res[I]))
  tyr.hh.sat <- length(pdbres) - length(levels(pdbres))
  print(tyr.hh.sat)  
  # both OH buried and both SAT
  I <- pdb.sasa0$aa=="TYR" & pdb.sasa0$hbond < -0.01 &
       (pdb.sasa0$atomname=="OH" | pdb.sasa0$atomname=="HH")
  pdbres <- as.factor(paste(as.character(pdb.sasa0$pdb[I]),pdb.sasa0$res[I]))
  tyr.ohhh.sat <- length(pdbres) - length(levels(pdbres))
  print(tyr.ohhh.sat)

  Pohsat <- tyr.oh.sat / tyr.ohhh.buried
  Phhsat <- tyr.hh.sat / tyr.ohhh.buried
  Pohhhsat <- tyr.ohhh.sat / tyr.ohhh.buried
  Pohhhindep <- Pohsat * Phhsat
  print(paste("indep:",Pohhhindep,"actual:",Pohhhsat))

}


assess.uns.independence <- function() {
  aa <- "PRO"; i <- 2; j <- 1
  iters <- 1
  aalist <- c(); name1list <- c(); name2list <- c()
  indeplist <- c(); actuallist <- c();
  Nbblist <- c(); N1ulist <- c(); N2ulist <- c(); Nbulist <- c(); 
  for(ii in 1:20) {    
    aa <- aanames3[ii]
    print(paste(ii,aa))
    tmphbond <- pdb.sasa0[pdb.sasa0$aa==aa,]
    resatoms <- levels(as.factor(as.character(tmphbond$atomname)))
   if(length(resatoms) > 1) {
    for(i in 2:length(resatoms)) {
      for(j in 1:(i-1)) {        
        name1 <- resatoms[i]
        name2 <- resatoms[j]
        # both buried
        I <- (tmphbond$atomname==name1 | tmphbond$atomname==name2) # and buried
        pdbres <- as.factor(paste(as.character(tmphbond$pdb[I]),tmphbond$res[I]))
        both.buried <- length(pdbres) - length(levels(pdbres))
        # both buried & name1 unsat
        I <- ((tmphbond$atomname==name1 & tmphbond$hbond>=-0.01)
              | tmphbond$atomname==name2) # and buried
        pdbres <- as.factor(paste(as.character(tmphbond$pdb[I]),tmphbond$res[I]))
        name1.unsat <- length(pdbres) - length(levels(pdbres))
        # both buried & name2 unsat
        I <- (tmphbond$atomname==name1 |
              (tmphbond$atomname==name2 & tmphbond$hbond>=-0.01)) # and buried
        pdbres <- as.factor(paste(as.character(tmphbond$pdb[I]),tmphbond$res[I]))
        name2.unsat <- length(pdbres) - length(levels(pdbres))
        # both buried & name2 unsat
        I <- ((tmphbond$atomname==name1 | tmphbond$atomname==name2)
              & tmphbond$hbond>=-0.01 ) # and buried
        pdbres <- as.factor(paste(as.character(tmphbond$pdb[I]),tmphbond$res[I]))
        both.unsat <- length(pdbres) - length(levels(pdbres))

        indep <-  round(name1.unsat/both.buried * name2.unsat/both.buried, 4 )
        actual <- round(both.unsat/both.buried , 4 )
        print(paste(iters,
                    aa,name1,name2,
                    "indep:",  round(indep, 4) ,
                    "actual:", round(actual,4) ,
                    "Nsamp:",
                    both.buried,
                    name1.unsat,
                    name2.unsat,
                    both.unsat))
        iters <- iters + 1
        aalist <- append(aalist,aa)
        name1list <- append(name1list,name1)
        name2list <- append(name2list,name2)

        #if(length(name2list) != iters-1)
        #  print(paste("AAH!",iters))
        
        indeplist <- append(indeplist,indep)
        actuallist <- append(actuallist,actual)
        Nbblist   <- append(Nbblist,both.buried)
        N1ulist   <- append(N1ulist,name1.unsat)
        N2ulist   <- append(N2ulist,name2.unsat)
        Nbulist   <- append(Nbulist,both.unsat)
      }
    }
   }
  }
  lr <- round(log2( actuallist / indeplist ),3)
  uns.indep <- as.data.frame(list(aalist,name1list,name2list,
                                  indeplist,actuallist, lr ,
                                  Nbblist,N1ulist,N2ulist,Nbulist))
  names(uns.indep) <- c("aa","atom1","atom2","Pindept","Pactual","LR",
                          "Nbthbur","Nat1uns","Nat2uns","Nbthuns",)

}
