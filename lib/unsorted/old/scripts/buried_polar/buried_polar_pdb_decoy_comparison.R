#################################################################3
#        with SS and more decoys
#################################################################

setwd('/users/sheffler/project/features/buried_polar/')


load('rdata/dp.decoyset2.polar.rdata')
load('rdata/pp.pdb.polar.rdata')

atomtype <- dp$type
type <- c()
type[dp$type== 8] <- "NH"
type[dp$type==13] <- "OH"
type[dp$type==14] <- "OA"
type[dp$type==15] <- "OC"
type[dp$type==20] <- "OB"
type[dp$type==22] <- "HP"
type[dp$type==25] <- "HB"
isbb <- rep(F,length(type))
isbb[dp$type==20] <- T
isbb[dp$type==25] <- T
type <- paste(as.character(dp$aa),as.character(type),sep='.')
dp$type <- as.factor(type)
dp <- cbind(dp,isbb)
dp <- cbind(dp,atomtype)
save(dp,file='rdata/dp.decoyset2.polar.rdata')

atomtype <- pp$type
type <- c()
type[pp$type== 8] <- "NH"
type[pp$type==13] <- "OH"
type[pp$type==14] <- "OA"
type[pp$type==15] <- "OC"
type[pp$type==20] <- "OB"
type[pp$type==22] <- "HP"
type[pp$type==25] <- "HB"
isbb <- rep(F,length(type))
isbb[pp$type==20] <- T
isbb[pp$type==25] <- T
type <- paste(as.character(pp$aa),as.character(type),sep='.')
pp$type <- as.factor(type)
pp <- cbind(pp,isbb)
pp <- cbind(pp,atomtype)
save(pp,file='rdata/pp.pdb.polar.rdata')

aanames3 <- c("ALA","ARG","ASN","ASP","CYS","GLN","GLU","GLY","HIS","ILE",
              "LEU","LYS","MET","PHE","PRO","SER","THR","TRP","TYR","VAL")

otypes.nosp <- unique(pp$type[!pp$isbb])[
                order(substr(unique(as.character(pp$type[!pp$isbb])),5,6))]

otypessc <- rev(c("ARG.HP","ASN.HP","GLN.HP","THR.HP","SER.HP","LYS.HP",
                  "TYR.HP","TRP.HP","HIS.HP","CYS.HP",NA,"HIS.NH",NA,"ASN.OA",
                  "GLN.OA",NA,"ASP.OC","GLU.OC",NA,"SER.OH","THR.OH","TYR.OH"))
otypesbbh <- rev(c(paste(aanames3[-15],'HB',sep='.')))
otypesbbo <- rev(c(paste(aanames3,'OB',sep='.')))
otypesbb  <- c(otypesbbh,NA,otypesbbo)
otypes    <- c(otypesbb,NA,otypessc)

isbb  <- substr(levels(pp$type),5,6) %in% c('OB','HB')
isbbh <- substr(levels(pp$type),5,6) %in% c('HB')
isbbo <- substr(levels(pp$type),5,6) %in% c('OB')

########## burial thresholds #############

pp.bt <- tapply(pp$anb, pp$type, median)
dp.bt <- tapply(dp$anb, dp$type, median)


par(cex=0.7)
par(mar=c(5,6,3,1))
par(las=2)

barplot(rbind(pp.bt[otypes],dp.bt[otypes]),
        beside=T,col=c(1,2),space=c(0.2,1),horiz=T,
        main="PDB Median Burial")
dev.print(file='fig/burial_median_pdb.ps',device=postscript,horizontal=F)


########## uns frac stuff ################

pp.ct <- tapply(pp$hbe,   pp$type, length)
pp.uc <- tapply(pp$hbe==0,pp$type, sum)
pp.uf <- pp.uc/pp.ct

dp.ct <- tapply(dp$hbe,   dp$type, length)
dp.uc <- tapply(dp$hbe==0,dp$type, sum)
dp.uf <- dp.uc/dp.ct

par(cex=0.7)
par(mar=c(5,6,3,1))
par(las=2)

barplot(rbind(pp.ct[otypes],pp.uc[otypes]),
        beside=T,col=c(1,2),space=c(0.2,1),horiz=T,xlim=c(0,5000),
        main="PDB Buried Polar Counts")
dev.print(file='fig/type_counts_pdb.ps',device=postscript,horizontal=F)

barplot(rbind(dp.ct[otypes],dp.uc[otypes]),
        beside=T,col=c(1,2),space=c(0.2,1),horiz=T,xlim=c(0,5000),
        main="Decoy Buried Polar Counts")
dev.print(file='fig/type_counts_decoy.ps',device=postscript,horizontal=F)

barplot(rbind(pp.uf[otypes],dp.uf[otypes]),
        beside=T,col=c(1,2),space=c(0.2,1),horiz=T,xlim=c(0,1),
        main="Buried Polar UNS Frac PDB=Black, Decoy=Red")
dev.print(file='fig/uns_frac_pdb_v_decoy.ps',device=postscript,horizontal=F)

barplot(dp.uf[otypes]/pp.uf[otypes],
        beside=T,col=c(1),space=c(0.2,1),horiz=T,
        main="Buried Polar UNS Frac Ratio Decoy/PDB")
lines(c(1,1),c(0,100),col=2)
dev.print(file='fig/uns_frac_ratio_pdb_v_decoy.ps',device=postscript,horizontal=F)



pp.ct.ss <- tapply(pp$hbe,    list(pp$type,pp$ss), length)
pp.uc.ss <- tapply(pp$hbe==0, list(pp$type,pp$ss), sum)
pp.uf.ss <- pp.uc.ss / pp.ct.ss

dp.ct.ss <- tapply(dp$hbe,    list(dp$type,dp$ss), length)
dp.uc.ss <- tapply(dp$hbe==0, list(dp$type,dp$ss), sum)
dp.uf.ss <- dp.uc.ss / dp.ct.ss

ur.ss <- dp.uf.ss / pp.uf.ss


par(cex=0.7)
par(mar=c(5,6,3,1))
par(las=2)

barplot(rbind(pp.uc.ss[,'L'][otypes],pp.uc.ss[,'E'][otypes],pp.uc.ss[,'H'][otypes]),
        beside=T,col=3:1,space=c(0.2,1),horiz=T,xlim=c(0,1000),
        main="PDB Buried Polar UNS Counts by SS")
dev.print(file='fig/type_counts_pdb_ss.ps',device=postscript,horizontal=F)

barplot(rbind(dp.uc.ss[,'L'][otypes],dp.uc.ss[,'E'][otypes],dp.uc.ss[,'H'][otypes]),
        beside=T,col=3:1,space=c(0.2,1),horiz=T,xlim=c(0,1000),
        main="Decoy Buried Polar UNS Counts by SS")
dev.print(file='fig/type_counts_decoy_ss.ps',device=postscript,horizontal=F)

barplot(rbind(pp.uf.ss[,'L'][otypes],pp.uf.ss[,'E'][otypes],pp.uf.ss[,'H'][otypes]),
        beside=T,col=3:1,space=c(0.2,1),horiz=T,
        main="PDB Buried Polar UNS Frac by SS")
dev.print(file='fig/uns_frac_decoy_ss.ps',device=postscript,horizontal=F)

barplot(rbind(dp.uf.ss[,'L'][otypes],dp.uf.ss[,'E'][otypes],dp.uf.ss[,'H'][otypes]),
        beside=T,col=3:1,space=c(0.2,1),horiz=T,
        main="Decoy Buried Polar UNS Frac by SS")
dev.print(file='fig/uns_frac_decoy_ss.ps',device=postscript,horizontal=F)

barplot(rbind(ur.ss[,'L'][otypes],
              ur.ss[,'E'][otypes],
              ur.ss[,'H'][otypes]),
        beside=T,col=3:1,space=c(0.2,1),horiz=T,
        main="Decoy Buried Polar UNS Frac Ratio Decoy/PDB by SS")
dev.print(file='fig/uns_frac_ratio_pdb_v_decoy_ss.ps',device=postscript,horizontal=F)





bth <- 300

pp.br <- c();
pp.br[pp$anb >bth] <- "buried";
pp.br[pp$anb<=bth] <- "surface"; 
pp.br <- as.factor(pp.br)
dp.br <- c();
dp.br[dp$anb >bth] <- "buried";
dp.br[dp$anb<=bth] <- "surface"; 
dp.br <- as.factor(dp.br)


pp.ct.br <- tapply(pp$hbe,    list(pp$type,pp.br), length)
pp.uc.br <- tapply(pp$hbe==0, list(pp$type,pp.br), sum)
pp.uf.br <- pp.uc.br / pp.ct.br

dp.ct.br <- tapply(dp$hbe,    list(dp$type,dp.br), length)
dp.uc.br <- tapply(dp$hbe==0, list(dp$type,dp.br), sum)
dp.uf.br <- dp.uc.br / dp.ct.br

ur.br <- dp.uf.br / pp.uf.br


par(cex=0.7)
par(mar=c(5,6,3,1))
par(las=2)

barplot(rbind(pp.uc.br[,'buried'][otypes],pp.uc.br[,'surface'][otypes]),
        beside=T,col=1:2,space=c(0.2,1),horiz=T,xlim=c(0,1000),
        main="PDB Buried Polar UNS Counts by Burial")
dev.print(file='fig/type_counts_pdb_burial.ps',device=postscript,horizontal=F)

barplot(rbind(dp.uc.br[,'buried'][otypes],dp.uc.br[,'surface'][otypes]),
        beside=T,col=1:2,space=c(0.2,1),horiz=T,xlim=c(0,1000),
        main="Decoy Buried Polar UNS Counts by Burial")
dev.print(file='fig/type_counts_decoy_burial.ps',device=postscript,horizontal=F)

barplot(rbind(pp.uf.br[,'buried'][otypes],pp.uf.br[,'surface'][otypes]),
        beside=T,col=1:2,space=c(0.2,1),horiz=T,
        main="PDB Buried Polar UNS Frac by Burial")
dev.print(file='fig/uns_frac_decoy_burial.ps',device=postscript,horizontal=F)

barplot(rbind(dp.uf.br[,'buried'][otypes],dp.uf.br[,'surface'][otypes]),
        beside=T,col=1:2,space=c(0.2,1),horiz=T,
        main="Decoy Buried Polar UNS Frac by Burial")
dev.print(file='fig/uns_frac_decoy_burial.ps',device=postscript,horizontal=F)

barplot(rbind(ur.br[,'buried'][otypes],ur.br[,'surface'][otypes]),
        beside=T,col=1:2,space=c(0.2,1),horiz=T,
        main="Decoy Buried Polar UNS Frac Ratio Decoy/PDB by Burial")
dev.print(file='fig/uns_frac_ratio_pdb_v_decoy_burial.ps',device=postscript,horizontal=F)









