r <- read.table('pdb_region.data',header=F)
pr <- r[,-1]
names(pr) <- c('prot','pdb','chain','res','aa','ss','nb',
               'atr1','res1','spk1','bup1',
               'r5', 'nb5', 'atr5', 'res5', 'spk5', 'bup5', 'rms5',
               'r7', 'nb7', 'atr7', 'res7', 'spk7', 'bup7', 'rms7',
               'r10','nb10','atr10','res10','spk10','bup10','rms10',
               'r15','nb15','atr15','res15','spk15','bup15','rms15'      )

save(pr,file="~/project/features/pr.pdb.region.rdata")
