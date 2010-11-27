load('pdb.sasa0.rdata')
load('pdbdata.rdata')
load('pdbresmap.rdata')
load('pdbinfo.rdata')


pdbdir <- "/data/sheffler/pdb/rosetta_pdb_chains/"
prot  <- "1c75A"
file  <- "~/project/packing_index/visualize/uns_hbond/test.rasmol"

rasmol.whole.chain <- function(prot,file,pdbdir,I) {
  rosres <- pdb.sasa0$res[I]  
  #pdbres <- pdbresmap$pdbresnum[pdbresmap$pdb==prot][rosres]
  pdbatom <- pdb.sasa0$atomname[I]

  # write rasmol script
  write('zap',file=file)
  write("stereo -5",file=file,append=T)
  write("set picking center",file=file,append=T)
  write(paste('load pdb "',pdbdir,prot,'.pdb"',sep=""),file=file,append=T)
  write("restrict none",file=file,append=T)
  write("select backbone",file=file,append=T)
  write("strands",file=file,append=T)
  write("select not hydrogen",file=file,append=T)
  write("wireframe",file=file,append=T)
  for(ii in 1:length(rosres)) {
    sel <- paste(paste("(",rosres[ii],' and ','*.',pdbatom[ii],')',sep=''))
    write(paste('select',sel),file=file,append=T)
    write("spacefill",file=file,append=T)
    write(paste("select within(4.0,(",sel,"))"),file=file,append=T)
    write("wireframe 0.2",file=file,append=T)
    write(paste("select within(3.0,(",sel,")) and not",rosres[ii]),file=file,append=T)
    write("wireframe 0.4",file=file,append=T)
  }
  
}

find.pdb.examples <- function() {
  hirespdb <- c("1gvkB", "1luqA", "1rtqA", "1i1wA", "1g66A", "1k5cA",
                "1n55A", "1yljA", "1muwA", "1n4wA", "1f9yA", "1c75A",
                "1sfdA", "1et1A", "1ok0A", "1x6zA", "1m1qA", "1tqgA",
                "1gweA", "1pjxA", "1l9lA", "1nkiA", "1v0lA", "1vbwA",
                "1g6xA", "1w0nA", "1ug6A", "1wy3A", "1p1xA", "1ufyA",
                "1kwfA", "1p9gA", "1xg0A", "1xg0C", "1ix9A", "1dy5A",
                "1unqA", "1ucsA", "1byi_", "3lzt_", "7a3hA", "1ejgA",
                "1vyrA", "1hj9A", "1mwqA", "1ixh_", "2fdn_", "1kthA",
                "1us0A", "1x8qA", "1nls_", "1g4iA", "1ic6A", "1gqvA",
                "1nwzA", "1iuaA", "1m40A", "1r6jA", "1j0pA", "1gci_",
                "2bf9A", "1f94A", "1mj5A", "1oewA", "1iqzA", "1k4iA",
                "1ssxA", "1mc2A", "1rb9_", "1pq7A", "1aho_", "1lugA",
                "2pvbA", "1ea7A", "1v6pA" )  
  pdbdir <- "/data/sheffler/pdb/rosetta_pdb_chains/"
  filepath  <- "~/project/packing_index/visualize/uns_hbond/whole_chain/"
  for(ii in 1:length(hirespdb)){
    print(ii)
    prot <- hirespdb[ii]
    file <- paste(filepath,prot,'_uns.rasmol',sep='')
    I <- pdb.sasa0$pdb==prot & pdb.sasa0$hbond > -0.01 #&TODO & buried 
    rasmol.whole.chain(prot,file,pdbdir,I)
  }
}


pdbdir <- "/data/sheffler/pdb/rosetta_pdb_chains/"
prot  <- "1w2wA"
file  <- "~/project/packing_index/visualize/uns_hbond/testatom.rasmol"
res   <- 59
atom  <- "HG"

rasmol.one.atom <- function(prot,res,atom,file,pdbdir) {
  #i <- which(pdb.sasa0$pdb==prot & pdb.sasa0$res==res & pdb.sasa0$atomname==atom)
  write('zap',file=file)
  write("stereo -5",file=file,append=T)
  write("set picking distance",file=file,append=T)
  write(paste('load pdb "',pdbdir,prot,'.pdb"',sep=""),file=file,append=T)
  write("restrict none",file=file,append=T)
  sel <- paste(paste("(",res,' and ','*.',atom,')',sep=''))
  write(paste('select',sel),file=file,append=T)
  write("spacefill",file=file,append=T)
  write("center selected",file=file,append=T)
  write(paste("select within(6.0,(",sel,"))"),file=file,append=T)
  write("wireframe 0.2",file=file,append=T)
  write("spacefill 0.3",file=file,append=T)
  write(paste("select within(3.0,(",sel,")) and not",rosres[ii]),file=file,append=T)
  write("wireframe 0.4",file=file,append=T)
  write("spacefill 0.5",file=file,append=T)
  write("select all",file=file,append=T)
  write("strands",file=file,append=T)
  write("zoom 300",file=file,append=T)
  
}


find.cys.examples <- function() {
  I <- pdb.sasa0$aa=="CYS" & pdb.sasa0$hbond >= -0.01 & pdb.sasa0$atomname=="HG"
  pdb <- pdb.sasa0$pdb[I]
  res <- pdb.sasa0$res[I]
  pdbdir <- "/data/sheffler/pdb/rosetta_pdb_chains/"
  filepath  <- "~/project/packing_index/visualize/uns_hbond/cys/"
  perm <- order(runif(length(pdb)))
  for(ii in 1:50){#length(pdb)) {
    print(ii)
    idx <- perm[ii]
    file <- paste(filepath,pdb[idx],'_cys_',res[idx],'_',idx,'.rasmol',sep="")
    rasmol.one.atom(pdb[idx],res[idx],"HG",file,pdbdir)
  }
}

find.tyr.examples <- function() {
  I <- pdb.sasa0$aa=="TYR" & pdb.sasa0$hbond >= -0.01 & pdb.sasa0$atomname=="OH"
  pdb <- pdb.sasa0$pdb[I]
  res <- pdb.sasa0$res[I]
  pdbdir <- "/data/sheffler/pdb/rosetta_pdb_chains/"
  filepath  <- "~/project/packing_index/visualize/uns_hbond/tyrOH/"
  perm <- order(runif(length(pdb)))
  for(ii in 1:50){#length(pdb)) {
    print(ii)
    idx <- perm[ii]
    file <- paste(filepath,pdb[idx],'_tyrOH_',res[idx],'_',idx,'.rasmol',sep="")
    rasmol.one.atom(pdb[idx],res[idx],"OH",file,pdbdir)
  }
}


find.buried.arg.bbo.examples <- function() {
  I <- pdb.sasa0$aa=="ARG" & pdb.sasa0$hbond >= -0.01 &
       pdb.sasa0$atomname=="O" & pdb.sasa0$nb == 30
  pdb <- pdb.sasa0$pdb[I]
  res <- pdb.sasa0$res[I]
  pdbdir <- "/data/sheffler/pdb/rosetta_pdb_chains/"
  filepath  <- "~/project/packing_index/visualize/uns_hbond/buried.arg.bbo/"
    for(ii in 1:length(pdb)){#length(pdb)) {
    print(ii)
    idx <- ii
    file <- paste(filepath,pdb[idx],'_arg_bbo_',res[idx],'_',idx,'.rasmol',sep="")
    rasmol.one.atom(pdb[idx],res[idx],"O",file,pdbdir)
  }

}

find.nb.burial.examples <- function() {
  hirespdb <- c("1gvkB", "1luqA", "1rtqA", "1i1wA", "1g66A", "1k5cA",
                "1n55A", "1yljA", "1muwA", "1n4wA", "1f9yA", "1c75A",
                "1sfdA", "1et1A", "1ok0A", "1x6zA", "1m1qA", "1tqgA",
                "1gweA", "1pjxA", "1l9lA", "1nkiA", "1v0lA", "1vbwA",
                "1g6xA", "1w0nA", "1ug6A", "1wy3A", "1p1xA", "1ufyA",
                "1kwfA", "1p9gA", "1xg0A", "1xg0C", "1ix9A", "1dy5A",
                "1unqA", "1ucsA", "1byi_", "3lzt_", "7a3hA", "1ejgA",
                "1vyrA", "1hj9A", "1mwqA", "1ixh_", "2fdn_", "1kthA",
                "1us0A", "1x8qA", "1nls_", "1g4iA", "1ic6A", "1gqvA",
                "1nwzA", "1iuaA", "1m40A", "1r6jA", "1j0pA", "1gci_",
                "2bf9A", "1f94A", "1mj5A", "1oewA", "1iqzA", "1k4iA",
                "1ssxA", "1mc2A", "1rb9_", "1pq7A", "1aho_", "1lugA",
                "2pvbA", "1ea7A", "1v6pA" )  
  pdbdir <- "/data/sheffler/pdb/rosetta_pdb_chains/"
  filepath  <- "~/project/packing_index/visualize/uns_hbond/nb_burial/"
  for(ii in 1:length(hirespdb)){
    print(ii)
    prot <- hirespdb[ii]
    file <- paste(filepath,prot,'_nb_burial.rasmol',sep='')
    I <- pdb.sasa0$pdb==prot & pdb.sasa0$nb > 25 #&TODO & sasa buried 
    rasmol.whole.chain(prot,file,pdbdir,I)
  }

}
