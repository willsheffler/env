gvk <- read.table('gvkb_anb5.data',header=F)[,-1]
names(gvk) <- c('pdb','res1','aa1','atom1','type1','anb',
                         'res2','aa2','atom2','type2','dis')

rasmol.surrounding <- function(pdb,file,cres,catom,sres,satom) {
  # write rasmol script
  write('zap',file=file)
  write("stereo -5",file=file,append=T)
  write("set picking center",file=file,append=T)
  write(paste('load pdb ',pdb,sep=""),file=file,append=T)
  write("restrict none",file=file,append=T)
  write("select backbone",file=file,append=T)
  write("strands",file=file,append=T)
  write("select all",file=file,append=T)
  write("wireframe",file=file,append=T)
  write(paste("select (",cres,' and *.',catom,')',sep=''),file=file,append=T)
  write("spacefill",file=file,append=T)
  write("color yellow",file=file,append=T)
  write("center selected",file=file,append=T)
  for(ii in 1:length(sres)) {
    sel <- paste(paste("(",sres[ii],' and ','*.',satom[ii],')',sep=''))
    write(paste('select',sel),file=file,append=T)
    write("spacefill",file=file,append=T)
  }  
}

dbbh <- density(gvk$anb[gvk$type1==25])
dbbo <- density(gvk$anb[gvk$type1==20])
plot(dbbo,type='n',xlim=c(min(gvk$anb),max(gvk$anb)),ylim=c(0,max(dbbo$y,dbbh$y)))
lines(dbbo,col=1)
lines(dbbh,col=2)
      
file   <- "~/project/features/buried_polar/visualizations/1gvkB_all_anb5_bbo.rasmol"
pdb <- "~/project/features/buried_polar/1gvk_0001.pdb"
centerres  <- 10
centeratom <- 'O'
I <- gvk$res1==centerres & gvk$atom1=='O'                   
rasmol.surrounding(pdb,file,centerres,'O',gvk$res2[I],gvk$atom2[I])
