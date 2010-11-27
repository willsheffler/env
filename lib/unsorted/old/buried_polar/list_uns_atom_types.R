# assume uns is a table of UNS data as output by decoystats

res1 <- c("A","C","D","E","F","G","H","I","K","L","M","N","P","Q","R","S","T","V","W","Y")
res3 <- c("ALA","CYS","ASP","GLU","PHE","GLY","HIS","ILE","LYS","LEU",
          "MET","ASN","PRO","GLN","ARG","SER","THR","VAL","TRP","TYR")

for(i in 1:20) {
  r <- res1[i]
  print(paste(res3[i],res1[i]))
  tmp <- as.factor(as.character(uns[uns[,4]==r,8]))
  print(levels( tmp ))
}
