write.data.frame <- function(data,file,pretty=F,len=4) {
  print(paste("writting data frame",file))
  n <- colnames(data)
  top <- as.data.frame(as.list(n))
  colnames(top) <- n
  for(ii in 1:length(n))
    top[,ii] <- as.character(top[,ii])
  data <- rbind(top,data)
  data <- cbind(rownames(data),data)
  if(pretty) {
    for(ii in 1:length(n)) {
      data[,ii] <- format(substr(as.character(data[,ii]),1,len),justify='left') 
    }
  }
  write.table(data,file=file,quote=F,row.names=F,col.names=F)
}

