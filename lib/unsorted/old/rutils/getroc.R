function(x,y) {
  predictions = c(x,y)
  labels      = c(rep(T,length(x)),rep(F,length(y)))
  pred <- prediction(predictions, labels)
  perf <- performance(pred, measure = "auc")
  return(slot(perf,"y.values")[[1]])
}
