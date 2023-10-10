# sample reading of result
setwd("/home/dominic/Documents/Thesis-SPTM")
result<-readRDS("Notebooks/Data/100doc_100word_0.9spar_m1.RDS")

library(ggplot2)

accuracy<-aggregate(result$class.stat[,4:7],
          by=list("method" = result$class.stat$method,"train"=result$class.stat$in_test),
          function(x){ mean(x, na.rm=TRUE)})

overall_cossim<-aggregate(result$method.stat[,5:6],
          by=list("method" = result$method.stat$method,"train"=result$method.stat$in_test),
          function(x){ mean(x, na.rm=TRUE)})

x.cossim<-aggregate(result$method.stat[,7:16],
          by=list("method" = result$method.stat$method,"train"=result$method.stat$in_test),
          function(x){ mean(x, na.rm=TRUE)})

b.cossim<-aggregate(result$method.stat[,17:26],
          by=list("method" = result$method.stat$method,"train"=result$method.stat$in_test),
          function(x){ mean(x, na.rm=TRUE)})

run_time<-aggregate(result$method.stat[,4],
          by=list("method" = result$method.stat$method,"train"=result$method.stat$in_test),
          function(x){ mean(x, na.rm=TRUE)})

accuracy[accuracy$train==1,]
accuracy[accuracy$train==0,]

x.cossim[x.cossim$train==1,]

b.cossim[x.cossim$train==1,]
b.cossim[x.cossim$train==0,]

overall_cossim[overall_cossim$train==1,]
overall_cossim[overall_cossim$train==0,]
