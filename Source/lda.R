library(topicmodels)
library(tidytext)

LDA2<-function(Y,set_seed=923,ntopics=10){
  start_time<-Sys.time()
  
  if(!is.null(set_seed)){
    est.lda<-LDA(t(Y),k=10,control=list(seed=set_seed),method="Gibbs")  
  }else{
    est.lda<-LDA(t(Y),k=10,method="Gibbs")
  }
  
  
  X.lda <- t(est.lda@beta)
  B.lda <- est.lda@gamma
  end_time<-Sys.time()
  
  
  colnames(X.lda)<-paste0("b",seq(from=1,to=ntopics))
  rownames(X.lda)<-rownames(Y)
  
  colnames(B.lda)<-paste0("b",seq(from=1,to=ntopics))
  rownames(B.lda)<-colnames(Y)
  
  time_duration<-difftime(end_time, start_time, units='mins')
  
  return(list("X.lda"=X.lda,"B.lda"=B.lda,"lda"=est.lda,"comp_time"=time_duration,
              "ntopics"=ntopics))
}

predict.LDA2<-function(Y,model,ntopics=10){
  
  start_time<-Sys.time()
  est.lda<-posterior(model$lda,t(Y))
  
  B<-est.lda$topics
  colnames(B)<-paste0("b",seq(from=1,to=model$ntopics))
  end_time<-Sys.time()
  
  time_duration<-difftime(end_time, start_time, units='mins')
  
  return(list("B.lda"=B,"comp_time"=time_duration))
}