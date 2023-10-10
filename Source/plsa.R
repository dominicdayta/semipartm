# Normalize on columns

# plsa train on in-sample dtm
PLSA<-function(Y,burn_iter=100,max_iter=1000,set_seed=923, ntopics=10, tol=1e-3,chi=10){
  
  #perform normalized factorization
  est.plsa<-mat.factor(Y,normalize=TRUE,
                       burn=burn_iter,maxiter=max_iter,seed=set_seed,tol=tol,chi=chi)
  
  B.plsa<-t(est.plsa$B)
  X.plsa<-est.plsa$X
  
  colnames(X.plsa)<-paste0("b",seq(from=1,to=ntopics))
  rownames(X.plsa)<-rownames(Y)
  
  colnames(B.plsa)<-paste0("b",seq(from=1,to=ntopics))
  rownames(B.plsa)<-colnames(Y)
  
  return(list(
    "X.plsa"=X.plsa,
    "B.plsa"=B.plsa,
    "errors"=est.plsa$errors,
    "comp_time"=est.plsa$time,
    "niter"=est.plsa$niter,
    "ntopics"=ntopics
  ))
}

#plsa train on out-sample dtm
predict.PLSA<-function(Y,model,burn_iter=100,max_iter=1000,set_seed=923,tol=1e-3,chi=10){
  
  #get data
  Y2<-mat.normalize(Y)
  
  #perform normalized factorization
  est.plsa<-mat.factor(Y2,X=model$X.plsa,normalize=TRUE,
                       burn=burn_iter,maxiter=max_iter,seed=set_seed,tol=tol,chi=chi)
  
  B.plsa<-t(est.plsa$B)
  colnames(B.plsa)<-paste0("b",seq(from=1,to=ntopics))
  rownames(B.plsa)<-colnames(Y)
  
  return(list(
    "B.plsa"=B.plsa,
    "errors"=est.plsa$errors,
    "comp_time"=est.plsa$time,
    "niter"=est.plsa$niter
  ))
}