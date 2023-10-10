# functions for training of latent semantic analysis

LSA<-function(Y, transform="none", ntopics=10){
  
  start_time<-Sys.time()
  
  if(transform=="tf-idf"){
    Y2<-Y
  }else{
    Y2<-Y
  }
  
  svd.dtm<-svd(Y2)
  
  X.lsa<-svd.dtm$u[,1:ntopics]
  B.lsa<-svd.dtm$v[,1:ntopics]
  
  end_time<-Sys.time()
  
  time_duration<-difftime(end_time, start_time, units='mins')
  
  colnames(X.lsa)<-paste0("b",seq(from=1,to=ntopics))
  rownames(X.lsa)<-rownames(Y)
  
  colnames(B.lsa)<-paste0("b",seq(from=1,to=ntopics))
  rownames(B.lsa)<-colnames(Y)
  
  
  return(
    list(
      "X.lsa"=X.lsa,
      "B.lsa"=B.lsa,
      "comp_time"=time_duration,
      "svd"=svd.dtm,
      "ntopics"=ntopics
    )
  )
}


predict.LSA2<-function(Y, model, transform="none",tol=1e-4){
  
  # get model
  ntopics = ncol(model$X.lsa)
  X = model$X.lsa %*% diag(model$svd$d[1:ntopics],nrow=ntopics,ncol=ntopics)
  
  #initialize matrix factors
  
  #make a vector of random seeds
  seeds<-round(100*runif(100))
  
  #define function to produce random initializations
  get_init_err<-function(x){
    set.seed(x)
    
    B<-matrix(rnorm(ncol(X)*ncol(Y),0,1),
                  nrow=ncol(X),ncol=ncol(Y),byrow=TRUE)
    
    err(Y,X,B)
  }
  
  #perform initialization for each random seed
  init_errs<-sapply(seeds,get_init_err,simplify=TRUE)
  
  #chosen seed and initialization is that with the lowest initialization error
  set.seed(seeds[init_errs==min(init_errs)])
  
  B.0<-matrix(rnorm(ncol(X)*ncol(Y),0,1),
                  nrow=ncol(X),ncol=ncol(Y),byrow=TRUE)

  #begin factorization
  errors<-numeric(maxiter)
  start_time<-Sys.time()
  
  for(i in 1:maxiter){
    
    #update
    B.denom<-(t(X) %*% X %*% B.0)
    B.denom<-B.denom + 1*(B.denom==0)
    B.1<-B.0 * (t(X) %*% Y) / B.denom
    
    errors[i]<-err(Y,X,B.1)
    niter<-i
    
    if(i>=burn){delta<-abs( errors[i]-errors[i-1] )/max(1,errors[i-1]) }else{ delta<-9999 }
    
    if(delta<tol){
      errors<-errors[1:niter]
      break
    }else{
      B.0<-B.1
    }
  }
  
  end_time<-Sys.time()
  time_duration<-end_time-start_time
  
  return(list(
    "B.plsa"=t(B.1),
    "errors"=errors,
    "comp_time"=time_duration,
    "niter"=niter
  ))
}


predict.LSA<-function(Y, model, transform="none",tol=1e-4,
                      burn_iter=100,max_iter=1000,set_seed=923){
  
  # get model
  ntopics = ncol(model$X.lsa)
  X2 = model$X.lsa %*% diag(model$svd$d[1:ntopics],nrow=ntopics,ncol=ntopics)
  
  #perform normal factorization
  est.lsa<-mat.factor(Y=Y,X=X2,nmf=FALSE,burn=burn_iter,maxiter=max_iter,seed=set_seed)
  
  B.lsa=t(est.lsa$B)
  colnames(B.lsa)<-paste0("b",seq(from=1,to=ntopics))
  
  return(list(
    "B.lsa"=B.lsa,
    "errors"=est.lsa$errors,
    "comp_time"=est.lsa$time,
    "niter"=est.lsa$niter
  ))
}
