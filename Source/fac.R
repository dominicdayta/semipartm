# codes for matrix operations performed within the codes

# error
err<-function(original,fac1,fac2,type=1,chi=1){
  
  original<-chi*original
  fac1<-chi*fac1
  fac2<-chi*fac2
  
  if(type==1){
    sum(((original - fac1 %*% fac2)**2)/(original+1*(original==0)))
  }else if(type==2){
    sum(((original - fac1 %*% fac2)**2)/(nrow(original)*ncol(original)))
  }
}

# normalize matrix
mat.normalize<-function(x, transpose=FALSE){
  
  if(transpose==FALSE){
    # if transpose = FALSE (default), col sums add to 1
    sums<-apply(x,2,sum)
    sums<-sums + 1*(sums==0)
    colsums<-diag(1/sums,nrow=ncol(x),ncol=ncol(x))  
    return(x%*%colsums)
    
  }else{
    # if tranpose = TRUE, row sums add to 1
    sums<-apply(x,1,sum)
    sums<-sums + 1*(sums==0)
    rowsums<-diag(1/sums,nrow=nrow(x),ncol=nrow(x))
    return(rowsums %*% x)
  }
  
}

# initialize matrices
mat.init<-function(normalize,seed=923,nmf=TRUE){
  
  set.seed(seed)
  
  if(nmf){
    if(normalize){
      X<-mat.normalize(matrix(runif(ntopics*nrow(Y2)),
                              nrow=nrow(Y2),ncol=ntopics,byrow=TRUE))
      B<-t(mat.normalize(matrix(runif(ntopics*ncol(Y2)),
                                nrow=ncol(Y2),ncol=ntopics,byrow=TRUE)))  
    }else{
      X<-matrix(runif(ntopics*nrow(Y2)),
                nrow=nrow(Y2),ncol=ntopics,byrow=TRUE)
      B<-t(matrix(runif(ntopics*ncol(Y2)),
                  nrow=ncol(Y2),ncol=ntopics,byrow=TRUE))
    }
  }else{
    X<-matrix(rnorm(ntopics*ncol(Y2),0,1),
                nrow=ncol(Y2),ncol=ntopics,byrow=TRUE)
    
    B<-t(matrix(rnorm(ntopics*ncol(Y2),0,1),
           nrow=ncol(Y2),ncol=ntopics,byrow=TRUE))
  }
  
  return(list("X"=X,"B"=B))
}

# factor a given matrix Y into Y = XB

mat.factor<-function(Y,X=NULL,nmf=TRUE,normalize=FALSE,epsilon=0,
                     burn=100,maxiter=1000,tol=1e-4,ntopics=10,chi=1){
  
  #--- normalization is not done alongside non-nmf
  if((!nmf) && normalize){
    normalize<-FALSE
    print("Warning: Normalization is only allowed for non-negative type factorization.")
  }
  
  
  #--- get matrix to factor
  if(normalize){
    Y2<-mat.normalize(Y)  
  }else{
    Y2<-Y
  }
  
  
  #--- initialize matrix factors
  #make a vector of random seeds
  seeds<-round(100*runif(100))  
  
  #define function to produce random initializations
  get_init_err<-function(x,X=NULL,chi=1,nmf=TRUE){
    
    init<-mat.init(normalize=normalize,seed=x,nmf=nmf)
    B = init$B
    
    if(is.null(X)){
      X = init$X
    }
    
    err(Y2,X,B,type=1,chi)
  }
  
  #perform initialization for each random seed
  init_errs<-sapply(seeds,function(x){
      get_init_err(x,X,chi,nmf=nmf)
    },simplify=TRUE)
  
  #chosen seed and initialization is that with the lowest initialization error
  init<-mat.init(normalize,seed=seeds[init_errs==min(init_errs)])
  B.0 = init$B
  if(is.null(X)){
    X.0 = init$X
  }else{
    X.0 = X
  }
  
  #begin factorization
  errors<-numeric(maxiter)
  start_time<-Sys.time()
  
  for(i in 1:maxiter){
    
    #update
    # not that liu et al's method is not well-suited for matrices
    # not are not constrained to the simplex range
    # I wonder why
    
    if(is.null(X)){
      X.denom<-(X.0 %*% B.0 %*% t(B.0))
      X.denom<-X.denom*(X.denom>0) + 1*(X.denom=0)
      X.1<-X.0 * (Y %*% t(B.0)) / X.denom
      if(normalize){ X.1<-mat.normalize(X.1) }  
    }else{
      X.1<-X
    }
    
    B.denom<-(t(X.1) %*% X.1 %*% B.0)
    B.denom<-B.denom*(B.denom>0) + 1*(B.denom=0)
    B.1<-B.0 * (t(X.1) %*% Y) / B.denom
    if(normalize){ B.1<-mat.normalize(B.1,transpose=TRUE) }
    
    errors[i]<-err(Y2,X.1,B.1,chi)
    if(i>=burn){delta<-abs( errors[i]-errors[i-1] )/max(1,errors[i-1]) }else{ delta<-9999 }
    niter<-i
    
    if(is.nan(delta)){
      break
      errors<-errors[1:niter]
    }

    if(delta<tol){
      errors<-errors[1:niter]
      break
    }else{
      X.0<-X.1
      B.0<-B.1
    }
  }
  
  end_time<-Sys.time()
  time_duration<-end_time-start_time
 
  return(list(
    "errors"=errors,
    "error"=errors[niter],
    "niter"=niter,
    "X" = X.1,
    "B" = B.1,
    "time"= time_duration
  )) 
}