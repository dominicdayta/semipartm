# codes for matrix operations performed within the codes
# experimental

matrixCos<-function(X,Y,average=TRUE){
  
  if(average){
    mean(
      sapply((1:ncol(X)),function(x){
        vectorCos(X[,x],Y[,x])
      })
    )  
  }else{
    sapply((1:ncol(X)),function(x){
      vectorCos(X[,x],Y[,x])
    })
  }
  
}

vectorCos<-function(x,y){
  d<-sqrt(sum(x**2))*sqrt(sum(y**2))
  d<-d*(d>=0)+1*(d==0)
  
  sum(x*y)/d
}

matrixDif<-function(X,Y){
  
  apply(
    sapply((1:nrow(X)),function(x){
      abs(X[,x]-Y[,x])
    }),1,mean)    
  
}

makenonNegative<-function(X){
  
  minX = min(X,na.rm=TRUE)
  
  if(minX>=0){
    minX2 = 0
  }else{
    minX2 = -1*minX    
  }
  
  X + minX2
}

accMeasure<-function(true,fitted){
  misclass = sum(fitted!=true)/max(1,length(true))
  precision = sum((fitted==1)*(true==1))/max(1,sum(fitted==1))
  recall = sum((fitted==1)*(true==1))/max(1,sum(true==1))
  specificity = sum((fitted==0)*(true==0))/max(1,sum(true==0))
  cbind(misclass,precision,recall,specificity)
}

# error
err<-function(original,fac1,fac2,chi=1,type=1){
  
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

# factor a given matrix Y into Y = XB
mat.factor<-function(Y,X=NULL,nmf=TRUE,normalize=FALSE,e=0,
                     burn=100,maxiter=1000,tol=1e-3,ntopics=10,chi=1, seed=923){
  
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
  if(!is.null(seed)){ set.seed(seed) }
  
  # initialize matrices
  if(nmf){
    if(normalize){
      X.0<-mat.normalize(matrix(runif(ntopics*nrow(Y2)),
                              nrow=nrow(Y2),ncol=ntopics,byrow=TRUE))
      B.0<-t(mat.normalize(matrix(runif(ntopics*ncol(Y2)),
                                nrow=ncol(Y2),ncol=ntopics,byrow=TRUE)))  
    }else{
      X.0<-matrix(runif(ntopics*nrow(Y2)),
                nrow=nrow(Y2),ncol=ntopics,byrow=TRUE)
      B.0<-t(matrix(runif(ntopics*ncol(Y2)),
                  nrow=ncol(Y2),ncol=ntopics,byrow=TRUE))
    }
  }else{
    X.0<-matrix(rnorm(ntopics*ncol(Y2),0,1),
              nrow=ncol(Y2),ncol=ntopics,byrow=TRUE)
    
    B.0<-t(matrix(rnorm(ntopics*ncol(Y2),0,1),
                nrow=ncol(Y2),ncol=ntopics,byrow=TRUE))
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
      X.denom<-X.denom*(1 + e) + 1*(X.denom==0)
      X.1<-X.0 * (Y %*% t(B.0)) / X.denom
      if(normalize){ X.1<-mat.normalize(X.1) }  
    }else{
      X.1<-X
    }
    
    B.denom<-(t(X.1) %*% X.1 %*% B.0)
    B.denom<-B.denom*(1 + e) + 1*(B.denom==0)
    B.1<-B.0 * (t(X.1) %*% Y) / B.denom
    if(normalize){ B.1<-mat.normalize(B.1,transpose=TRUE) }
    
    errors[i]<-err(Y2,X.1,B.1,chi)
    if(i>=burn){delta<-abs( errors[i]-errors[i-1] )/max(1,errors[i-1]) }else{ delta<-9999 }
    niter<-i
    
    if(delta<tol){
      errors<-errors[1:niter]
      break
    }else{
      X.0<-X.1
      B.0<-B.1
    }
  }
  
  end_time<-Sys.time()
  time_duration<-time_duration<-difftime(end_time, start_time, units='mins')
  
  return(list(
    "errors"=errors,
    "error"=errors[niter],
    "niter"=niter,
    "X" = X.1,
    "B" = B.1,
    "time"= time_duration
  )) 
}
