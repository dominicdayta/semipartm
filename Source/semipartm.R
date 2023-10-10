# functions for training of semiparametric analysis

library(splines)
library(MASS)

# estimate splines model on a set of dependent variables B against a set of Z
get_model<-function(B,Z,method="splines"){
  n_ind = ncol(Z)
  n_dep = ncol(B)
  Zmeans<-apply(Z,2,mean)
  Zdevs<-apply(Z,2,sd)
  
  Z2<-data.frame(sapply(1:n_ind,function(x){
    dnorm(Z[,x],mean=Zmeans[x],sd=Zdevs[x])
  }))
  
  #list to return all models
  models<-vector(mode = "list", length = n_dep)
  
  #matrix of fitted values
  B.reg<-matrix(0,nrow=nrow(B),ncol=ncol(B))
  
  # if(method=="splines")
  for(i in 1:n_dep){
    b<-B[,i]
    colnames(Z2)<-paste0("z",1:n_ind)
    
    formula<-"b ~"
    for(j in 1:n_ind){
      formula<-paste0(formula," bs(z",j,", knots=c(from=-4,to=4,by=0.5), Boundary.knots=c(-6, 6), degree=3)")
      if(j < n_ind) formula<-paste0(formula," +")
    }
    b.model<-lm(as.formula(formula),data=data.frame(b,Z2))
    
    B.reg[,i]<-b.model$fitted.values
    models[[i]]<-b.model
  }
  
  return(list("normalize"=list("mu"=Zmeans,"sig"=Zdevs),
              "models"=models,"B.reg"=B.reg))
}

# predict regression models
predict_model<-function(Z,model){
  n_dep<-length(model$models)
  n_ind<-length(model$normalize$mu)
  
  Z2<-data.frame(sapply(1:n_ind,function(x){
    dnorm(Z[,x],mean=model$normalize$mu[x],sd=model$normalize$sig[x])
  }))
  colnames(Z2)<-paste0("z",1:n_ind)
  
  B<-matrix(0,nrow=nrow(Z),ncol=n_dep)
  for(i in 1:n_dep){
    B[,i]<-predict(model$models[[i]],newdata=Z2)
    B[,i]<-B[,i]*(B[,i]>=0) + 0
  }
  
  return(B)
}

# perform semiparametric topic modelling
SPTM<-function(Y,Z=NULL,epsilon=0,ntopics=10,reg_type="splines",
               tol=1e-3,max_iter=1000,burn_iter=100,set_seed=923){
  
  start.time<-Sys.time()
  
  # Catch data matrix
  if(is.null(Z)){
    stop("No data matrix Z provided")
  }else if(is.null(nrow(Z)) || nrow(Z) != ncol(Y)){
    stop("Corpus and data matrix dimensions do not match.")
  }
  
  # Get corpus properties
  Y2<-Y
  ndoc<-ncol(Y2)
  nword<-nrow(Y2)
  nvar<-ncol(Z)
  
  # Catch epsilon entry
  if(!is.numeric(epsilon)){
    if(epsilon=="cv"){
      to.cv<-TRUE
    }else{
      stop("Penalty parameter provided not recognized.")
    }
  }else{
    if(epsilon < 0){
      stop("Penalty parameter provided not recognized.")
    }else{
      to.cv<-FALSE
    }
  }
  
  # Now perform cross-validation
  if(to.cv){
    
    message("Performing semiparametric topic modelling with cross validation. This might take a while.")
    
    # cv.epsilon<-seq(from=0,to=5,by=0.1)
    cv.epsilon<-seq(from=4,to=5,by=0.1)
    cv.errors<-sapply(cv.epsilon,function(x){
      
      errs<-sapply(1:20,function(y){
        
        samp<-rbinom(ndoc,1,prob=0.7)
        train<-(1:ndoc)[samp==1]
        test<-(1:ndoc)[samp==0]
        
        Y.1<-Y2[,train]
        Y.2<-Y2[,test]
        
        Z.1<-Z[train,]
        Z.2<-Z[test,]
        
        # matrix factorization
        cv.nmf<-mat.factor(Y.1,X=NULL,nmf=TRUE,e=x,
                           burn=burn_iter,maxiter=max_iter,tol=tol,ntopics=ntopics, seed=set_seed)
        
        # regression
        cv.reg<-get_model(t(cv.nmf$B),Z.1,method=reg_type)
        
        # predict on test set
        B.cv<-predict_model(Z.2,cv.reg)
        
        # get errors
        err(Y.2,cv.nmf$X,t(B.cv),type=2)
      })
      
      mean(errs,na.rm=TRUE)
    })
    
    set_eps<-cv.epsilon[cv.errors==min(cv.errors, na.rm=TRUE)]
    df.cv<-data.frame(cv.epsilon,cv.errors)
    
  }else{
    
    message(paste0("Performing semiparametric topic modelling with epsilon = ",epsilon,"."))
    set_eps<-epsilon
    df.cv<-NULL
    
  }
  
  # matrix factorization
  est.nmf<-mat.factor(Y,X=NULL,nmf=TRUE,e=set_eps,
             burn=burn_iter,maxiter=max_iter,tol=tol,ntopics=ntopics, seed=set_seed)
  
  B.sptm<-t(est.nmf$B)
  X.sptm<-est.nmf$X
  
  colnames(B.sptm)<-paste0("b",1:ntopics)
  colnames(X.sptm)<-paste0("b",1:ntopics)
  
  rownames(B.sptm)<-colnames(Y)
  rownames(X.sptm)<-rownames(Y)
  
  # train semiparametric model
  est.reg<-get_model(B.sptm,Z,method=reg_type)
  B.reg<-est.reg$B.reg
  
  colnames(B.reg)<-paste0("b",1:ntopics)
  rownames(B.reg)<-colnames(Y)
  
  end.time<-Sys.time()
  
  return(list(
    "X.sptm"=X.sptm,
    "B.sptm"=B.sptm,
    "B.reg"=B.reg,
    "reg"=est.reg,
    "errors"=est.nmf$errors,
    "comp_time"=difftime(end.time, start.time, units='mins'),
    "niter"=est.nmf$niter,
    "ntopics"=ntopics,
    "penalty"=set_eps,
    "with.cv"=to.cv,
    "cv.errors"=df.cv
  ))
}

predict.SPTM<-function(Z,model){
  
  start.time<-Sys.time()
  est.sptm<-predict_model(Z,model$reg)
  end.time<-Sys.time()
  
  colnames(est.sptm)<-paste0("b",1:ntopics)
  
  return(list(
    "B.sptm"=est.sptm,
    "comp_time"=difftime(end.time, start.time, units='mins')
  )) 
}
