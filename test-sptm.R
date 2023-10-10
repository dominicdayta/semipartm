####################################################################
# Combined Codes for Thesis, Version 1
# By Dominic Dayta
#
# beta values refer to topic expressions within matrices
# z values refer to external variables as defined
# x values refer to word importance per topics (set as fixed)
#
# Generated data are run through four methods:
# 1. LSA
# 2. PLSA
# 3. LDA
# 4. Semiparametric TM (proposed)
#
# Each method produces two important values: X and B, where
# X is W x T, a matrix of topic structures, while
# B is T x D, a matrix of topic expressions in the documents
# Y is then XB -> W x D
#
####################################################################

#new, for getting of new packages
#suppressWarnigs(options(warn=-1))
#
#install.packages(c("topicmodels","tidytext","e1071","splines","MASS"))

################# Set working directory, load necessary directories
setwd("/Users/dominic/Documents/semipartm_thesis")
#setwd("~/Dropbox/SemiparTM")
library(topicmodels)
library(e1071)
library(splines)
library(MASS)

################ Load codes for specialized functions
source("Source/fac-2.R") #matrix operations
source("Source/gen_dtm.R") #generate data
source("Source/lsa.R") #latent semantic analysis
source("Source/plsa.R") #probabilistic LSA
source("Source/lda.R") #latent dirichlet allocation
source("Source/semipartm.R") #semiparametric topic modelling

################ Constants
seed<-923 #seed
k<-10 #number of topics
ntopics<-10

################ Test Proper
test_sptm<-function(nword,ndoc,spar,m,niter, fastmode = FALSE){
  
  start.time<-Sys.time()
  run_name<-paste(ndoc,"doc_",nword,"word_",spar,"spar_","m",m,sep="")
  
  #list of selected cv
  cv.list<-numeric(niter)
  
  #dataframe of accuracy statistics for classification
  class.stat<-data.frame(
    iter=numeric(0),
    method=numeric(0), # 0 : lsa; 1: plsa; 2: lda; 3: sptm1; 4: sptm3; sptm-cv
    in_test=numeric(0), # 0: on training set; 1: test set
    misclass = numeric(0), #misclassification rate
    precision = numeric(0),
    recall = numeric(0),
    specificity = numeric(0),
    stringsAsFactors = FALSE
  )
  
  #dataframe for method performance
  method.stat<-data.frame(
    iter=numeric(0),
    method=numeric(0), # same as above
    in_test = numeric(0), # same as above
    run_time=numeric(0), # total runtime
    xcos = numeric(0), # average cosine similarity on the X matrices
    bcos = numeric(0), #average cosine similarity on the B matrices
    
    #inidividual topics in the x matrix
    x01 = numeric(0),
    x02 = numeric(0),
    x03 = numeric(0),
    x04 = numeric(0),
    x05 = numeric(0),
    x06 = numeric(0),
    x07 = numeric(0),
    x08 = numeric(0),
    x09 = numeric(0),
    x10 = numeric(0),
    
    #individual topics in the b matrix
    b01 = numeric(0),
    b02 = numeric(0),
    b03 = numeric(0),
    b04 = numeric(0),
    b05 = numeric(0),
    b06 = numeric(0),
    b07 = numeric(0),
    b08 = numeric(0),
    b09 = numeric(0),
    b10 = numeric(0),
    stringsAsFactors = FALSE
  )
  
  twolevel.B<-data.frame(
    iter = numeric(0),
    epsilon = numeric(0), # as 1, 3, or -9 for cv
    b01 = numeric(0),
    b02 = numeric(0),
    b03 = numeric(0),
    b04 = numeric(0),
    b05 = numeric(0),
    b06 = numeric(0),
    b07 = numeric(0),
    b08 = numeric(0),
    b09 = numeric(0),
    b10 = numeric(0),
    stringsAsFactors = FALSE
  )
  
  set.seed(seed)
  
  for(i in 1:niter){
    
    message(paste0("[Iteration ",i,"]"))
    
    #----generate document term matrix
    # X is the nword x 10
    # Y is nword x ndoc
    # Z is ndoc x 5
    dtm<-gen_dtm(ndoc,nword,spar,m)
    
    #----train LSA
    #in-sample
    lsa<-LSA(Y=dtm$Y, transform="none")
    svm.lsa.is<-svm(factor(gamma)~b1 + b2 + b3 + b4 + b5 + b6 + b7 + b8 + b9 + b10,
                    data=data.frame(gamma=dtm$gam,lsa$B.lsa))
    
    #out-sample
    lsa.os<-predict.LSA(dtm$Yp,lsa)
    svm.lsa.os<-predict(svm.lsa.is,newdata=data.frame(lsa.os$B.lsa))
    
    #transform resulting X and B matrices for comparison
    X.lsa.is<-mat.normalize(makenonNegative(lsa$X.lsa))
    B.lsa.is<-mat.normalize(makenonNegative(lsa$B.lsa))
    B.lsa.os<-mat.normalize(makenonNegative(lsa.os$B.lsa))

 
    #----train PLSA
    plsa<-PLSA(Y=dtm$Y,
               burn_iter=100,max_iter=1000,set_seed=NULL, ntopics=10, tol=1e-3,chi=10)
    svm.plsa.is<-svm(factor(gamma)~b1 + b2 + b3 + b4 + b5 + b6 + b7 + b8 + b9 + b10,
                     data=data.frame(gamma=dtm$gam,plsa$B.plsa))
    
    #out-sample
    plsa.os<-predict.PLSA(dtm$Yp,plsa,burn_iter=100,max_iter=1000,set_seed=NULL,tol=1e-3,chi=10)
    svm.plsa.os<-predict(svm.plsa.is,newdata=data.frame(plsa.os$B.plsa))
    
    #transform resulting X and B matrices for comparison
    X.plsa.is<-mat.normalize(makenonNegative(plsa$X.plsa))
    B.plsa.is<-mat.normalize(makenonNegative(plsa$B.plsa))
    B.plsa.os<-mat.normalize(makenonNegative(plsa.os$B.plsa))
    
    
    #----train LDA
    lda<-LDA2(dtm$Y,set_seed=NULL,ntopics=10)
    svm.lda.is<-svm(factor(gamma)~b1 + b2 + b3 + b4 + b5 + b6 + b7 + b8 + b9 + b10,
                     data=data.frame(gamma=dtm$gam,lda$B.lda))
    
    #out-sample
    lda.os<-predict.LDA2(dtm$Yp,lda)
    svm.lda.os<-predict(svm.lda.is,newdata=data.frame(lda.os$B.lda))

    #transform resulting X and B matrices for comparison
    X.lda.is<-mat.normalize(makenonNegative(lda$X.lda))
    B.lda.is<-mat.normalize(makenonNegative(lda$B.lda))
    B.lda.os<-mat.normalize(makenonNegative(lda.os$B.lda))
    
    
    #----train SemiparTM epsilon = 1
    sptm1<-SPTM(dtm$Y,dtm$Z,epsilon=1.0,ntopics=10,reg_type="splines",
               tol=1e-3,max_iter=1000,burn_iter=100,set_seed=923)
    
    svm.sptm1.is<-svm(factor(gamma)~b1 + b2 + b3 + b4 + b5 + b6 + b7 + b8 + b9 + b10,
                    data=data.frame(gamma=dtm$gam,sptm1$B.sptm))
    
    #out-sample
    sptm1.os<-predict.SPTM(dtm$Zp,sptm1)
    svm.sptm1.os<-predict(svm.sptm1.is,newdata=data.frame(sptm1.os$B.sptm))
    
    #transform resulting X and B matrices for comparison
    X.sptm1.is<-mat.normalize(makenonNegative(sptm1$X.sptm))
    B.sptm1.is<-mat.normalize(makenonNegative(sptm1$B.sptm))
    B.sptm1.os<-mat.normalize(makenonNegative(sptm1.os$B.sptm))
    
    
    #----train SemiparTM epsilon = 3
    sptm3<-SPTM(dtm$Y,dtm$Z,epsilon=3,ntopics=10,reg_type="splines",
                tol=1e-3,max_iter=1000,burn_iter=100,set_seed=923)
    
    svm.sptm3.is<-svm(factor(gamma)~b1 + b2 + b3 + b4 + b5 + b6 + b7 + b8 + b9 + b10,
                      data=data.frame(gamma=dtm$gam,sptm3$B.sptm))
    
    #out-sample
    sptm3.os<-predict.SPTM(dtm$Zp,sptm3)
    svm.sptm3.os<-predict(svm.sptm3.is,newdata=data.frame(sptm3.os$B.sptm))
    
    #transform resulting X and B matrices for comparison
    X.sptm3.is<-mat.normalize(makenonNegative(sptm3$X.sptm))
    B.sptm3.is<-mat.normalize(makenonNegative(sptm3$B.sptm))
    B.sptm3.os<-mat.normalize(makenonNegative(sptm3.os$B.sptm))
    
    
    #----train SemiparTM cross-validated
    #for fastmode = TRUE
    if(fastmode){
      
      if(i == 1){
        cv.cur = "cv"
        cv.set = "cv"
        cv.message = "Cross-validation selected penalty of "
        
      }else{
        cv.set = cv.cur
        cv.message = "Fast Mode: Retaining previous cross-validation penalty of "
      }
      
    }else{
      cv.set = "cv"
      cv.message = "Cross-validation selected penalty of "
    }
    
    sptm.cv<-SPTM(dtm$Y,dtm$Z,epsilon=cv.set,ntopics=10,reg_type="splines",
                tol=1e-3,max_iter=1000,burn_iter=100,set_seed=923)
    cv.list[i]<-sptm.cv$penalty
    cv.cur <- sptm.cv$penalty
    
    message(paste0(cv.message,cv.list[i]))
    
    svm.sptm.cv.is<-svm(factor(gamma)~b1 + b2 + b3 + b4 + b5 + b6 + b7 + b8 + b9 + b10,
                      data=data.frame(gamma=dtm$gam,sptm.cv$B.sptm))
    
    #out-sample
    sptm.cv.os<-predict.SPTM(dtm$Zp,sptm.cv)
    svm.sptm.cv.os<-predict(svm.sptm.cv.is,newdata=data.frame(sptm.cv.os$B.sptm))
    
    #transform resulting X and B matrices for comparison
    X.sptm.cv.is<-mat.normalize(makenonNegative(sptm.cv$X.sptm))
    B.sptm.cv.is<-mat.normalize(makenonNegative(sptm.cv$B.sptm))
    B.sptm.cv.os<-mat.normalize(makenonNegative(sptm.cv$B.sptm))
    
    
    #----save results to dataframe
    class.stat<-rbind(class.stat,
                      
          #lsa
          cbind(iter=i,method=0,in_test=1,accMeasure(dtm$gam,svm.lsa.is$fitted)),
          cbind(iter=i,method=0,in_test=0,accMeasure(dtm$gamp,svm.lsa.os)),
            
          #plsa
          cbind(iter=i,method=1,in_test=1,accMeasure(dtm$gam,svm.plsa.is$fitted)),
          cbind(iter=i,method=1,in_test=0,accMeasure(dtm$gamp,svm.plsa.os)),
          
          #lda
          cbind(iter=i,method=2,in_test=1,accMeasure(dtm$gam,svm.lda.is$fitted)),
          cbind(iter=i,method=2,in_test=0,accMeasure(dtm$gamp,svm.lda.os)),
        
          #sptm1
          cbind(iter=i,method=3,in_test=1,accMeasure(dtm$gam,svm.sptm1.is$fitted)),
          cbind(iter=i,method=3,in_test=0,accMeasure(dtm$gamp,svm.sptm1.os)),
                      
          #sptm3
          cbind(iter=i,method=4,in_test=1,accMeasure(dtm$gam,svm.sptm3.is$fitted)),
          cbind(iter=i,method=4,in_test=0,accMeasure(dtm$gamp,svm.sptm3.os)),
      
          #sptm-cv
          cbind(iter=i,method=5,in_test=1,accMeasure(dtm$gam,svm.sptm.cv.is$fitted)),
          cbind(iter=i,method=5,in_test=0,accMeasure(dtm$gamp,svm.sptm.cv.os))
                      
    )
    
    method.stat<-rbind(method.stat,
                       
         #lsa
         cbind(iter=i,method=0,in_test=1,run_time=lsa$comp_time,
               xcos=matrixCos(mat.normalize(dtm$X),X.lsa.is),
               bcos=matrixCos(mat.normalize(dtm$B),B.lsa.is),
               t(matrixCos(mat.normalize(dtm$X),X.lsa.is, average=FALSE)),
               t(matrixCos(mat.normalize(dtm$B),B.lsa.is, average=FALSE)) 
         ),
         cbind(iter=i,method=0,in_test=0,run_time=lsa.os$comp_time,
             xcos=NA,
             bcos=matrixCos(mat.normalize(dtm$Bp),B.lsa.os),
              t(rep(NA,times=10)),
              t(matrixCos(mat.normalize(dtm$Bp),B.lsa.os, average=FALSE))
         ),
                       
         #plsa
         cbind(iter=i,method=1,in_test=1,run_time=plsa$comp_time,
            xcos=matrixCos(mat.normalize(dtm$X),X.plsa.is),
            bcos=matrixCos(mat.normalize(dtm$B),B.plsa.is),
            t(matrixCos(mat.normalize(dtm$X),X.plsa.is, average=FALSE)),
            t(matrixCos(mat.normalize(dtm$B),B.plsa.is, average=FALSE)) 
         ),
         cbind(iter=i,method=1,in_test=0,run_time=plsa.os$comp_time,
             xcos=NA,
             bcos=matrixCos(mat.normalize(dtm$Bp),B.plsa.os),
             t(rep(NA,times=10)),
             t(matrixCos(mat.normalize(dtm$Bp),B.plsa.os, average=FALSE))
         ),
                       
         #lda
         cbind(iter=i,method=2,in_test=1,run_time=lda$comp_time,
             xcos=matrixCos(mat.normalize(dtm$X),X.lda.is),
             bcos=matrixCos(mat.normalize(dtm$B),B.lda.is),
             t(matrixCos(mat.normalize(dtm$X),X.lda.is, average=FALSE)),
             t(matrixCos(mat.normalize(dtm$B),B.lda.is, average=FALSE)) 
         ),
         cbind(iter=i,method=2,in_test=0,run_time=lda.os$comp_time,
             xcos=NA,
             bcos=matrixCos(mat.normalize(dtm$Bp),B.lda.os),
             t(rep(NA,times=10)),
             t(matrixCos(mat.normalize(dtm$Bp),B.lda.os, average=FALSE))
         ),
                       
         #sptm1
         cbind(iter=i,method=3,in_test=1,run_time=sptm1$comp_time,
             xcos=matrixCos(mat.normalize(dtm$X),X.sptm1.is),
             bcos=matrixCos(mat.normalize(dtm$B),B.sptm1.is),
             t(matrixCos(mat.normalize(dtm$X),X.sptm1.is, average=FALSE)),
             t(matrixCos(mat.normalize(dtm$B),B.sptm1.is, average=FALSE)) 
         ),
         cbind(iter=i,method=3,in_test=0,run_time=sptm1.os$comp_time,
             xcos=NA,
             bcos=matrixCos(mat.normalize(dtm$Bp),B.sptm1.os),
             t(rep(NA,times=10)),
             t(matrixCos(mat.normalize(dtm$Bp),B.sptm1.os, average=FALSE))
         ),
                   
         #sptm3
         cbind(iter=i,method=4,in_test=1,run_time=sptm1$comp_time,
             xcos=matrixCos(mat.normalize(dtm$X),X.sptm3.is),
             bcos=matrixCos(mat.normalize(dtm$B),B.sptm3.is),
             t(matrixCos(mat.normalize(dtm$X),X.sptm3.is, average=FALSE)),
             t(matrixCos(mat.normalize(dtm$B),B.sptm3.is, average=FALSE)) 
         ),
         cbind(iter=i,method=4,in_test=0,run_time=sptm1.os$comp_time,
             xcos=NA,
             bcos=matrixCos(mat.normalize(dtm$Bp),B.sptm1.os),
             t(rep(NA,times=10)),
             t(matrixCos(mat.normalize(dtm$Bp),B.sptm1.os, average=FALSE))
         ),
                     
                     
         #sptm-cv
        cbind(iter=i,method=5,in_test=1,run_time=sptm.cv$comp_time,
             xcos=matrixCos(mat.normalize(dtm$X),X.sptm.cv.is),
             bcos=matrixCos(mat.normalize(dtm$B),B.sptm.cv.is),
             t(matrixCos(mat.normalize(dtm$X),X.sptm.cv.is, average=FALSE)),
             t(matrixCos(mat.normalize(dtm$B),B.sptm.cv.is, average=FALSE)) 
         ),
        cbind(iter=i,method=5,in_test=0,run_time=sptm.cv.os$comp_time,
             xcos=NA,
             bcos=matrixCos(mat.normalize(dtm$Bp),B.sptm.cv.os),
             t(rep(NA,times=10)),
             t(matrixCos(mat.normalize(dtm$Bp),B.sptm.cv.os, average=FALSE))
         )
    )
    
    twolevel.B[nrow(twolevel.B)+1,]<-c(i,sptm1$penalty,matrixCos(dtm$B, sptm1$B.reg, average=FALSE)/matrixCos(dtm$B, sptm1$B.sptm, average=FALSE))
    twolevel.B[nrow(twolevel.B)+1,]<-c(i,sptm3$penalty,matrixCos(dtm$B, sptm3$B.reg, average=FALSE)/matrixCos(dtm$B, sptm3$B.sptm, average=FALSE))
    twolevel.B[nrow(twolevel.B)+1,]<-c(i,sptm.cv$penalty,matrixCos(dtm$B, sptm.cv$B.reg, average=FALSE)/matrixCos(dtm$B, sptm.cv$B.sptm, average=FALSE))
    
  }
  
  comp.time<-difftime(Sys.time(),start.time,units="mins")
  
  # a little housekeeping
  colnames(method.stat)[7:16]<-paste("x",1:10,sep="")
  colnames(method.stat)[17:26]<-paste("b",1:10,sep="")
  
   saveRDS(list("class.stat"=class.stat,"method.stat"=method.stat,"twolevel.B"=twolevel.B,
   "cv.list"=cv.list,"nword"=nword,"ndoc"=ndoc,
   "spar"=spar,"m"=m,"comp.time"=comp.time),
      paste0("Notebooks/Data/",run_name,".RDS"))
   
   # saveRDS(list("class.stat"=class.stat,"method.stat"=method.stat,"twolevel.B"=twolevel.B,
   #              "cv.list"=cv.list,"nword"=nword,"ndoc"=ndoc,
   #              "spar"=spar,"m"=m,"comp.time"=comp.time),
   #         paste0("Outputs/",run_name,".RDS"))
  
  # return(list("class.stat"=class.stat,"method.stat"=method.stat,"twolevel.B"=twolevel.B,
  #             "cv.list"=cv.list,"nword"=nword,"ndoc"=ndoc,
  #             "spar"=spar,"m"=m,"comp.time"=comp.time))
}


#r1<-test_sptm(nword=300,ndoc=150,spar=0.99,niter=1,cv.switch="on")
  
#---- test proper run
#words 500 1500 3500
#docs 150 1000 3000
#spar 0.7, 0.9, 0.99

#test_sptm(nword=100,ndoc=100,spar=0.9,m=1,niter=2)

# m=1
# nwords=c(500)
# ndocs=c(150,1500,3000)
# spar=c(0.7,0.9,0.99)
# 
# for(w in nwords){
#   for(d in ndocs){
#     for(s in spar){
#       
#       test_sptm(nword=w,ndoc=d,spar=s,m=1,niter=100)
#       
#     }
#   }
# }
