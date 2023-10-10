# Generates document-term-matrix

# adds words into documents/topics to ensure no topic and no word have zero scores

word_insert<-function(Y,fudge=1,transpose=FALSE){
  if(transpose==TRUE){
    Y = t(Y)
  }
  
  nword<-nrow(Y)
  ndoc<-ncol(Y)
  word_sums<-apply(Y,2,sum)
  
  Y_new = Y
  
  #index where the random 1 word will be placed per document
  insert_index<-round(runif(ndoc,min=1,max=nword),0)
  for(x in 1:ndoc){
    Y_new[insert_index[x],x] =  Y[insert_index[x],x] + fudge*(word_sums[x] == 0)
  }
  
  if(transpose==TRUE){
    return(t(Y_new))
  }else{
    return(Y_new) 
  }
}

# function for generating the corpus Y

gen_dtm<-function(ndoc, nword, spar=0.9,m=1){
  #note, a total of ndoc + 50 documents will be generated, where the
  #additional 50 are for out-of-sample tests
  
  #----simulation of z variables, z1 to z5, dimension is (ndoc+100)x5
  z1<-rpois(ndoc+100,1)
  z2<-rnorm(ndoc+100,20,7)
  z3<-rbinom(ndoc+100,1,0.8)
  z4<-rbeta(ndoc+100,6,2)
  z5<-rbeta(ndoc+100,10,2)
  
  Z = cbind("z1"=z1,"z2"=z2,"z3"=z3,"z4"=z4,"z5"=z5)
  
  #----topic functions
  
  s = rbinom(ndoc+100,1,spar)
  e = rnorm(ndoc+100,0,1)
  
  # linear models  
  b1 <- (1-s)*(-1 + z1 + 0.2*z2 + z3 - 0.9*z4 - 2*z5 + m*e) 
  b2 <- (1-s)*(3 + 1.5*z1 + 0.15*z2 - 5*z3 - 5*z5 + m*e)
  b3 <- (1-s)*(2 + 0.2*z2 - 1.4*z1 + m*e)
  b4 <- (1-s)*(1.6*z1 + 8*z3 - 9*z4 + m*e)
  
  b1<-b1*(b1>=0) + 0
  b2<-b2*(b2>=0) + 0
  b3<-b3*(b3>=0) + 0
  b4<-b4*(b4>=0) + 0
  
  # nonlinear models
  b5<- (1-s)*(z1**2/(5*z5) + 3 + m*e)
  b6<- (1-s)*(6 + sin(z5*z1) + m*e)
  b7<- (1-s)*(2 + 3*z1*z4 - 2*z3+ m*e)
  
  b5<-b5*(b5>=0) + 0
  b6<-b6*(b6>=0) + 0
  b7<-b7*(b7>=0) + 0 
  
  # correlation models
  b8<- (1-s)*(1 + 10*z4 - 2*b3 + m*e)
  b9<- (1-s)*(0.2*z2 + 0.2*b7 + m*e)
  b10<- (1-s)*(5 + 0.9*b1 - 1.2*b7 + m*e)
  
  b8<-b8*(b8>=0) + 0
  b9<-b9*(b9>=0) + 0
  b10<-b10*(b10>=0) + 0
  
  #transform to topic matrix, dimension is (ndoc + 100)x10
  B<-cbind(b1,b2,b3,b4,b5,b6,b7,b8,b9,b10)
  
  #----simulation of document classification
  
  gam_prob<-1.2*b1 + 4.2*b2 - 3*b3 + 3.4*b4 + 3.1*b5 - 1.5*b6 + 2.4*b7 + 3.2*b8 + 2*b9 - 1*b10 + rnorm(ndoc+100,0,1)
  gam_prob<-exp(gam_prob)/(1 + exp(gam_prob))
  gam<-1*(gam_prob>=0.60)
  
  
  #-----dictionary
  X<-word_insert(matrix((1-s)*rpois(10*nword,100)/90,ncol=10,nrow=nword), transpose=TRUE)
  colnames(X)<-c("b1","b2","b3","b4","b5","b6","b7","b8","b9","b10")
  
  
  #----generate raw word-document matrix and adjust, dimension is nword x (ndoc + 100)
  Y<-word_insert(round(X%*%t(B),0))
  
  #----some housekeeping
  colnames(Y)<-paste0("d",seq(from=1,to=ndoc+100))
  rownames(Y)<-paste0("w",seq(from=1,to=nword))
  
  rownames(X)<-paste0("w",seq(from=1,to=nword))
  
  rownames(B)<-paste0("d",seq(from=1,to=ndoc+100))
  
  return(list("gam"=gam[1:ndoc],"X"=X,"Z"=Z[1:ndoc,],"Y"=Y[,1:ndoc],
              "B"=B[1:ndoc,],"Bp"=B[(ndoc+1):(ndoc+100),],
              "Zp"=Z[(ndoc+1):(ndoc+100),],"gamp"=gam[(ndoc+1):(ndoc+100)],
              "Yp"=Y[,(ndoc+1):(ndoc+100)]))
}
