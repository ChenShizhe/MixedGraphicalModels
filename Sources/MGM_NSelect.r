###--------------------------###
### Neighbourhood selection
### Last updated: Mar.14th 2014
###--------------------------###
## This file contains the following functions:
##  1) combine_select: combination using the selection rule (for Poisson-binary)
##  2) combine: a general combination function
##  3) combine_type: recontruction of edges for each edge type


##########
## Required library:


library(glmnet)

###-------------------------------------###
## neighbour ####
## neighbourhood selection on Gaussian-binary networks
## Input: 
##      dat: the data matrix, where the first p columns are Gaussian nodes and the 
##           last q are binary
##      p: number of Gaussian nodes
##      q: number of binary nodes
##      clambda: a vector of constant c
##      ratio: the ratio between the tuning parameters between 
##             the binary nodes and the Gaussian nodes
##      pf: whether to use weights
##      pw: the weights of penalty on each coefficients
## Output: 
##      N: the whole path of estimated neighbourhoods
##      A: the whole path of estimated intercepts for all regression
##      nlg: negative log-likelihood
##      bic: the whole path of BIC for all regression
###-------------------------------------###

neighbour<-function(dat, p, q, clambda, ratio=1, pf=FALSE, pw){
  n<-dim(dat)[1]
  
  #determine the tuning parameters for Gaussian nodes and binary nodes
  lambda_g<-clambda*max(c(sqrt(2*log(p)/n),sqrt(2*log(q)/n)))  
  lambda_b<- ratio*lambda_g
  
  
  nlambda_g<-length(lambda_g)
  N<-replicate(nlambda_g, matrix(0,p+q,p+q))    
  A<-nlg<-bic<-matrix(0,p+q, nlambda_g)
  for(i in 1:p){  
    #setting different weights for each coefficients, using pw.
    if(pf==FALSE){
      penalty_weight=rep(1,p+q-1)
    } else{
      penalty_weight=pw[-i] 
    }
  
    # need to scale the covariate because glmnet will scale the tuning parameters
    
    scale_fac<-(p+q-1)/sum(penalty_weight) 
    penalty_weight<-scale_fac*penalty_weight
    reg<-glmnet(x=scale_fac*dat[,-i], y=(dat[,i]), family="gaussian", 
                lambda=lambda_g, standardize=F,penalty.factor=penalty_weight )
    
    N[i,-i,]<-rbind(-(as.matrix(reg$beta))[1:(p-1),], (as.matrix(reg$beta))[p:(p+q-1),])*scale_fac
    A[i,]<-as.matrix(reg$a0)
    
    for(j in 1: nlambda_g ){
      nlg[i,j]<-N2Loglklh(beta=as.matrix(reg$beta)[,j], a0=as.matrix(reg$a0)[j],
                          y=dat[,i], x=scale_fac*dat[,-i],dist='G')
      bic[i,j]<-BIC(nlg[i,j],n,as.matrix(reg$beta)[,j])
    }
    
    
  }
  for(i in p+(1:q)){
    if(pf==FALSE){
      penalty_weight=rep(1,p+q-1)
    } else{
       penalty_weight=pw[-i] 
    }
    scale_fac<-(p+q-1)/sum(penalty_weight) #this is the scaling factor
    penalty_weight<-scale_fac*penalty_weight
     reg<-glmnet(x=scale_fac*dat[,-i], y=(dat[,i]>mean(dat[,i])), family="binomial", 
                lambda=lambda_b, standardize=F,penalty.factor=penalty_weight )
    
    N[i,-i,]<- (as.matrix(reg$beta))*scale_fac  
    A[i,]<-as.matrix(reg$a0)
    for(j in 1: nlambda_g ){
      nlg[i,j]<-N2Loglklh(beta=as.matrix(reg$beta)[,j], a0=as.matrix(reg$a0)[j],
                          y=dat[,i], x=scale_fac*dat[,-i],dist='B')
      bic[i,j]<-BIC(nlg[i,j],n,as.matrix(reg$beta)[,j])
      
    }
    
  }
  return(list(N=N,Intercept=A,nlg=nlg,bic=bic))  
}


###-------------------------------------###
## neighbour_PB ####
## neighbourhood selection on Poisson-binary networks
## Input: 
##      dat: the data matrix, where the first p columns are Gaussian nodes and the 
##           last q are binary
##      p: number of Gaussian nodes
##      q: number of binary nodes
##      clambda: a vector of constant c
##      ratio: the ratio between the tuning parameters between 
##             the binary nodes and the Gaussian nodes
##      pf: whether to use weights
##      pw: the weights of penalty on each coefficients
##      maxit: the maximums number of interation for glmnet
## Output: 
##      N: the whole path of estimated neighbourhoods
##      A: the whole path of estimated intercepts for all regression
##      nlg: negative log-likelihood
##      bic: the whole path of BIC for all regression
###-------------------------------------###

neighbour_PB<-function(dat, p, q, clambda, ratio=FALSE, pf=FALSE, pw,maxit=10000){
  n<-dim(dat)[1]
  
  #determine the tuning parameters for Gaussian nodes and binary nodes
  lambda_g<-clambda*max(c(sqrt(2*log(p)/n),sqrt(2*log(q)/n)))  
  lambda_b<- ratio*lambda_g
  
  
  nlambda_g<-length(lambda_g)
  N<-replicate(nlambda_g, matrix(0,p+q,p+q))    
  A<-nlg<-bic<-matrix(0,p+q, nlambda_g)
  
  for(i in 1:p){
    
    #setting different weights for each coefficients, using standard deviations.
    if(pf==FALSE){
      penalty_weight=rep(1,p+q-1)
    } else{
     penalty_weight=pw[-i] 
    }
    scale_fac<-(p+q-1)/sum(penalty_weight) #this is the scaling factor
    penalty_weight<-scale_fac*penalty_weight
    reg<-glmnet(x=scale_fac*dat[,-i], y=(dat[,i]), family="poisson", 
                lambda=lambda_g, standardize=F,penalty.factor=penalty_weight,maxit=maxit )
    
    limi<-dim(as.matrix(reg$beta))[2]
    
    N[i,-i,1:limi]<-(as.matrix(reg$beta))*scale_fac  
    A[i,1:limi]<-as.matrix(reg$a0)
    
    for(j in 1: limi ){
      nlg[i,j]<-N2Loglklh(beta=as.matrix(reg$beta)[,j], a0=as.matrix(reg$a0)[j],
                          y=dat[,i], x=scale_fac*dat[,-i],dist="P")
      bic[i,j]<-BIC(nlg[i,j],n,as.matrix(reg$beta)[,j])
    }
    
    
  }
  for(i in p+(1:q)){
    if(pf==FALSE){
      penalty_weight=rep(1,p+q-1)
    } else{
    penalty_weight=pw[-i] 
    }
    scale_fac<-(p+q-1)/sum(penalty_weight)  
    penalty_weight<-scale_fac*penalty_weight
    reg<-glmnet(x=scale_fac*dat[,-i], y=(dat[,i]>mean(dat[,i])), family="binomial", 
                lambda=lambda_b, standardize=F,penalty.factor=penalty_weight,maxit=maxit )

    limi<-dim(as.matrix(reg$beta))[2]
    N[i,-i,1:limi]<-(as.matrix(reg$beta))*scale_fac  
    A[i,1:limi]<-as.matrix(reg$a0)
    for(j in 1: limi ){
      nlg[i,j]<-N2Loglklh(beta=as.matrix(reg$beta)[,j], a0=as.matrix(reg$a0)[j],
                          y=dat[,i], x=scale_fac*dat[,-i],dist='B')
      bic[i,j]<-BIC(nlg[i,j],n,as.matrix(reg$beta)[,j])
      
    }
    
  }
  return(list(N=N,Intercept=A,nlg=nlg,bic=bic))  
}



###-------------------------------------###
## neighbour_Gaussian ####
## neighbourhood selection on Gaussian graphical models
## Input: 
##      dat: the data matrix, where the first p columns are Gaussian nodes and the 
##           last q are binary
##      p: number of Gaussian nodes
##      q: number of binary nodes (not really useful)
##      clambda: a vector of constant c
##      ratio: the ratio between the tuning parameters between 
##             the binary nodes and the Gaussian nodes
##      pf: whether to use weights
##      pw: the weights of penalty on each coefficients
## Output: 
##      N: the whole path of estimated neighbourhoods
###-------------------------------------###
## Note:
## This function is modified based on neighbour(),
## and some trunks of code are unnecessary
###-------------------------------------###
neighbour_Gaussian<-function(dat, p, q, clambda,  pf=FALSE, pw){
  n<-dim(dat)[1]
  lambda_g<-clambda*max(c(sqrt(2*log(p)/n),sqrt(2*log(q)/n)))  
  nlambda_g<-length(lambda_g)
  
  N<-replicate(nlambda_g, matrix(0,p+q,p+q))    
  
  for(i in 1:p){
    
    if(pf==FALSE){
      penalty_weight=rep(1,p+q-1)
    } else{
    penalty_weight=pw[-i] 
    }
    
    #The penalty factors is scaled to p+q-1
    scale_fac<-(p+q-1)/sum(penalty_weight) 
    penalty_weight<-scale_fac*penalty_weight
    reg<-glmnet(x=scale_fac*dat[,-i], y=(dat[,i]), family="gaussian", 
                lambda=lambda_g, standardize=F,penalty.factor=penalty_weight)
    N[i,-i,]<- rbind(-(as.matrix(reg$beta))[1:(p-1),], (as.matrix(reg$beta))[p:(p+q-1),])*scale_fac  
  }
  for(i in p+(1:q)){
    #setting different weights for each coefficients, using standard deviations.
    if(pf==FALSE){
      penalty_weight=rep(1,p+q-1)
    } else{
      penalty_weight=pw[-i] 
      
    }
    scale_fac<-(p+q-1)/sum(penalty_weight) #this is the scaling factor
    penalty_weight<-scale_fac*penalty_weight
    reg<-glmnet(x=dat[,-i]*scale_fac, y=(dat[,i]),family="gaussian",
                lambda=lambda_g, standardize=F,penalty.factor=penalty_weight)
    
    N[i,-i,]<- rbind((as.matrix(reg$beta))[1:p,], (as.matrix(reg$beta))[ -c(1:p),])*scale_fac
  }
  N  
}


###-------------------------------------###
## neighbour_Ising ####
## neighbourhood selection on the Ising models
## Input: 
##      dat: the data matrix, where the first p columns are Gaussian nodes and the 
##           last q are binary
##      p: number of Gaussian nodes
##      q: number of binary nodes (not really useful)
##      clambda: a vector of constant c
##      ratio: the ratio between the tuning parameters between 
##             the binary nodes and the Gaussian nodes
##      pf: whether to use weights
##      pw: the weights of penalty on each coefficients
## Output: 
##      N: the whole path of estimated neighbourhoods
###-------------------------------------###
## Note:
## This function is modified based on neighbour(),
## and some trunks of code are unnecessary
###-------------------------------------###
neighbour_Ising<-function(dat, p, q, clambda,  pf=FALSE, pw){

  n<-dim(dat)[1]
  lambda_<-clambda*max(c(sqrt(2*log(p)/n),sqrt(2*log(q)/n)))  
  nlambda_g<-length(lambda_g)
  
  N<-replicate(nlambda_g, matrix(0,p+q,p+q))    
  
  for(i in 1:p){
    
    #setting different weights for each coefficients, using standard deviations.
    if(pf==FALSE){
      penalty_weight=rep(1,p+q-1)
    } else{
      penalty_weight=pw[-i] 
    }
    scale_fac<-(p+q-1)/sum(penalty_weight) #this is the scaling factor
    penalty_weight<-scale_fac*penalty_weight
    reg<-glmnet(x=scale_fac*dat[,-i], y=(dat[,i]), family="binomial", 
                lambda=lambda_g, standardize=F,penalty.factor=penalty_weight )    
    N[i,-i,]<- rbind(-(as.matrix(reg$beta))[1:(p-1),], (as.matrix(reg$beta))[p:(p+q-1),])*scale_fac  
  }
  for(i in p+(1:q)){
    if(pf==FALSE){
      penalty_weight=rep(1,p+q-1)
    } else{ 
      penalty_weight=pw[-i] 
    }
    scale_fac<-(p+q-1)/sum(penalty_weight) #this is the scaling factor
    penalty_weight<-scale_fac*penalty_weight
    reg<-glmnet(x=scale_fac*dat[,-i], y=(dat[,i]), family="binomial", 
                lambda=lambda_g, standardize=F,penalty.factor=penalty_weight )
    
    N[i,-i,]<- (as.matrix(reg$beta))*scale_fac  
  }
  N  
}
