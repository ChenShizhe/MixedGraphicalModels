###########################
#### Graphical models for mixed data
###########################
## Created: Mar 30th. 2013
## Updated: Feb.22nd 2014
###########################
## Functions:
##  1, generator: a function that generate a Erdos-Renyi graph, with parameters
##     generator.2df: a function that generates a three neighbour graph
##     generator.5df: a function that generates a five neighbour graph
##  2, sampler.old: a function that generates samples from a graph
##     sampler: a function that generates samples from discrete samples and the graph
##  3, neighbour: estimating the neighbour for all node
##  4, combine: reconstruct graph from neighbours
##  5, perf: assessing performance of estimated graph.
##  6, is.same: testing whether the sample is degenerated
##  7, wholegraph: draw a heatmap for the whole graph
############################
##
## Required package: MASS, glmnet
##

###################################################
### 1,Graph generator:
###
### Input: p, q, weight=c(0.01,0.01,0.01); lwb: 0.5; upb: 1; neg=True
### Output: B= A p by p matrix, p.d.;
###         P= A p by q matrix;
###         Phi= A q by q matrix, symmetric;

# The huge matrix = rbind(cbind(B,P),cbind(t(P),Phi))
library(glasso)
generator<-function(p, q, weight, lwb, upb, neg=TRUE){
  #B.cc
  B<-matrix(0,p,p)
  B[upper.tri(B)]<-rbinom(n= (p*(p-1)/2), size=1, prob=weight[1])
  B[B!=0]<-(2*rbinom(n=sum(B==1),size=1, prob=0.5) -1)*runif(n=sum(B==1),min=lwb, max=upb)
  B[lower.tri(B)] <- t(B)[lower.tri(B)] 
  diag(B)<-1
  if(neg==FALSE) {B<-abs(B)}
  #Normalize B:
  if(min(eigen(B)$values)<0) {diag(B)<- diag(B)+(0.1-min(eigen(B)$values))}
  B<-B/diag(B)[1]

  #P
  P<-rbinom(n= p*q, size=1, prob=weight[2])
  dim(P)<-c(p,q)
  P[P!=0]<-(2*rbinom(n=sum(P==1),size=1, prob=0.5) -1)*runif(n=sum(P==1),min=lwb, max=upb)
  if(neg==FALSE) {P<-abs(P)}
  
  #Phi
  Phi<-matrix(0,q,q)
  Phi[upper.tri(Phi)]<-rbinom(n= (q*(q-1)/2), size=1, prob=weight[3])
  Phi[Phi!=0]<-(2*rbinom(n=sum(Phi==1),size=1, prob=0.5) -1)*runif(n=sum(Phi==1),min=lwb, max=upb)
  Phi[lower.tri(Phi)] <- t(Phi)[lower.tri(Phi)] 
  diag(Phi)<-1
  if(neg==FALSE) {Phi<-abs(Phi)}
  
  list(B=B,P=P, Phi=Phi)
}

##############################
## Anothe graph generator
## The new graph looks like:
##  O - O - O
##  |   |   |
##  D - D - D
##############################

generator.3df<-function(p, lwb, upb, neg=TRUE){
  #B.cc
  B<-matrix(0,p,p)
  for(i in 1:(p-1)){
    B[i,i+1]<-B[i+1,i]<-1
  }    
  B[B!=0]<-(2*rbinom(n=sum(B==1)/2,size=1, prob=0.5) -1)*runif(n=sum(B==1)/2,min=lwb, max=upb)
  B[lower.tri(B)] <- t(B)[lower.tri(B)] 
  diag(B)<-1
  if(neg==FALSE) {B<-abs(B)}
  #Normalize B:
  #only do this if the minimum eigen value is negative.
  if(min(eigen(B)$values)<0) {diag(B)<- diag(B)+(0.1-min(eigen(B)$values))}
  B<-B/diag(B)[1]
  
  #P
  P<-diag(p)
  P[P!=0]<-(2*rbinom(n=sum(P==1),size=1, prob=0.5) -1)*runif(n=sum(P==1),min=lwb, max=upb)
  if(neg==FALSE) {P<-abs(P)}
  
  #Phi    
  Phi<-matrix(0,p,p)
  for(i in 1:(p-1)){
      Phi[i,i+1]<-Phi[i+1,i]<-1
    }    
  Phi[Phi!=0]<-(2*rbinom(n=sum(Phi==1)/2,size=1, prob=0.5) -1)*runif(n=sum(Phi==1)/2,min=lwb, max=upb)
  Phi[lower.tri(Phi)] <- t(Phi)[lower.tri(Phi)] 
  diag(Phi)<-1
  if(neg==FALSE) {B<-abs(Phi)}
  list(B=B,P=P, Phi=Phi)
}

# 
# generator_PB<-function(p, lwb=c(0.1, 0.4), upb=c(0.2,0.8), neg=TRUE){
#   #B.cc
#   B<-matrix(0,p,p)
#   for(i in 1:(p-1)){
#     B[i,i+1]<-B[i+1,i]<-1
#   }    
#   B[B!=0]<-(2*rbinom(n=sum(B==1)/2,size=1, prob=0.5) -1)*
#     {runif(n=sum(B==1)/2,min=lwb[1], max=upb[1])}
#   B[lower.tri(B)] <- t(B)[lower.tri(B)] 
#   diag(B)<-1
#   if(neg==FALSE) {B<-abs(B)}
#   
#   
#   #P
#   P<-diag(p)
#   P[P!=0]<-(2*rbinom(n=sum(P==1),size=1, prob=0.5) -1)*
#   {runif(n=sum(P==1),min=lwb[1], max=upb[1])+rbinom(n=sum(P==1),size=1,prob=0.5)*runif(n=sum(P==1),min=lwb[2], max=upb[2])   }
#   if(neg==FALSE) {P<-abs(P)}
#   
#   #Phi    
#   Phi<-matrix(0,p,p)
#   for(i in 1:(p-1)){
#     Phi[i,i+1]<-Phi[i+1,i]<-1
#   }    
#   Phi[Phi!=0]<-(2*rbinom(n=sum(Phi==1)/2,size=1, prob=0.5) -1)*
#   {runif(n=sum(Phi==1)/2,min=lwb[1], max=upb[1])}
#   Phi[lower.tri(Phi)] <- t(Phi)[lower.tri(Phi)] 
#   diag(Phi)<-1
#   if(neg==FALSE) {B<-abs(Phi)}
#   list(B=B,P=P, Phi=Phi)
# }

##



generator_PB<-function(p, lwb=0.1, upb=0.8, neg=TRUE){
  #B.cc
  B<-matrix(0,p,p)
  for(i in 1:(p-1)){
    B[i,i+1]<-B[i+1,i]<-1
  }    
  B[B!=0]<-(2*rbinom(n=sum(B==1)/2,size=1, prob=0.5) -1)*runif(n=sum(B==1)/2,min=lwb, max=upb)
  B[lower.tri(B)] <- t(B)[lower.tri(B)] 
  B<- -abs(B)
  #P
P<-diag(rep(1,p))
# for(i in 1:(p-1)){
#   P[i+1,i]<-1
# }    

P[P!=0]<-(2*rbinom(n=sum(P==1),size=1, prob=0.5) -1)*runif(n=sum(P==1),min=lwb, max=upb)
if(neg==FALSE) {P<-abs(P)}

#Phi    
  Phi<-matrix(0,p,p)
  for(i in 1:(p-1)){
    Phi[i,i+1]<-Phi[i+1,i]<-1
  }    
  Phi[Phi!=0]<-(2*rbinom(n=sum(Phi==1)/2,size=1, prob=0.5) -1)*
  {runif(n=sum(Phi==1)/2,min=lwb, max=upb)}
  Phi[lower.tri(Phi)] <- t(Phi)[lower.tri(Phi)] 

list(B=B,P=P, Phi=Phi)
}
##############################
## Anothe graph generator with more edges
##############################
# 
# generator.5df<-function(p, lwb, upb, neg=TRUE){
#   #B.cc
#   B<-matrix(0,p,p)
#   for(i in 1:(p-1)){
#     B[i,p+1-i]<-B[i,p-i]<-B[i+1,p-i+1]<-B[i,i+1]<-B[i+1,i]<-1
#   }    
#   B[B!=0]<-(2*rbinom(n=sum(B>0),size=1, prob=0.5) -1)*runif(n=sum(B>0),min=lwb, max=upb)
#   B[lower.tri(B)] <- t(B)[lower.tri(B)] 
#   diag(B)<-1
#   if(neg==FALSE) {B<-abs(B)}
#   #Normalize B:
#   if(min(eigen(B)$values)<0) {diag(B)<- diag(B)+(0.1-min(eigen(B)$values))}
#   B<-B/diag(B)[1]
#   
#   #P
#   P<-diag(p)
#   P[P!=0]<-(2*rbinom(n=sum(P==1),size=1, prob=0.5) -1)*runif(n=sum(P==1),min=lwb, max=upb)
#   if(neg==FALSE) {P<-abs(P)}
#   
#   #Phi    
#   Phi<-matrix(0,p,p)
#   for(i in 1:(p-1)){
#     Phi[i,p+1-i]<-Phi[i,p-i]<-Phi[i+1,p-i+1]<-Phi[i,i+1]<-Phi[i+1,i]<-1
#   }    
#   Phi[Phi!=0]<-(2*rbinom(n=sum(Phi>0),size=1, prob=0.5) -1)*runif(n=sum(Phi>0),min=lwb, max=upb)
#   Phi[lower.tri(Phi)] <- t(Phi)[lower.tri(Phi)] 
#   diag(Phi)<-1
#   if(neg==FALSE) {B<-abs(Phi)}
#   
#   list(B=B,P=P, Phi=Phi)
# }
################################
## 2,Sampler.old (old version): 
## Input: sample size: n; Parameter matrices, B, P and Phi;
##        seedmultiplier: 12345; Magic number: Gibbs.n
## Output: A matrix of samples, with each row being one observation

##########
## Required library:
library(MASS)

sampler.old<-function(n, B, P, Phi, seedmultiplier=12345, Gibbs.n=100,burnin=200){
  p<-dim(B)[1]
  q<-dim(Phi)[1]

  Sig<-solve(B)
  M<- t(P)%*%Sig%*%P+Phi

  res<-matrix(0,n,p+q)

  set.seed(seedmultiplier)
  y<-2*rbinom(q,1,0.5)-1
  for(j in 1:burnin){
    for(i in 1:q){
      ###
#       y[i]<-1
#       odd.pos<-exp(t(y)%*%M%*%y/2)
#       y[i]<- -1
#       odd.neg<-exp(t(y)%*%M%*%y/2)
#       #alternatively
      fac<-t(y[-i])%*%M[-i,i]
      odd.pos<-exp(1*fac )
      odd.neg<-exp(-1*fac )
      y[i]<-2*rbinom(1, 1, prob=odd.pos/(odd.neg+odd.pos) )-1  
    }      
  }
  for(k in 1:n){    
    #Draw samples from p(y)
    for(j in 1:Gibbs.n){
      for(i in 1:q){
        fac<-t(y[-i])%*%M[-i,i]
        odd.pos<-exp(1*fac )
        odd.neg<-exp(-1*fac )
        y[i]<-2*rbinom(1, 1, prob=odd.pos/(odd.neg+odd.pos) )-1  
      }  
    }  
   res[k,(p+1):(p+q)]<-y  
  }
      #draw samples from p(x|y)
    WM<-Sig%*%P
    for(k in 1:n){
    res[k,1:p]<-mvrnorm(n=1,mu=WM%*%res[k,(p+1):(p+q)] , Sigma=Sig)
    } 
  res
}

#take around 80 seconds to run p=80, q=40, n=200


################################
## Sampler (new version): 
## Input: A matrix of discrete samples: Y; Parameter matrices, B, P and Phi;
##        seedmultiplier: 12345;
## Output: A matrix of samples, with each row being one observation 

##########
## Required library:
library(MASS)

sampler<-function(Y, B, P, Phi, seedmultiplier=12345){
  p<-dim(B)[1]
  q<-dim(Phi)[1]
  
  Sig<-solve(B)
  
  res<-matrix(0,dim(Y)[1],p+q)
  res[, p+(1:q)]<-Y
  set.seed(seedmultiplier)
  
  for(k in 1:dim(Y)[1]){
  #draw samples from p(x|y)
    y<-Y[k,]
    x<-mvrnorm(n=1,mu=Sig%*%P%*%y, Sigma=Sig)
    res[k,1:p]<-x  
  }
  res
}

#############################
# sampler_poisson_binary
# Input:
#       the same as the previous one
# Output:
#       a matrix of samples

# the offset a0 is set so that the information from a Poisson node is not too strong
sampler_poisson_binary<-function(n, B, P, Phi, a0= NULL, seedmultiplier=12345, Gibbs.n=100,burnin=200){
  p<-dim(B)[1]
  q<-dim(Phi)[1]
  
  res<-matrix(0,n,p+q)
  
  set.seed(seedmultiplier)
  x<-rpois(p,lambda=1)
  y<-2*rbinom(q,1,0.5)-1
  
  for(j in 1:burnin){
    
    for(i in 1:q){
      fac<- t(x)%*%P[,i]+t(y[-i])%*%Phi[-i,i]  
      odd.pos<-exp(1*fac )
      odd.neg<-exp(-1*fac )
      
      y[i]<-2*rbinom(1, 1, prob=odd.pos/(odd.neg+odd.pos) )-1  
    }      
    for(i in 1:p){
      fac<- t(x[-i])%*%B[-i,i]+ t(y)%*%P[i,] +a0[i]
      
      x[i]<-rpois(1,lambda=exp(fac)) 
      
    }
  }
  for(k in 1:n){    
    #Draw samples from p(y)
    for(j in 1:Gibbs.n){
      for(i in 1:q){
        fac<- t(x)%*%P[,i]+t(y[-i])%*%Phi[-i,i]  
        odd.pos<-exp(1*fac )
        odd.neg<-exp(-1*fac )
        
        y[i]<-2*rbinom(1, 1, prob=odd.pos/(odd.neg+odd.pos) )-1  
      }      
      for(i in 1:p){
        fac<- t(x[-i])%*%B[-i,i]+ t(y)%*%P[i,]+a0[i] 
        
        x[i]<-rpois(1,lambda=exp(fac)) 
        
      }
    }
    res[k,]<-c(x,y)  
  }
  res
}
##############################################
### 3,Neighbour: 
### Input: dat; p; q; coefficient for lambda: these coefficients are put in a vector form
### Output: An array , with the third dimension being indicator, each row being a vector indicates the neighbours
library(glmnet)

neighbour.value<-function(dat, p, q, clambda.g, ratio=1, pf=FALSE, pw,cf=T){
  ##determine lambda.g and lambda.b
  n<-dim(dat)[1]
  #lambda.g and lambda.b could be different.
  lambda.g<-clambda.g*max(c(sqrt(2*log(p)/n),sqrt(2*log(q)/n)))  
  n.lambda.g<-length(lambda.g)
  
  lambda.b<- ratio*lambda.g
  N<-replicate(n.lambda.g, matrix(0,p+q,p+q))    
  A<-nlg<-bic<-matrix(0,p+q, n.lambda.g)
  # the third dimension: (1,1) ,(2,1) ,..., (n.lambda.g,n.lambda.b)
  for(i in 1:p){
    
    #setting different weights for each coefficients, using standard deviations.
    if(pf==FALSE){
      penalty.weight=rep(1,p+q-1)
    } else{
      if(cf==T){
        penalty.weight=pw[-i]*pw[i]
      }else{ penalty.weight=pw[-i] }
    }
    scale.fac<-(p+q-1)/sum(penalty.weight) #this is the scaling factor
    penalty.weight<-scale.fac*penalty.weight
    reg<-glmnet(x=scale.fac*dat[,-i], y=(dat[,i]), family="gaussian", 
                lambda=lambda.g, standardize=F,penalty.factor=penalty.weight )
    
    N[i,-i,]<- rbind(-(as.matrix(reg$beta))[1:(p-1),], (as.matrix(reg$beta))[p:(p+q-1),])*scale.fac
    A[i,]<-as.matrix(reg$a0)
    
    for(j in 1: n.lambda.g ){
      nlg[i,j]<-N2Loglklh(beta=as.matrix(reg$beta)[,j], a0=as.matrix(reg$a0)[j],
                          y=dat[,i], x=scale.fac*dat[,-i],logistic=F)
      bic[i,j]<-BIC(nlg[i,j],n,as.matrix(reg$beta)[,j])
    }
    
    
  }
  for(i in p+(1:q)){
    #setting different weights for each coefficients, using standard deviations.
    if(pf==FALSE){
      penalty.weight=rep(1,p+q-1)
    } else{
      
      if(cf==T){
        penalty.weight=pw[-i]*pw[i]
      }else{ penalty.weight=pw[-i] }
    }
    scale.fac<-(p+q-1)/sum(penalty.weight) #this is the scaling factor
    # because glmnet will automatically scale the penalty
    # we scale the covariates to obtain the same penalty 
    penalty.weight<-scale.fac*penalty.weight
    reg<-glmnet(x=scale.fac*dat[,-i], y=(dat[,i]>mean(dat[,i])), family="binomial", 
                lambda=lambda.b, standardize=F,penalty.factor=penalty.weight )
    
    N[i,-i,]<- (as.matrix(reg$beta))*scale.fac  
    A[i,]<-as.matrix(reg$a0)
    for(j in 1: n.lambda.g ){
      nlg[i,j]<-N2Loglklh(beta=as.matrix(reg$beta)[,j], a0=as.matrix(reg$a0)[j],
                          y=dat[,i], x=scale.fac*dat[,-i],logistic=T)
      bic[i,j]<-BIC(nlg[i,j],n,as.matrix(reg$beta)[,j])
      
    }
    
  }
  return(list(N=N,Intercept=A,nlg=nlg,bic=bic))  
}


#####################
## Poisson-Binary networks



neighbour_PB<-function(dat, p, q, clambda.g, ratio=FALSE, pf=FALSE, pw,cf=T,maxit){
  ##determine lambda.g and lambda.b
  n<-dim(dat)[1]
  #lambda.g and lambda.b could be different.
  lambda.g<-clambda.g*max(c(sqrt(2*log(p)/n),sqrt(2*log(q)/n)))  
  n.lambda.g<-length(lambda.g)
  
  lambda.b<- ratio*lambda.g
  N<-replicate(n.lambda.g, matrix(NA,p+q,p+q))    
  A<-nlg<-bic<-matrix(NA,p+q, n.lambda.g)

  # the third dimension: (1,1) ,(2,1) ,..., (n.lambda.g,n.lambda.b)
  for(i in 1:p){
    
    #setting different weights for each coefficients, using standard deviations.
    if(pf==FALSE){
      penalty.weight=rep(1,p+q-1)
    } else{
      if(cf==T){
        penalty.weight=pw[-i]*pw[i]
      }else{ penalty.weight=pw[-i] }
    }
    scale.fac<-(p+q-1)/sum(penalty.weight) #this is the scaling factor
    penalty.weight<-scale.fac*penalty.weight
    reg<-glmnet(x=scale.fac*dat[,-i], y=(dat[,i]), family="poisson", 
                lambda=lambda.g, standardize=F,penalty.factor=penalty.weight,maxit=maxit )
    
      # it gives the whole solution path
      # when solutions are not available for all lambda
      limi<-dim(as.matrix(reg$beta))[2]
    
     
    N[i,-i,1:limi]<-(as.matrix(reg$beta))*scale.fac  
    A[i,1:limi]<-as.matrix(reg$a0)
    
    for(j in 1: limi ){
      nlg[i,j]<-N2Loglklh(beta=as.matrix(reg$beta)[,j], a0=as.matrix(reg$a0)[j],
                          y=dat[,i], x=scale.fac*dat[,-i],logistic="P")
      bic[i,j]<-BIC(nlg[i,j],n,as.matrix(reg$beta)[,j])
    }
    
    
  }
  for(i in p+(1:q)){
    #setting different weights for each coefficients, using standard deviations.
    if(pf==FALSE){
      penalty.weight=rep(1,p+q-1)
    } else{
      
      if(cf==T){
        penalty.weight=pw[-i]*pw[i]
      }else{ penalty.weight=pw[-i] }
    }
    scale.fac<-(p+q-1)/sum(penalty.weight) #this is the scaling factor
    # because glmnet will automatically scale the penalty
    # we scale the covariates to obtain the same penalty 
    penalty.weight<-scale.fac*penalty.weight
    reg<-glmnet(x=scale.fac*dat[,-i], y=(dat[,i]>mean(dat[,i])), family="binomial", 
                lambda=lambda.b, standardize=F,penalty.factor=penalty.weight,maxit=maxit )
    # it gives the whole solution path
    # when solutions are not available for all lambda
    limi<-dim(as.matrix(reg$beta))[2]
    
    
    N[i,-i,1:limi]<-(as.matrix(reg$beta))*scale.fac  
    A[i,1:limi]<-as.matrix(reg$a0)
    for(j in 1: limi ){
      nlg[i,j]<-N2Loglklh(beta=as.matrix(reg$beta)[,j], a0=as.matrix(reg$a0)[j],
                          y=dat[,i], x=scale.fac*dat[,-i],logistic=T)
      bic[i,j]<-BIC(nlg[i,j],n,as.matrix(reg$beta)[,j])
      
    }
    
  }
  return(list(N=N,Intercept=A,nlg=nlg,bic=bic))  
}



###### M-B#####
#cf=fommon factor

MB<-function(dat, p, q, clambda.g,  pf=FALSE, pw, cf=T){
  ##determine lambda.g and lambda.b
  n<-dim(dat)[1]
  #lambda.g and lambda.b could be different.
  lambda.g<-clambda.g*max(c(sqrt(2*log(p)/n),sqrt(2*log(q)/n)))  
  n.lambda.g<-length(lambda.g)

  N<-replicate(n.lambda.g, matrix(0,p+q,p+q))    
  
  # the third dimension: (1,1) ,(2,1) ,..., (n.lambda.g,n.lambda.b)
  for(i in 1:p){
    
    #setting different weights for each coefficients, using standard deviations.
    if(pf==FALSE){
      penalty.weight=rep(1,p+q-1)
    } else{
      if(cf==T){
        penalty.weight=pw[-i]*pw[i]
      }else{ penalty.weight=pw[-i] }
      
      }
    #The penalty factors is scaled to p+q-1
    scale.fac<-(p+q-1)/sum(penalty.weight) #this is the scaling factor
    penalty.weight<-scale.fac*penalty.weight
    scale.sq<-sqrt(scale.fac)
    reg<-glmnet(x=scale.fac*dat[,-i], y=(dat[,i]), family="gaussian", 
                lambda=lambda.g, standardize=F,penalty.factor=penalty.weight)
    N[i,-i,]<- rbind(-(as.matrix(reg$beta))[1:(p-1),], (as.matrix(reg$beta))[p:(p+q-1),])*scale.fac  
  }
  for(i in p+(1:q)){
    #setting different weights for each coefficients, using standard deviations.
    if(pf==FALSE){
      penalty.weight=rep(1,p+q-1)
    } else{
      
      if(cf==T){
        penalty.weight=pw[-i]*pw[i]
      }else{ penalty.weight=pw[-i] }
      
    }
    scale.fac<-(p+q-1)/sum(penalty.weight) #this is the scaling factor
    penalty.weight<-scale.fac*penalty.weight
    scale.sq<-sqrt(scale.fac)
    reg<-glmnet(x=dat[,-i]*scale.fac, y=(dat[,i]),family="gaussian",
                lambda=lambda.g, standardize=F,penalty.factor=penalty.weight)
    
    N[i,-i,]<- rbind((as.matrix(reg$beta))[1:p,], (as.matrix(reg$beta))[ -c(1:p),])*scale.fac
  }
  N  
}


Ising<-function(dat, p, q, clambda.g,  pf=FALSE, pw, cf=T){
  ##determine lambda.g and lambda.b
  n<-dim(dat)[1]
  #lambda.g and lambda.b could be different.
  lambda.g<-clambda.g*max(c(sqrt(2*log(p)/n),sqrt(2*log(q)/n)))  
  n.lambda.g<-length(lambda.g)
  
  N<-replicate(n.lambda.g, matrix(0,p+q,p+q))    
  
  # the third dimension: (1,1) ,(2,1) ,..., (n.lambda.g,n.lambda.b)
  for(i in 1:p){
    
    #setting different weights for each coefficients, using standard deviations.
    if(pf==FALSE){
      penalty.weight=rep(1,p+q-1)
    } else{
      
      if(cf==T){
        penalty.weight=pw[-i]*pw[i]
      }else{ penalty.weight=pw[-i] }
    }
    scale.fac<-(p+q-1)/sum(penalty.weight) #this is the scaling factor
    penalty.weight<-scale.fac*penalty.weight
    scale.sq<-sqrt(scale.fac)
    reg<-glmnet(x=scale.fac*dat[,-i], y=(dat[,i]), family="binomial", 
                lambda=lambda.g, standardize=F,penalty.factor=penalty.weight )    
    N[i,-i,]<- rbind(-(as.matrix(reg$beta))[1:(p-1),], (as.matrix(reg$beta))[p:(p+q-1),])*scale.fac  
  }
  for(i in p+(1:q)){
    #setting different weights for each coefficients, using standard deviations.
    if(pf==FALSE){
      penalty.weight=rep(1,p+q-1)
    } else{
      
      if(cf==T){
        penalty.weight=pw[-i]*pw[i]
      }else{ penalty.weight=pw[-i] }
    }
    scale.fac<-(p+q-1)/sum(penalty.weight) #this is the scaling factor
    penalty.weight<-scale.fac*penalty.weight
    scale.sq<-sqrt(scale.fac)
    reg<-glmnet(x=scale.fac*dat[,-i], y=(dat[,i]), family="binomial", 
                lambda=lambda.g, standardize=F,penalty.factor=penalty.weight )
    
    N[i,-i,]<- (as.matrix(reg$beta))*scale.fac  
  }
  N  
}

#############################################
#### 4,combine results
#### Input: A neighbourhood matrix, with each row being a vector indicates the neighbours
####        cc = "and", "or"
####        cd = "and", "or", "gaussian"
####        dd = "and", "or"
####        neg= TRUE
####        p, q
#### Output: E.cc, E.cd and Phi.dd


combine<-function(N, p, q, cc,cd,dd, neg=TRUE){
  N<-sign(N)
  E.cc<-N[1:p,1:p]
  
  if(cc=="and"){
    E.cc<-(E.cc==t(E.cc))*E.cc
  } else if (cc=="or"){
    E.cc<-E.cc+t(E.cc)    
    E.cc<-sign(E.cc)
  }
  
  E.dd<-N[p+1:q,  p+1:q] 
  if(dd=="and"){
    E.dd<-(E.dd==t(E.dd))*E.dd
  } else if (dd=="or"){
    E.dd<-E.dd+t(E.dd)    
    E.dd<-sign(E.dd)
  }
  
  
  Ec<-N[1:p, p+1:q ]
  Ed<-N[p+1:q, 1:p ]
  if(cd=="and"){
    E.cd<-(Ec==t(Ed))*Ec
  } else if (cd=="or"){
    E.cd<-Ec+t(Ed)    
    E.cd<-sign(E.cd)
  } else if(cd=="gaussian"){
    E.cd<-Ec
  }
  b<-c(E.cc[upper.tri(E.cc)])
  pC<-c(E.cd)  
  phi<-c(E.dd[upper.tri(E.dd)] )  
  pD<-pC
 list(b=b,pC=pC,pD=pD,phi=phi )
}

discrete<-function(N, dd, neg=TRUE){
  
  E.dd<-N
  if(dd=="and"){
    E.dd<-(E.dd==t(E.dd))*E.dd
  } else if (dd=="or"){
    E.dd<-E.dd+t(E.dd)    
    E.dd<-sign(E.dd)
  }
  list(E.dd=E.dd )
}
# used in probability
combine.type<-function(N, p, q, cc,dd, neg=TRUE){
  N<-sign(N)
  E.cc<-N[1:p,1:p]
  
  if(cc=="and"){
    E.cc<-(E.cc==t(E.cc))*E.cc
  } else if (cc=="or"){
    E.cc<-E.cc+t(E.cc)    
    E.cc<-sign(E.cc)
  }
  
  E.dd<-N[p+1:q,  p+1:q] 
  if(dd=="and"){
    E.dd<-(E.dd==t(E.dd))*E.dd
  } else if (dd=="or"){
    E.dd<-E.dd+t(E.dd)    
    E.dd<-sign(E.dd)
  }
  
  Ec<-N[1:p, p+1:q ]
  Ed<-t(N[p+1:q, 1:p ])
  list(E.cc=E.cc, Ec=Ec,Ed=Ed, E.dd=E.dd )
}
combine.value<-function(N, p, q,corr=F,pw){
  if(corr==T){
    left<-diag(1/pw)
    right<-diag(1/pw)
    N<-left%*%N%*%right
  }
  
  E.cc<-N[1:p,1:p]
  
  Ec<-N[1:p, p+1:q ]
  Ed<-N[p+1:q, 1:p ]
  
  E.dd<-N[p+1:q,  p+1:q] 
  b<-c({E.cc+t(E.cc)}[upper.tri(E.cc)])/2
  p<-c(Ec+t(Ed))/2
  phi<-c({E.dd+t(E.dd)}[upper.tri(E.dd)] )/2
  list(b=b,p=p,phi=phi)
}

combine.node<-function(N, p, q, corr=F, pw){
  
  if(corr==T){
  left<-diag(1/pw)
  right<-diag(1/pw)
  N<-left%*%N%*%right
  }
  
  E.cc<-N[1:p,1:p] 
  Ec<-N[1:p, p+1:q ]
  Ed<-N[p+1:q, 1:p ]
  E.dd<-N[p+1:q,  p+1:q] 
  
  
  b<-c({E.cc+t(E.cc)}[upper.tri(E.cc)])
  pC<-c(Ec)

  
  phi<-c({E.dd+t(E.dd)}[upper.tri(E.dd)] )

  pD<-c(t(Ed))
  list(b=b,pC=pC,pD=pD,phi=phi)
}
#########################
### 6, Testing whether the sample is degenerate.
### Input: a data set
### Output: T or F
is.same<-function(dat){
  std<-apply(dat,2,sd) 
  return( (min(std)==0))  
}


###################################################
#### 5,Assessing performance
#### Input: a list of estimated (edges) matrix Est, [[1]] [[2]] [[3]]
####        a list of True parameter matrix: E
####        signed edges: neg=True
####        Option: op= "graph", "edges".
#### Output: op=="graph", 0, 1 for each part of edges
####         op="edges", total (estimated) edges: total, and true edges
# 
# 
# perf<-function(Est, E, neg=TRUE, op="graph"  ){
#   
#   if(op=="graph"){
#     diag(Est[[1]])<-1
#     diag(Est[[3]])<-1
#     
#     if(neg==TRUE){ #number of disagreements
#       cc<-sum(sign(Est[[1]])!=sign(E[[1]]))==0 
#       cd<-sum(sign(Est[[2]])!=sign(E[[2]]))==0
#       dd<-sum(sign(Est[[3]])!=sign(E[[3]]))==0  
#     } else if(neg==FALSE){
#       cc<-sum( (Est[[1]]!=0) - (E[[1]]!=0) )==0 
#       cd<-sum( (Est[[2]]!=0) - (E[[2]]!=0) )==0
#       dd<-sum( (Est[[3]]!=0) - (E[[3]]!=0) )==0      
#     }
#     return(c(cc,cd,dd))
#   } else if(op=="edges") {
#     total.cc<- sum(Est[[1]]!=0)/2
#     total.cd<-sum(Est[[2]]!=0)
#     total.dd<-sum(Est[[3]]!=0)/2
#     if(neg==TRUE){
#       cc<-sum( ( sign(Est[[1]]) ==sign(E[[1]]) )*( Est[[1]]!=0 ) )  /2
#       cd<-sum((sign(Est[[2]])==sign(E[[2]]))*(Est[[2]]!=0))
#       dd<-sum((sign(Est[[3]])==sign(E[[3]]))*(Est[[3]]!=0))/2
#     }else if(neg==FALSE ){  
#       cc<-sum( (Est[[1]]!=0) * (E[[1]]!=0)  ) /2
#       cd<-sum( (Est[[2]]!=0) * (E[[2]]!=0)  )/2
#       dd<-sum( (Est[[3]]!=0) * (E[[3]]!=0)  )/2
#     } 
#     total<-c(total.cc,total.cd, total.dd)
#     true<-c(cc,cd,dd)
#     return(rbind(total,true))
#   } 
# }
# perf.disc<-function(Est,E){
#   total.dd<-(sum(Est!=0)-dim(Est)[1])/2
#   dd<-(sum((sign(Est)==sign(E))*(Est!=0))-dim(Est)[1])/2
#   c(dd,total.dd)
# }

# #########################
# ### 7, Generating a heatmap for graphs
# ### Input: a huge parameter matrix, p+q by p+q
# ### Output: a picture.
# 
# wholegraph<-function(graph,col=c("white","black")){
#   m<-rbind(cbind(graph[[1]],graph[[2]]),cbind(t(graph[[2]]), graph[[3]] ) )
#   diag(m)<-0
#   m <- t(m)[,nrow(m):1]
#   image(m,col=col,xaxt='n',yaxt='n')
# }


############################
### 8, Create a huge vector of the absolute value of matrix

# #library(Kendall)
# huge<-function(graph){  
#   b<-graph[[1]]
#   b<-c(b[upper.tri(b)])
#   p<-graph[[2]]
#   p<-c(p)
#   phi<-graph[[3]]
#   phi<-c( phi[upper.tri(phi)])
#   list(b=b,p=p,phi=phi)
# }
# 

##New measure:
## Concordance and disconcordance for non-zero true edges;
## times the discovery rate!
# concord<-function(x,y){
#   indi<-(1:length(y))[y!=0]
#   c(cor(x[indi],y[indi],method="kendall"), sum(x[indi]!=0))
# }
