###-------------------------------------###
## Data generation
## Last updated: Mar.14th 2014
###-------------------------------------###
## This file contains the following functions:
##  1) sample: draw samples from Gaussian-binary network
##  2) combine: a general combination function
##  3) combine_type: recontruction of edges for each edge type

## Required library:
library(MASS)


###-------------------------------------###
## sampler ####
## generating data from Gaussian-Binary network
## Input: 
##      n: sample size
##      B: a p by p matrix of edge potentials for Gaussian-Gaussian edges
##      P: a p by q matrix of edge potentials for Gaussian-binary edges
##      Phi: a q by q matrix of edge potentials for binary-binary edges
##      seedmultiplier: seed
##      Gibbs.n: Gibbs.n-1 is the number of samples to skip
##      burnin: the length of burn-in period
## Output: 
##      res: a data matrix where the first p columns are Gaussian and the last q are binary (-1,1)
###-------------------------------------###
sampler<-function(n, B, P, Phi, seedmultiplier=12345, Gibbs.n=100,burnin=200){
  p<-dim(B)[1]
  q<-dim(Phi)[1]
  
  Sig<-solve(B)
  M<- t(P)%*%Sig%*%P+Phi
  
  res<-matrix(0,n,p+q)
  
  set.seed(seedmultiplier)
  y<-2*rbinom(q,1,0.5)-1
  for(j in 1:burnin){
    for(i in 1:q){
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




###-------------------------------------###
## sampler_poisson_binary ####
## generating data from Gaussian-Binary network
## Input: 
##      n: sample size
##      B: a p by p matrix of edge potentials for Gaussian-Gaussian edges
##      P: a p by q matrix of edge potentials for Gaussian-binary edges
##      Phi: a q by q matrix of edge potentials for binary-binary edges
##      a0: the alpha_{1s} for Poisson nodes
##      seedmultiplier: seed
##      Gibbs.n: Gibbs.n-1 is the number of samples to skip
##      burnin: the length of burn-in period
## Output: 
##      res: a data matrix where the first p columns are Poisson and the last q are binary (-1,1)
###-------------------------------------###
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

