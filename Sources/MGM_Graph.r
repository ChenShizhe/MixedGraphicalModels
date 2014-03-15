###--------------------------###
## Graph Generation
## Last updated: Mar.14th 2014
###-------------------------------------###
## This file contains the following functions:
##  1) generator: generating graphs 


###-------------------------------------###
## The graph looks like:
##  O - O - O
##  |   |   |
##  D - D - D
###-------------------------------------###


###-------------------------------------###
## generator ####
## generating graphs 
## Input:
##       p: number of Gaussian nodes
##       q: number of binary nodes
##       lwb: lower bound of edge potentials
##       upb: upper bound of edge potentials
##       PB: T means the P-B networks; F means the G-B networks
## Output:
##       B: the true parameter matrix of Gaussian-Gaussian interactions
##       P: the true parameter matrix of Gaussian-Binary interactions
##       Phi: the true parameter matrix of Binary-Binary interactions
###-------------------------------------###


generator<-function(p, q,lwb, upb, PB=F){
  #B
  B<-matrix(0,p,p)
  for(i in 1:(p-1)){
    B[i,i+1]<-B[i+1,i]<-1
  }    
  B[B!=0]<-(2*rbinom(n=sum(B==1)/2,size=1, prob=0.5) -1)*runif(n=sum(B==1)/2,min=lwb, max=upb)
  B[lower.tri(B)] <- t(B)[lower.tri(B)] 
  diag(B)<-1 
  
  if( PB == F){
  #Normalize B:
  #only do this if the minimum eigen value is negative.
  if(min(eigen(B)$values)<0) {
  B<- B+diag(rep(1,p))*(0.1-min(eigen(B)$values))
  B_norm<-sqrt(diag(1/diag(B)))
  B<-B_norm%*%B%*%B_norm
  }
  } else {
    B<- -abs(B)
  }
  #P
  P<-matrix(0,p,q)
  for(i in 1:p){
    P[i,i]<-1
  }
  P[P!=0]<-(2*rbinom(n=sum(P==1),size=1, prob=0.5) -1)*runif(n=sum(P==1),min=lwb, max=upb)

  
  #Phi    
  Phi<-matrix(0,q,q)
  for(i in 1:(q-1)){
    Phi[i,i+1]<-Phi[i+1,i]<-1
  }    
  Phi[Phi!=0]<-(2*rbinom(n=sum(Phi==1)/2,size=1, prob=0.5) -1)*runif(n=sum(Phi==1)/2,min=lwb, max=upb)
  Phi[lower.tri(Phi)] <- t(Phi)[lower.tri(Phi)] 
  diag(Phi)<-1
  
  
  list(B=B,P=P, Phi=Phi)
}


# ###-------------------------------------###
# ## generator_PB ####
# ## generating graphs for Poisson-binary networks
# ## Input:
# ##       p: number of Gaussian nodes
# ##       q: number of binary nodes
# ##       lwb: lower bound of edge potentials
# ##       upb: upper bound of edge potentials
# ## Output:
# ##       B: the true parameter matrix of Gaussian-Gaussian interactions
# ##       P: the true parameter matrix of Gaussian-Binary interactions
# ##       Phi: the true parameter matrix of Binary-Binary interactions
# ###-------------------------------------###
# 
# generator_PB<-function(p, q, lwb=0.1, upb=0.8){
#   #B
#   B<-matrix(0,p,p)
#   for(i in 1:(p-1)){
#     B[i,i+1]<-B[i+1,i]<-1
#   }    
#   B[B!=0]<-(2*rbinom(n=sum(B==1)/2,size=1, prob=0.5) -1)*runif(n=sum(B==1)/2,min=lwb, max=upb)
#   B[lower.tri(B)] <- t(B)[lower.tri(B)] 
#   B<- -abs(B)
#   
#   
#   #P
#   P<-matrix(0,p,q)
#   for(i in 1:p){
#     P[i,i]<-1
#   }
#   P[P!=0]<-(2*rbinom(n=sum(P==1),size=1, prob=0.5) -1)*runif(n=sum(P==1),min=lwb, max=upb)
#   if(neg==FALSE) {P<-abs(P)}
#   
#   #Phi    
#   Phi<-matrix(0,q,q)
#   for(i in 1:(p-1)){
#     Phi[i,i+1]<-Phi[i+1,i]<-1
#   }    
#   Phi[Phi!=0]<-(2*rbinom(n=sum(Phi==1)/2,size=1, prob=0.5) -1)*
# {runif(n=sum(Phi==1)/2,min=lwb, max=upb)}
# Phi[lower.tri(Phi)] <- t(Phi)[lower.tri(Phi)] 
# 
# list(B=B,P=P, Phi=Phi)
# }

