###-------------------------------------###
## Functions that turn estimates to a graph
## Last updated: Mar.14th 2014
###-------------------------------------###
## This file contains the following functions:
##  1) combine_select: combination using the selection rule (for Poisson-binary)
##  2) combine: a general combination function
##  3) combine_type: recontruction of edges for each edge type


###-------------------------------------###
## Combine_select ####
## combination using the selection rule (for Poisson-binary)
## Input: 
##      N: the estimated neighbourhoods (a (p+q) by (p+q) matrix), used to reconstruct the graph
##      p: number of Poisson nodes
##      q: number of binary nodes
##      parameter: a list containing $a0 (the intercepts)
##                and $coef (the edge potentials)
## Output: 
##      Est: the estimated graph
###-------------------------------------###


combine_select<-function(N, p, q,parameter){
  support<-(N!=0)
  
  a0<-parameter$a0
  coef<-parameter$coef
  bp<-numeric(p)
  
  for(i in 1:p){
    bp[i]<- exp( a0[i]+sum( abs( coef[i,p+(1:q)]) )    )      
  }
  
  Est<-support
  Est[1:p, p+(1:q) ]<- t(Est[p+(1:q), 1:p ]) # set the neighbours from binary nodes by default
  for( i in 1:p){
    if( bp[i] <1 ){
      Est[i,p+(1:q) ]<- support[i,p+(1:q) ]
    } else {
      
    }
    
  }
  Est[p+(1:q), 1:p ] <- t(Est[1:p, p+(1:q) ]) #force the edge to be the same
  Est<-(Est+t(Est))/2
  
  return(Est)
}


###-------------------------------------###
## Combine ####
## a general combination function
## Input: 
##      N: the estimated neighbourhoods (a (p+q) by (p+q) matrix), used to reconstruct the graph
##      p: number of Poisson nodes
##      q: number of binary nodes
##      cc: rule to construct the continous-continous edges ('and', 'or' )
##      cd: rule to construct the continous-discrete edges ('and', 'or', 'gaussian')
##      dd: rule to construct the discrete-discrete edges ('and','or')
## Output: 
##      Est: the estimated graph
###-------------------------------------###
combine<-function(N, p, q, cc,cd,dd, neg=TRUE){
  N_I<-(N!=0)
  E.cc<-N_I[1:p,1:p]
  
  if(cc=="and"){
    E.cc<-(E.cc==t(E.cc))*(N[1:p,1:p]+t(N[1:p,1:p]))/2
  } else if (cc=="or"){
    E.cc<-(N[1:p,1:p]+t(N[1:p,1:p]))/2
  }
  
  E.dd<-N_I[p+1:q,  p+1:q] 
  if(dd=="and"){
    E.dd<-(E.dd==t(E.dd))*(N[p+1:q,  p+1:q]+t(N[p+1:q,  p+1:q]))/2
  } else if (dd=="or"){
    E.dd<-(N[p+1:q,  p+1:q]+t(N[p+1:q,  p+1:q]))/2    
  }
  
  Ec<-N_I[1:p, p+1:q ]
  Ed<-N_I[p+1:q, 1:p ]
  if(cd=="and"){
    E.cd<-(Ec==t(Ed))*(N[1:p, p+1:q ]+t(N[p+1:q, 1:p ]))/2
    E.dc<-t(E.cd)
  } else if (cd=="or"){
    E.cd<-(N[1:p, p+1:q ] +t(N[p+1:q, 1:p ]))/2    
    E.dc<-t(E.cd)
  } else if(cd=="gaussian"){
    E.cd<-N[1:p, p+1:q ]
    E.dc<-t(N[1:p, p+1:q ])
  } else {
    E.cd<-N[1:p, p+1:q ]
    E.dc<-N[p+1:q, 1:p ]
  }
  
  Est<-rbind(cbind(E.cc, E.cd), cbind(E.dc, E.dd))
  return(Est)
}




###-------------------------------------###
## Combine_type ####
## combine_type: recontruction of edges for each edge type
## Input: 
##      N: the estimated neighbourhoods (a (p+q) by (p+q) matrix), used to reconstruct the graph
##      p: number of Poisson nodes
##      q: number of binary nodes
##      cc: rule to construct the continous-continous edges ('and', 'or' )
##      dd: rule to construct the discrete-discrete edges ('and','or')
## Output: 
##      A list of estimated subgraphs: 
##      E.cc (c-c edges), Ec (c-d edges estimated from Gaussian nodes)
##      Ed (c-d edges estimated from binary nodes), E.dd (d-d edges)
###-------------------------------------###
combine_type<-function(N, p, q, cc,dd, neg=TRUE){
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
