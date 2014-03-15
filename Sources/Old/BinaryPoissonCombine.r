###-------------------------------------###
## Functions related with the Poisson-Binary networks
## Created: Feb.16th 2014
## Updated: Feb.25th 2014
###-------------------------------------###
### This file contains the following functions:
## 1) generator_PB: a function that generates the graph for Poisson-Binary data
## 2) sampler_poisson_binary
## 3) neighbour_PB
## 4) combine_smart
## 5) edge_eval
## 6) edge_aver
## 7) combine.node
## 8) combine


###
  # The input parameter should contain both the intercept and the estimated parameters
  
  # parameter$a0
  # parameter$coef
combine_smart<-function(estimate, p, q,parameter){
  support<-(estimate!=0)
  
  a0<-parameter$a0
  coef<-parameter$coef
    bp<-numeric(p)
  
  
  for(i in 1:p){
      bp[i]<- exp( a0[i]+sum( abs( coef[i,p+(1:q)]) )    )      
  }
  
N<-support
 N[1:p, p+(1:q) ]<- t(N[p+(1:q), 1:p ]) # set the neighbours from binary nodes by default
  for( i in 1:p){
    if( bp[i] <1 ){
      N[i,p+(1:q) ]<- support[i,p+(1:q) ]
      } else {
          
    }
      
  }
  N[p+(1:q), 1:p ] <- t(N[1:p, p+(1:q) ]) #force the edge to be the same
  N<-(N+t(N))/2

  return(N)
}




combine_node<-function(N, p, q, corr=F, pw){
  
  if(corr==T){
    left<-diag(1/pw)
    right<-diag(1/pw)
    N<-left%*%N%*%right
  }
  
  E.cc<-N[1:p,1:p] 
  Ec<-N[1:p, p+1:q ]
  Ed<-N[p+1:q, 1:p ]
  E.dd<-N[p+1:q,  p+1:q] 
  
  
  E.cc<-{E.cc+t(E.cc)}/2
  E.dd<-{E.dd+t(E.dd)}/2
  
  Est<-rbind(cbind(E.cc, Ec), cbind(Ed, E.dd))
  return(Est)
}



combine_gen<-function(N, p, q, cc,cd,dd, neg=TRUE){
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
