###########################################  
## Evaluation function
## Created: Feb.14th 2014
## Updated: Feb.14th 2014
###-------------------------------------###
### This file contains the following functions:
### 1) Eva.total: a function that calculate the total number of estimated edges
### 2) Eva.count: a funtion that evaluate the TP and FP for one huge matrix
###


# Eva.total
# Input: 
#       B: the true parameter matrix of Gaussian-Gaussian interactions
#       P: the true parameter matrix of Gaussian-Binary interactions
#       Phi: the true parameter matrix of Binary-Binary interactions
#       Est: the estimated parameters
#       reverse: whether the binary node preceeds the continuous
# Output:
#       TE: number of estimated edges
# Note:
#       PC and PD are Rho matrices from the continuous nodes and discrete nodes, resp.
Eva.total<-function(B, P, Phi, Est, reverse=F){
  p<-dim(B)[1]
  q<-dim(Phi)[1]
  
  if (reverse==T){
    Est.Phi<-Est[1:p,1:p]
    Est.PC<-Est[(p+1):(p+q),1:p ]
    Est.PD<-t(Est[1:p, (p+1):(p+q)])
    Est.B<-Est[(p+1):(p+q), (p+1):(p+q)]
  }else{
    Est.B<-Est[1:p,1:p]
    Est.PC<-Est[1:p, (p+1):(p+q)]
    Est.PD<-t(Est[(p+1):(p+q),1:p ])
    Est.Phi<-Est[(p+1):(p+q), (p+1):(p+q)]
    
  }  
  
  diag(Est.B)<-0
  diag(Est.Phi)<-0
  
  TE<-numeric(4)
  TE[1]<-sum( (Est.B!=0),na.rm=T )/2
  TE[2]<-sum( (Est.PC!=0),na.rm=T )
  TE[3]<-sum((Est.PD!=0),na.rm=T )
  TE[4]<-sum( (Est.Phi!=0),na.rm=T)/2
  return(TE)  
}



# Eva.count
# Input: 
#       B: the true parameter matrix of Gaussian-Gaussian interactions
#       P: the true parameter matrix of Gaussian-Binary interactions
#       Phi: the true parameter matrix of Binary-Binary interactions
#       Est: the estimated parameters
# Output:
#       TP: number of true positive edges
Eva.count<-function(B, P, Phi, Est, reverse=F){
  p<-dim(B)[1]
  q<-dim(Phi)[1]
  if (reverse==T){
    Est.Phi<-Est[1:p,1:p]
    Est.PC<-Est[(p+1):(p+q),1:p ]
    Est.PD<-t(Est[1:p, (p+1):(p+q)])
    Est.B<-Est[(p+1):(p+q), (p+1):(p+q)]
  }else{
    Est.B<-Est[1:p,1:p]
    Est.PC<-Est[1:p, (p+1):(p+q)]
    Est.PD<-t(Est[(p+1):(p+q),1:p ])
    Est.Phi<-Est[(p+1):(p+q), (p+1):(p+q)]
  }  
  
  #make B and Phi symmetric
  #Est.B<-(Est.B+t(Est.B))/2
  #Est.Phi<-(Est.Phi+t(Est.Phi))/2
  diag(Est.B)<-0
  diag(Est.Phi)<-0
  diag(B)<-0
  diag(Phi)<-0

  TP<-numeric(4)
  TP[1]<-sum( (Est.B!=0) & (B!=0),na.rm=T)/2
  TP[2]<-sum( (Est.PC!=0) & (P!=0),na.rm=T)
  TP[3]<-sum( (Est.PD!=0) & (P!=0),na.rm=T)
  TP[4]<-sum( (Est.Phi!=0) & (Phi!=0),na.rm=T)/2
  return(TP)  
}


edge_eval<-function(est, p,q, B, P, Phi){
  res<-numeric(8)
  
  res[1]<-concord(abs(est$b ),abs(truth$b))[2]
  res[2]<-concord(abs(est$pC ),abs(truth$p))[2]
  res[3]<-concord(abs(est$pD ),abs(truth$p))[2]
  res[4]<-concord(abs(est$phi),abs(truth$phi))[2]
  
  res[5]<-ifelse( is.na(res[1]), NA, sum(est$b!=0))
  res[6]<-ifelse( is.na(res[2]), NA, sum(est$pC!=0))
  res[7]<-ifelse( is.na(res[3]), NA, sum(est$pD!=0))
  res[8]<-ifelse( is.na(res[4]), NA, sum(est$phi!=0))
  
  
  
  return(res)
}


edge_aver<-function(edges, B=NULL, P=NULL, Phi=NULL,count){
  M2<-dim(edges)[1]
  bic_aver<-apply(edges,c(2,3), mean,na.rm=T)*(M2-apply(is.na(edges),c(2,3), sum))/(count-apply(is.na(edges),c(2,3), sum))  
  return(t(bic_aver))
}
