###########################################  
## Evaluation functions
## Last updated: Mar.14th 2014
###-------------------------------------###
## This file contains the following functions:
##  1) Eva_total: counting the total number of estimated edges
##  2) Eva_count: Counting the total number of true positives
##  3) BIC_eval: evaluating the performance of three BIC rules for G-B networks
##  4) BIC_aver: taking average of the performance
##  5) edge_aver: taking average of the performance



###-------------------------------------###
## Eva_total ####
## Counting the total number of estimated edges
## Input: 
##       B: the true parameter matrix of Gaussian-Gaussian interactions
##       P: the true parameter matrix of Gaussian-Binary interactions
##       Phi: the true parameter matrix of Binary-Binary interactions
##       Est: the estimated parameters
##       reverse: whether the binary node preceeds the continuous
## Output:
##       TE: total number of estimated edges
###-------------------------------------###

Eva_total<-function(B, P, Phi, Est, reverse=F){
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



###-------------------------------------###
## Eva_count ####
## Counting the total number of true positives
## Input: 
##       B: the true parameter matrix of Gaussian-Gaussian interactions
##       P: the true parameter matrix of Gaussian-Binary interactions
##       Phi: the true parameter matrix of Binary-Binary interactions
##       Est: the estimated parameters
##       reverse: whether the binary node preceeds the continuous
## Output:
##       TP: number of true positive edges
###-------------------------------------###
Eva_count<-function(B, P, Phi, Est, reverse=F){
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
  diag(B)<-0
  diag(Phi)<-0

  TP<-numeric(4)
  TP[1]<-sum( (Est.B!=0) & (B!=0),na.rm=T)/2
  TP[2]<-sum( (Est.PC!=0) & (P!=0),na.rm=T)
  TP[3]<-sum( (Est.PD!=0) & (P!=0),na.rm=T)
  TP[4]<-sum( (Est.Phi!=0) & (Phi!=0),na.rm=T)/2
  return(TP)  
}




###-------------------------------------###
## BIC_eval ####
## Evaluating the performance of three BIC rules for G-B networks
## Input:
##       N: the whole path of estimated neighbours
##       bic: BIC associated with each tuning parameters
##       B: the true parameter matrix of Gaussian-Gaussian interactions
##       P: the true parameter matrix of Gaussian-Binary interactions
##       Phi: the true parameter matrix of Binary-Binary interactions
## Output:
##       Evak: a 3 by 8 matrix that contains the 
##          number of true positives and  total number of each type of edges
##          for the three approaches (each row)
###-------------------------------------###
BIC_eval<-function(N,bic,p,q, B, P, Phi){
  Eval<-matrix(0,3,8)
  
  est<-combine(N[,,BIC_all(bic)],p=p,q=q,cd='node', cc='or',dd='or')  
  Eval[1,]<-c(Eva_count(B=B,P=P,Phi=Phi,Est=est,reverse=F), Eva_total(B=B,P=P,Phi=Phi,Est=est,reverse=F)) 
  
  bytype<-BIC_type(bic,p=p,q=q)
  temp<- rbind(N[1:p,,bytype[1]],N[(p+1):(p+q),,bytype[2]])
  est<-combine(temp,p=p,q=q,cd='node', cc='or',dd='or')
  Eval[2,]<-c(Eva_count(B=B,P=P,Phi=Phi,Est=est,reverse=F), Eva_total(B=B,P=P,Phi=Phi,Est=est,reverse=F)) 
  
  bynode<-BIC_each(bic)
  temp<-N[,,1]
  for(te in 1:(p+q)){
    temp[te,]<- N[te,,bynode[te]]
  }
  est<-combine(temp,p=p,q=q,cd='node', cc='or',dd='or')
  Eval[3,]<-c(Eva_count(B=B,P=P,Phi=Phi,Est=est,reverse=F), Eva_total(B=B,P=P,Phi=Phi,Est=est,reverse=F)) 
  
  return(Eval)
}


###-----------------------------------------###
## BIC_aver ####
## Taking average of the performance
## Input: 
##        edges: output from BIC_eval
##        B: the true Beta matrix (Gaussian-Gaussian)
##        P: the true Rho matrix (Gaussian-Binary)
##        Phi: the true Phi matrix (Binary-Binary)
##        count: the number of datasets tht have been fitted (if it is nondegenerate)
## Output:
##         a 3 by 8 matrix contains
##         Average rank correlation, average number of true edges, and average number of total edges.
##        Or
##         a 3 by 2 matrix contains
##         precision and recall over the whole graph
###-----------------------------------------###
BIC_aver<-function(edges, B, P, Phi,count, type='number'){
  M2<-dim(edges)[1]
  if(type=='number'){
    
    bic_aver<-matrix(0,3,8)
    bic_aver<-apply(edges,c(2,3), mean)*M2/count  
    return(bic_aver)
  } else{
    
    
    diag(B)<-0
    diag(Phi)<-0
    # P= TP+FN, estimated from edges
    Tedges <- sum( B!=0 )/2+ sum( Phi!=0 )/2 + sum(P!=0)  
    
    # Total number of non-edges
    p<-dim(B)[1]; q<-dim(Phi)[1];
    TNedges<- (p+q)*(p+q-1)/2 - Tedges;
    
    # estimate the number of true positives: from edges[,1:4,]
    TP <- apply(edges[,,c(1,2,4)],c(1,2),sum)
    
    # estimate the number of total positives
    TotP<-apply(edges[,,c(5,6,8)],c(1,2),sum)
    
    # number of True negatives = total negatives - false negative
    # False negatives = true edges - true positive
    
    TN <-  TNedges - (Tedges-TP)
    
    #   #calculate the sensitivity
    #   sens<-apply(TP/Tedges,2,mean)*M2/count
    #   spec<-apply(TN/TNedges,2,mean)*M2/count
    
    #calculate the precision and recall
    pres<-apply(TP/TotP,2,mean)*M2/count
    rec<-apply(TP/Tedges,2,mean)*M2/count
    return(t(rbind(pres,rec)))
  }
}


###-----------------------------------------###
## edge_aver ####
## Taking average of the performance
## Input: 
##        edges: a M2 by 8 matrix that contains 
##               the true positives and total number of estimated edges
##        B: the true Beta matrix (Gaussian-Gaussian)
##        P: the true Rho matrix (Gaussian-Binary)
##        Phi: the true Phi matrix (Binary-Binary)
##        count: the number of datasets tht have been fitted (if it is nondegenerate)
## Output:
##         a 3 by 8 matrix contains
##         average number of true edges, and average number of total edges.
###-----------------------------------------###
edge_aver<-function(edges, B=NULL, P=NULL, Phi=NULL,count){
  M2<-dim(edges)[1]
  aver<-apply(edges,c(2,3), mean,na.rm=T)*(M2-apply(is.na(edges),c(2,3), sum))/(count-apply(is.na(edges),c(2,3), sum))  
  return(t(aver))
}


