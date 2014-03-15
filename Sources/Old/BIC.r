###-------------------------------------###
## Evaluation functions for BIC
## Created: Feb.16th 2014
## Updated: Feb.25nd 2014
###-------------------------------------###
### This file contains the following functions:
### 1) N2loglklh: a function that calculate the -2loglikelihood for Gaussian and Binary nodes
### 2) BIC: a function that returns the BIC for each regression
### 3) BIC_all
### 4) BIC_type
### 5) BIC_each
### 6) BIC_eval


###-----------------------------------------###
## N2Loglklh
## Input:
##       beta: estimated coefficients
##       a0: estimated intercept
##       y: outcomes (n times 1)
##       x: covariate matrix (n times p-1 )
##       logistic: an indicator of the distribution (T: binary, F: Gaussian)
## Output: The -2 loglikelihood 
## updated this to include the Poisson case

N2Loglklh<-function(beta, a0, y, x,logistic=T){
  if( logistic==T) {
    y <- (y > mean(y))
    Est_Prob<-expit(x%*%beta+a0)
    res<- -2*sum(log(Est_Prob^y)+log( (1-Est_Prob)^(1-y)))

  } else if (logistic == "P") {
    eta<-x%*%beta+a0
    res<- -2*sum(y*eta - log(factorial(y))- exp(eta) )
  } else  {
    res<- length(y)*log(sum((y-x%*%beta-a0)^2))
  }
  
 
  return(res)
}


expit<-function(z){
   1/(1+exp(-z))
 }


###-----------------------------------------###
## BIC: calculating the BIC using loglikelihood and estimated coefficient
## Input:
##       beta: estimated coefficients
##       nlg: 2 times negative loglikelihood (output from N2Loglklh)
##       n: number of observations
## Output: The BIC for this node
BIC<-function(nlg,n,beta){
  res<-sum(beta!=0)*log(n)+nlg
  return(res)
}


###-----------------------------------------###
## BIC_all: Finding the tuning parameter that minimizes the sum of BIC across all nodes 
## Input:
##       bic: a matrix of the BIC for every node (rows) and every tuning parameter (columns)
## Output: The index of the tuning parameter that minimizes the overall BICs
BIC_all<-function(bic){
  res<-which.min(apply(bic,2,sum))
  return(res)
}

###-----------------------------------------###
## BIC_type: Finding two tuning parameter that minimize the sum of BIC across nodes of eithe types  
## Input:
##       bic: a matrix of the BIC for every node (rows) and every tuning parameter (columns)
##       p: the number of Gaussian nodes
##       q: the number of binary nodes
## Output: The indices of the two tuning parameter for Gaussian nodes and binary nodes
BIC_type<-function(bic,p,q){
  type1<-BIC_all(bic[1:p,])
  type2<-BIC_all(bic[(p+1):(p+q),])
  res<-c(type1, type2)
  return(res)
}

###-----------------------------------------###
## BIC_each: Finding one tuning parameter for each node 
## Input:
##       bic: a matrix of the BIC for every node (rows) and every tuning parameter (columns)
## Output: Indices for all the nodes
BIC_each<-function(bic){
  res<-apply(bic,1,which.min)
  return(res)
}

###-----------------------------------------###
## BIC_eval: a function that evaluates the performance of three BIC rules.
## Note: this function needs to be cleaned
BIC_eval<-function(N,bic,p,q, B, P, Phi){
  Eval<-matrix(0,3,12)
    
  est<-combine_gen(N[,,BIC_all(bic)],p=p,q=q,cd='node', cc='or',dd='or')  
  Eval[1,]<-c(Eva.count(B=B,P=P,Phi=Phi,Est=est,reverse=F), Eva.total(B=B,P=P,Phi=Phi,Est=est,reverse=F),
              Eva.rank(B=B,P=P,Phi=Phi,Est=est,reverse=F)) 
  
  bytype<-BIC_type(bic,p=p,q=q)
  temp<- rbind(N[1:p,,bytype[1]],N[(p+1):(p+q),,bytype[2]])
  est<-combine_gen(temp,p=p,q=q,cd='node', cc='or',dd='or')
  Eval[2,]<-c(Eva.count(B=B,P=P,Phi=Phi,Est=est,reverse=F), Eva.total(B=B,P=P,Phi=Phi,Est=est,reverse=F),
              Eva.rank(B=B,P=P,Phi=Phi,Est=est,reverse=F)) 
  
  bynode<-BIC_each(bic)
  temp<-N[,,1]
  for(te in 1:(p+q)){
    temp[te,]<- N[te,,bynode[te]]
  }
  est<-combine_gen(temp,p=p,q=q,cd='node', cc='or',dd='or')
  Eval[3,]<-c(Eva.count(B=B,P=P,Phi=Phi,Est=est,reverse=F), Eva.total(B=B,P=P,Phi=Phi,Est=est,reverse=F),
              Eva.rank(B=B,P=P,Phi=Phi,Est=est,reverse=F)) 

  return(Eval)
}

BIC_eval_PB<-function(N,bic,p,q, B, P, Phi,parameter){
  Eval<-matrix(0,3,8)
  
  est<-combine_smart(estimate=N[,,BIC_all(bic)],p=p,q=q,parameter=parameter)
  Eval[1,]<-c(Eva.count(B=B,P=P,Phi=Phi,Est=est,reverse=F), Eva.total(B=B,P=P,Phi=Phi,Est=est,reverse=F)) 
  
  bytype<-BIC_type(bic,p=p,q=q)
  temp<- rbind(N[1:p,,bytype[1]],N[(p+1):(p+q),,bytype[2]])
  est<-combine_smart(estimate=temp,p=p,q=q,parameter=parameter)
  Eval[2,]<-c(Eva.count(B=B,P=P,Phi=Phi,Est=est,reverse=F), Eva.total(B=B,P=P,Phi=Phi,Est=est,reverse=F)) 
  
  bynode<-BIC_each(bic)
  temp<-N[,,1]
  for(te in 1:(p+q)){
    temp[te,]<- N[te,,bynode[te]]
  }
  est<-combine_smart(estimate=temp,p=p,q=q,parameter=parameter)
  Eval[3,]<-c(Eva.count(B=B,P=P,Phi=Phi,Est=est,reverse=F), Eva.total(B=B,P=P,Phi=Phi,Est=est,reverse=F)) 
  
  return(Eval)
}


###-----------------------------------------###
## BIC_aver
## Input: 
##        edges: output from BIC_eval
##        B: the true Beta matrix (Gaussian-Gaussian)
##        P: the true Rho matrix (Gaussian-Binary)
##        Phi: the true Phi matrix (Binary-Binary)
##        count: the number of datasets tht have been fitted (if it is nondegenerate)
## Output:
#         a 12 by 3 matrix contains
#         Average rank correlation, average number of true edges, and average number of total edges.
#        Or
#         a 3 by 2 matrix contains
#         sensitivity and specificity of the whole graph
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











 