###-------------------------------------###
## Tuning parameter selection using BIC
## Last updated: Mar.14th 2014
###-------------------------------------###
## This file contains the following functions:
## 1) N2loglklh: a function that calculate the -2loglikelihood
## 2) BIC: it returns the BIC for each node-wise regression
## 3) BIC_all: it returns the tuning parameter that minize the overall BIC 
## 4) BIC_type: it returns one tuning parameters that minize BIC for each time
## 5) BIC_each: it returns one turning parameter for each node-wise regression



###-----------------------------------------###
## N2Loglklh ####
## Input:
##       beta: estimated coefficients
##       a0: estimated intercept
##       y: outcomes (n times 1)
##       x: covariate matrix (n times p-1 )
##       dist: an indicator of the distribution (B: binary, P: Poisson, G: Gaussian)
## Output: The -2 loglikelihood 
###-----------------------------------------###
N2Loglklh<-function(beta, a0, y, x,dist=T){
  if( dist=='B') {
    y <- (y > mean(y))
    Est_Prob<-expit(x%*%beta+a0)
    res<- -2*sum(log(Est_Prob^y)+log( (1-Est_Prob)^(1-y)))
  } else if (dist == "P") {
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
## BIC ####
## calculating the BIC using loglikelihood and estimated coefficient
## Input:
##       beta: estimated coefficients
##       nlg: 2 times negative loglikelihood (output from N2Loglklh)
##       n: number of observations
## Output: The BIC for this node
###-----------------------------------------###
BIC<-function(nlg,n,beta){
  res<-sum(beta!=0)*log(n)+nlg
  return(res)
}


###-----------------------------------------###
## BIC_all ####
## Finding the tuning parameter that minimizes the sum of BIC across all nodes 
## Input:
##       bic: a matrix of the BIC for every node (rows) and every tuning parameter (columns)
## Output: The index of the tuning parameter that minimizes the overall BICs
###-----------------------------------------###
BIC_all<-function(bic){
  res<-which.min(apply(bic,2,sum))
  return(res)
}

###-----------------------------------------###
## BIC_type ####
##  Finding two tuning parameter which are minimzers of the two types of nodes, repsectively  
## Input:
##       bic: a matrix of the BIC for every node (rows) and every tuning parameter (columns)
##       p: the number of Gaussian nodes
##       q: the number of binary nodes
## Output: The indices of the two tuning parameter for Gaussian nodes and binary nodes
###-----------------------------------------###
BIC_type<-function(bic,p,q){
  type1<-BIC_all(bic[1:p,])
  type2<-BIC_all(bic[(p+1):(p+q),])
  res<-c(type1, type2)
  return(res)
}

###-----------------------------------------###
## BIC_each  ####
##  Finding one tuning parameter for each node-wise regression 
## Input:
##       bic: a matrix of the BIC for every node (rows) and every tuning parameter (columns)
## Output: Indices for all the nodes
###-----------------------------------------###
BIC_each<-function(bic){
  res<-apply(bic,1,which.min)
  return(res)
}










 