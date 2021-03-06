###--------------------------------------------------------------###
## This file contains the following functions:
## 1) Neighbourhood selection on the Gaussian-binary networks
## 2) Neighbourhood selection on the Ising models
## 3) Neighbourhood selection on the Gaussian graphical models
## 4) The graphical lasso
## Last updated: Mar.27 2014
###--------------------------------------------------------------###


## Arguments of the functions
# M1: the number of graphs
# M2: the number of datasets for each graph
# size: the sample size
# total: the number of tuning parameters 
# p: the number of Gaussian nodes, which equals to the number of binary nodes.



GB_Ising<-function(M1,M2,size,total=100,p){
  q<-p
  lambda<-exp(-seq(from= 4, to= -6, length.out=total))
  
  for( iterg in 1:M1){
    
    # Beta: the matrix of edge potentials for G-G edges
    B<-as.matrix(read.csv(file=paste("./Graph/graphGB",2*p, " B", iterg,".csv",sep=""),header=F))
    # Rho: the matrix of edge potentials for G-B edges
    P<-as.matrix(read.csv(file=paste("./Graph/graphGB",2*p, " P", iterg,".csv",sep=""),header=F))
    # Phi: the matrix of edge potentials for B-B edges
    Phi<-as.matrix(read.csv(file=paste("./Graph/graphGB",2*p, " Phi", iterg,".csv",sep=""),header=F))
    
    p<-dim(P)[1]
    q<-dim(P)[2]
    count<-0
    evaluation_or<-replicate(total,matrix(0,M2,8))
    for(iter in 1:M2){
      print(paste("Fitting the Ising model on Dataset ",iter, " for Graph ",iterg, sep=""))
      
      sample_all<-as.matrix(read.table(file=paste("./Data/sample",2*p, "GB", iterg, "N", iter,".txt",sep="" )))
      sample_limited<-sample_all[1:size,]
      
      if(is_same(sample_limited)==TRUE){
        
      } else {
        count<-count+1
        sample_limited<-sample_limited-rep(1,size)%*%t(apply(sample_limited,2,mean))
        sample_limited<-sign(sample_limited)
        sds<-as.numeric(apply(sample_limited, 2, sd))
        est_path<-neighbour_Ising(dat=sample_limited,p=p,q=q, clambda=lambda, pf=T, pw=sds)       
        for(j in 1: total){
          ggraph<-est_path[,,j]
          
          est<-combine(ggraph,p=p,q=q, cd='or', cc='or',dd='or')
          evaluation_or[iter,,j]<-c(Eva_count(B=B,P=P,Phi=Phi,Est=est,reverse=F), Eva_total(B=B,P=P,Phi=Phi,Est=est,reverse=F))         
        }
      }
    }
    result<-edge_aver(edges=evaluation_or, count=count)
    write.table(result, file=paste("./Estimates/Ising/G", iterg, "I_or.txt",sep=""))
  }
  
}


GB_Gaussian<-function(M1,M2,size,total=100,p){
  q<-p
  lambda<-exp(-seq(from= 4, to= -6, length.out=total))
  for( iterg in 1:M1){
    # Beta: the matrix of edge potentials for G-G edges
    B<-as.matrix(read.csv(file=paste("./Graph/graphGB",2*p, " B", iterg,".csv",sep=""),header=F))
    # Rho: the matrix of edge potentials for G-B edges
    P<-as.matrix(read.csv(file=paste("./Graph/graphGB",2*p, " P", iterg,".csv",sep=""),header=F))
    # Phi: the matrix of edge potentials for B-B edges
    Phi<-as.matrix(read.csv(file=paste("./Graph/graphGB",2*p, " Phi", iterg,".csv",sep=""),header=F))
    
    p<-dim(P)[1]
    q<-dim(P)[2]
    
    count<-0
    evaluation_or<-replicate(total,matrix(0,M2,8))
    for(iter in 1:M2){
      
      print(paste("Fitting the Gaussian graphical model on Dataset ",iter, " for Graph ",iterg, sep=""))
      
      sample_all<-as.matrix(read.table(file=paste("./Data/sample",2*p, "GB", iterg, "N", iter,".txt",sep="" )))
      sample_limited<-sample_all[1:size,]
      if(is_same(sample_limited)==TRUE){
        
      } else {
        count<-count+1
        sds<-as.numeric(apply(sample_limited, 2, sd))
        est_path<-neighbour_Gaussian(dat=sample_limited,p=p,q=q, 
                                     clambda=lambda, pf=T, pw=sds) 
        
        for(j in 1: total){
          ggraph<-est_path[,,j]
          
          est<-combine(ggraph,p=p,q=q, cd='or', cc='or',dd='or')
          evaluation_or[iter,,j]<-c(Eva_count(B=B,P=P,Phi=Phi,Est=est,reverse=F), Eva_total(B=B,P=P,Phi=Phi,Est=est,reverse=F))       
          
        }
      }
    }
    result<-edge_aver(edges=evaluation_or, count=count)
    write.table(result, file=paste("./Estimates/MB/G", iterg, "MB_or.txt",sep=""))
  }
}

GB_Glasso<-function(M1,M2,size,total=100,p){
  q<-p
  lambda<-exp(-seq(from= 4, to= -6, length.out=total))
  
  for( iterg in 1:M1){
    # Beta: the matrix of edge potentials for G-G edges
    B<-as.matrix(read.csv(file=paste("./Graph/graphGB",2*p, " B", iterg,".csv",sep=""),header=F))
    # Rho: the matrix of edge potentials for G-B edges
    P<-as.matrix(read.csv(file=paste("./Graph/graphGB",2*p, " P", iterg,".csv",sep=""),header=F))
    # Phi: the matrix of edge potentials for B-B edges
    Phi<-as.matrix(read.csv(file=paste("./Graph/graphGB",2*p, " Phi", iterg,".csv",sep=""),header=F))
    
    p<-dim(P)[1]
    q<-dim(P)[2]
  
    rho.vec<-sqrt(2*log(p)/size)*lambda/2
    
    count<-0
    evaluation<-replicate(total,matrix(0,M2,8))
    for(iter in 1:M2){
      
      print(paste("Fitting the Glasso on Dataset ",iter, " for Graph ",iterg, sep=""))
      
      sample_all<-as.matrix(read.table(file=paste("./Data/sample",2*p, "GB", iterg, "N", iter,".txt",sep="" )))
      sample_limited<-sample_all[1:size,]
      if(is_same(sample_limited)==TRUE){
      } else {
        count<-count+1
        est_path<-glassopath(s=cor(sample_limited)*(size-1)/size, rholist=rho.vec, approx=F, penalize.diagonal=F,thr=1e-07)
        for(j in 1: total){
          ggraph<-abs(est_path$wi[,,j]) 
          est<-combine(ggraph,p=p,q=q,cd='node', cc='or',dd='or')  
          evaluation[iter,,j]<-c(Eva_count(B=B,P=P,Phi=Phi,Est=est,reverse=F), Eva_total(B=B,P=P,Phi=Phi,Est=est,reverse=F))     
        }
      }
    }
    result<-edge_aver(edges=evaluation, count=count)
    write.table(result, file=paste("./Estimates/Glasso/G", iterg, "GR.txt",sep="")) 
  }
  
}


GB_Select<-function(M1,M2,size,total=100,p){  
  q<-p
  lambda<-exp(-seq(from= 4, to= -6, length.out=total))
  
  for( iterg in 1:M1){
    # Beta: the matrix of edge potentials for G-G edges
    B<-as.matrix(read.csv(file=paste("./Graph/graphGB",2*p, " B", iterg,".csv",sep=""),header=F))
    # Rho: the matrix of edge potentials for G-B edges
    P<-as.matrix(read.csv(file=paste("./Graph/graphGB",2*p, " P", iterg,".csv",sep=""),header=F))
    # Phi: the matrix of edge potentials for B-B edges
    Phi<-as.matrix(read.csv(file=paste("./Graph/graphGB",2*p, " Phi", iterg,".csv",sep=""),header=F))
    
    count<-0
    evaluation_ratio<-evaluation_node<-replicate(total,matrix(0,M2,8))
    edges_BIC<-replicate(8,matrix(0,M2,3))
    for(iter in 1:M2){
      
      print(paste("Fitting the MGM on Dataset ",iter, " for Graph ",iterg, sep=""))
      
      sample_all<-as.matrix(read.table(file=paste("./Data/sample",2*p, "GB", iterg, "N", iter,".txt",sep="" )))
      sample_limited<-sample_all[1:size,]
      if(is_same(sample_limited)==TRUE){
        
      } else {
        count<-count+1
        sds<-as.numeric(apply(sample_limited, 2, sd))
        est_path<-neighbour(dat=sample_limited,p=p,q=q, 
                            clambda=lambda, ratio=1, pf=T, pw=sds) 
        
        edges_BIC[iter,,]<-BIC_eval(N=est_path$N,bic=est_path$bic,p=p,q=q, B=B, P=P, Phi=Phi)
        
        
        lambdas<-BIC_type(est_path$bic, p = p, q=q)
        ratio_bic<-lambda[lambdas[1]]/lambda[lambdas[2]]
        est_ratio<-neighbour(dat=sample_limited,p=p,q=q, 
                             clambda=lambda, ratio=ratio_bic, pf=T, pw=sds) 
        for(j in 1: total){
          ggraph<-est_ratio$N[,,j]
          est<-combine(ggraph,p=p,q=q,cd='node', cc='or',dd='or')  
          evaluation_ratio[iter,,j]<-c(Eva_count(B=B,P=P,Phi=Phi,Est=est,reverse=F), Eva_total(B=B,P=P,Phi=Phi,Est=est,reverse=F)) 
        }
      }
    }
    result<-edge_aver(edges=evaluation_ratio, count=count)
    write.table(result, file=paste("./Estimates/MGM/G", iterg, "MGM_ratio.txt",sep=""))
    
    temp<-BIC_aver(edges=edges_BIC, B=B, P=P, Phi=Phi,count=count, type='number')
    write.table(temp, file=paste("./Estimates/BIC/G", iterg, "BICcount.txt",sep=""))   
    
    temp<-BIC_aver(edges=edges_BIC, B=B, P=P, Phi=Phi,count=count, type='rate')
    write.table(temp, file=paste("./Estimates/BIC/G", iterg, "BICrate.txt",sep=""))   
  }  
}
