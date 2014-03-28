###--------------------------------------------------------------###
### Neighbourhood selection on the Poisson-binary networks
### Last updated: Mar.27 2014
###--------------------------------------------------------------###
### Note:
###     Need annotation.



PB_Select<-function(M1,M2,size,low=-3,high=0,total=100,maxit=10000, p){
  
  q=p;
  meanlow<- low
  meanhigh <- high
  
  lambda<-exp(-seq(from= 6, to= -4, length.out=total))
  
  #calculate the number we need in order to get 1/2
  gap<-ceiling(log(2)/ log(lambda[2]/lambda[1])) #~7
  
  #####
  ## Also need to pick the correct level of samples
  for( iterg in 1:M1){
    
    
    B<-as.matrix(read.csv(file=paste("./Graph/graphPB", 2*p, " B",iterg, ".csv",sep=""),header=F))
    P<-as.matrix(read.csv(file=paste("./Graph/graphPB", 2*p, " P",iterg, ".csv",sep=""),header=F))
    Phi<-as.matrix(read.csv(file=paste("./Graph/graphPB", 2*p, " Phi",iterg, ".csv",sep=""),header=F))
    
    
    count<-0
    edges_BIC<-array(0,c(M2,8,3))
    evaluation_bic2<-evaluation_or<-evaluation_and<-evaluation<-replicate(total,matrix(0,M2,8))
    for(iter in 1:M2){
      print(paste("Fitting the MGM on Dataset ",iter, " for Graph ",iterg, sep=""))
      
      sample_all<-as.matrix(read.table(file=paste("./Data/sample", 2*p, "PB", iterg, "N", iter,".txt",sep="" )))
      sample_limited<-sample_all[1:size, ]
      if(is_same(sample_limited)==TRUE){
      } else {
        count<-count+1
        
        sds<-as.numeric(apply(sample_limited, 2, sd))
        est_path<-neighbour_PB(dat=sample_limited,p=p,q=q, 
                               clambda=lambda, ratio=1, pf=T, pw=sds,maxit=maxit) 
        #preparation for applying the seleciton rules
        
        parameter_bic2<-parameter<-list(coef=1, a0=1)
        
        parameter$coef<- rbind(cbind(B,P), cbind(t(P),Phi))
        parameter$a0 <- c(rep(meanlow,p/2),rep(meanhigh,p/2))
        
        bytype<-BIC_type(bic=est_path$bic,p=p,q=q)      
        
        parameter_bic2$a0 <- est_path$Intercept[1:p,bytype[1]+gap] 
        parameter_bic2$coef<- rbind(est_path$N[1:p,,bytype[1]+gap],est_path$N[p+(1:q),,bytype[2]+7])
        
        for(j in 1:total){ 
          
          ggraph<-abs(est_path$N[,,j]) 
          
          est<-combine_select(N=ggraph,p=p,q=q,parameter=parameter)
          evaluation[iter,,j]<-c(Eva_count(B=B,P=P,Phi=Phi,Est=est,reverse=F), Eva_total(B=B,P=P,Phi=Phi,Est=est,reverse=F))
          
          est<-combine_select(N=ggraph,p=p,q=q,parameter=parameter_bic2)
          evaluation_bic2[iter,,j]<-c(Eva_count(B=B,P=P,Phi=Phi,Est=est,reverse=F), Eva_total(B=B,P=P,Phi=Phi,Est=est,reverse=F))
          
          est<-combine(ggraph,p=p,q=q, cd='or', cc='or',dd='or')
          evaluation_or[iter,,j]<-c(Eva_count(B=B,P=P,Phi=Phi,Est=est,reverse=F), Eva_total(B=B,P=P,Phi=Phi,Est=est,reverse=F))
          
          est<-combine(ggraph,p=p,q=q, cd='and', cc='or',dd='or')
          evaluation_and[iter,,j]<-c(Eva_count(B=B,P=P,Phi=Phi,Est=est,reverse=F), Eva_total(B=B,P=P,Phi=Phi,Est=est,reverse=F))
          
        }
        
      }
    }
    result<-edge_aver(edges=evaluation, count=count)
    
    write.table(result, file=paste("./Estimates/G", iterg, "PB_S.txt",sep=""))
    
    result<-edge_aver(edges=evaluation_and, count=count)
    
    write.table(result, file=paste("./Estimates/G", iterg, "PB_AND.txt",sep=""))
    
    result<-edge_aver(edges=evaluation_or, count=count)
    
    write.table(result, file=paste("./Estimates/G", iterg, "PB_OR.txt",sep=""))
    
    result<-edge_aver(edges=evaluation_bic2, count=count)
    
    write.table(result, file=paste("./Estimates/G", iterg, "PB_BIC2.txt",sep=""))
    
  }
  
}