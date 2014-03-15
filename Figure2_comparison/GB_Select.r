###--------------------------------------------------------------###
### Neighbourhood selection on the Gaussian-binary networks
### Last updated: Mar.14th 2014
###--------------------------------------------------------------###


source("../Sources/MGM_Evaluation.r")
source("../Sources/MGM_Combine.r")
source("../Sources/MGM_BIC.r")
source("../Sources/MGM_NSelect.r")
source("../Sources/MGM_misc.r")


###################################################################################
####   Counting the true edges and total edges
###################################################################################



total<-100
lambda<-exp(-seq(from= 4, to= -6, length.out=total))
size<-200

M1<-5
M2<-3

for( iterg in 1:M1){
  B<-as.matrix(read.csv(file=paste("./Graph/graphGB",  "20 B",iterg, ".csv",sep=""),header=F))
  P<-as.matrix(read.csv(file=paste("./Graph/graphGB",  "20 P",iterg, ".csv",sep=""),header=F))
  Phi<-as.matrix(read.csv(file=paste("./Graph/graphGB", "20 Phi",iterg, ".csv",sep=""),header=F))

  p<-dim(P)[1]
  q<-dim(P)[2]
  count<-0
  evaluation_ratio<-evaluation_node<-replicate(total,matrix(0,M2,8))
  edges_BIC<-replicate(8,matrix(0,M2,3))
for(iter in 1:M2){
  sample_all<-as.matrix(read.table(file=paste("./Data/sample20GB", iterg, "N400", iter,".txt",sep="" )))
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
write.table(temp, file=paste("./Estimates/BIC/G", iterg, "BICcount200.txt",sep=""))   

temp<-BIC_aver(edges=edges_BIC, B=B, P=P, Phi=Phi,count=count, type='rate')
write.table(temp, file=paste("./Estimates/BIC/G", iterg, "BICrate200.txt",sep=""))   
}