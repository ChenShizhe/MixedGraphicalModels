###--------------------------------------------------------------###
### Probability of successful recovery on Gaussian-binary network
### Graph with 120 nodes
### Last updated: Mar.21st 2014
###--------------------------------------------------------------###



B<-as.matrix(read.csv(file=paste("./Graph/graphGB", "60 B.csv",sep=""),header=F))
P<-as.matrix(read.csv(file=paste("./Graph/graphGB", "60 P.csv",sep=""),header=F))
Phi<-as.matrix(read.csv(file=paste("./Graph/graphGB", "60 Phi.csv",sep=""),header=F))

###
p<-dim(P)[1]
q<-dim(P)[2]

total<-50
lambdas<-exp(-seq(from= 2, to= -2, length.out=total))

M2<-5
samplesize<-seq(from=200, to=6000, by=200)
nsize<-length(samplesize)

### Number of replicates
success.all<-replicate(4,matrix(0,total,M2))
all.rate<-replicate(4,matrix(0,total,nsize))

for(size in 1:nsize){
  count<-0
  for(iter in 1:M2){
    
    sample_limited<-as.matrix(read.table(file=paste("./Data/sample60N", iter, ".txt",sep="" )))
    sample_limited<-sample_limited[1:samplesize[size],]
    
    if( is_same(sample_limited[,p+(1:q) ])==T){ success.all[iter,]<-0
    }else{
      count<-count+1
      sds<-as.numeric(apply(sample_limited, 2, sd))
      est_path<-neighbour(dat=sample_limited,p=p,q=q,
                                clambda=lambdas,ratio=1, pf=T, pw=sds)
      for(j in 1: total){
        est<-combine_type(est_path$N[,,j],p=p,q=q, cc="and",dd="and", neg=TRUE)
        diag(est[[1]])<-1
        diag(est[[4]])<-1
        cc<-sum(sign(est[[1]])!=sign(B))==0 
        cd<-sum(sign(est[[2]])!=sign(P))==0
        dc<-sum(sign(est[[3]])!=sign(P))==0  
        dd<-sum(sign(est[[4]])!=sign(Phi))==0  
        success.all[j, iter,1]<-cc
        success.all[j, iter,2]<-cd
        success.all[j, iter,3]<-dc
        success.all[j, iter,4]<-dd
      }
    }  
  }
  
  for(k in 1:4){
    for(j in 1:total){
      all.rate[j, size,k]<-mean(success.all[j,,k])*M2/count
    } 
  }
}
write.table(all.rate[,,1], file=paste("./Estimates/120/","prob120cc.txt",sep=""))

write.table(all.rate[,,2], file=paste("./Estimates/120/","prob120cd.txt",sep=""))

write.table(all.rate[,,3], file=paste("./Estimates/120/","prob120dc.txt",sep=""))

write.table(all.rate[,,4], file=paste("./Estimates/120/","prob120dd.txt",sep=""))

