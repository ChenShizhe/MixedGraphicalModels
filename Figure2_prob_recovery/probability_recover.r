###--------------------------------------------------------------###
### Probability of successful recovery on Gaussian-binary network
### Last updated: Mar.27 2014
###--------------------------------------------------------------###

# p: number of Gaussian nodes
# total: number of tuning parameters
# M2: number of data sets
# size: the maximum sample size
# step: increment of sample size

P_recover<-function(p,total=50,M2,size,step=200){
  
  # Beta: the matrix of edge potentials for G-G edges
  B<-as.matrix(read.csv(file=paste("./Graph/graphGB",2*p, " B.csv",sep=""),header=F))
  # Rho: the matrix of edge potentials for G-B edges
  P<-as.matrix(read.csv(file=paste("./Graph/graphGB",2*p, " P.csv",sep=""),header=F))
  # Phi: the matrix of edge potentials for B-B edges
  Phi<-as.matrix(read.csv(file=paste("./Graph/graphGB",2*p, " Phi.csv",sep=""),header=F))
  
  #note: p=q in this experiment
  q<-p;
  
  # Range of tuning parameters to try
  lambdas<-exp(-seq(from= 2, to= -2, length.out=total))
  
  # Examinie the following sample sizes
  samplesize<-seq(from=200, to=size, by=step)
  nsize<-length(samplesize)
  
  
  success.all<-replicate(4,matrix(0,total,M2))
  all.rate<-replicate(4,matrix(0,total,nsize))
  
  for(size in 1:nsize){
    count<-0
    print(paste("Sample size = ", samplesize[size], sep=""))
    for(iter in 1:M2){
      
      sample_limited<-as.matrix(read.table(file=paste("./Data/sample",2*p,"N",iter, ".txt",sep="" )))
      sample_limited<-sample_limited[1:samplesize[size],]
      
      if( is_same(sample_limited[,p+(1:q) ])==T){ success.all[iter,]<-0
      }else{
        count<-count+1
        sds<-as.numeric(apply(sample_limited, 2, sd))
        est_path<-neighbour(dat=sample_limited,p=p,q=q,
                            clambda=lambdas,ratio=1, pf=T, pw=sds)
        for(j in 1: total){
          # assess the performance
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
  write.table(all.rate[,,1], file=paste("./Estimates/",2*p,"/","prob",2*p,"cc.txt",sep=""))
  
  write.table(all.rate[,,2], file=paste("./Estimates/",2*p,"/","prob",2*p,"cd.txt",sep=""))
  
  write.table(all.rate[,,3], file=paste("./Estimates/",2*p,"/","prob",2*p,"dc.txt",sep=""))
  
  write.table(all.rate[,,4], file=paste("./Estimates/",2*p,"/","prob",2*p,"dd.txt",sep=""))
  print("Finished.")
}
