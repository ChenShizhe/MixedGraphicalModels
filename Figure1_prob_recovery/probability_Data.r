###--------------------------------------------------------------###
### Generating random samples for Fig 1
### Last updated: Mar.13th 2014
###--------------------------------------------------------------###
### Note:
###   Data generation can be very slow depending on the choices of 
###   parameters of Gibbs sampler
###--------------------------------------------------------------###
source("../Sources/MGM_Sampler.r")


P_Data<-function(M2,Gibbs.n=20,burnin=200){
  B<-as.matrix(read.csv(file=paste("./Graph/graphGB", "30 B.csv",sep=""),header=F))
  P<-as.matrix(read.csv(file=paste("./Graph/graphGB", "30 P.csv",sep=""),header=F))
  Phi<-as.matrix(read.csv(file=paste("./Graph/graphGB", "30 Phi.csv",sep=""),header=F))
  
  
  for(iter in 1:M2){
    set.seed(19+23*iter)  
    rs<-sampler(10000, B=B, P=P, Phi=Phi, seedmultiplier=(19+23*iter),
                Gibbs.n=Gibbs.n,burnin=burnin)
    write.table(rs,file=paste("sample30N",  iter, ".txt",sep=""))
  }
  
  
  
  B<-as.matrix(read.csv(file=paste("./Graph/graphGB", "60 B.csv",sep=""),header=F))
  P<-as.matrix(read.csv(file=paste("./Graph/graphGB", "60 P.csv",sep=""),header=F))
  Phi<-as.matrix(read.csv(file=paste("./Graph/graphGB", "60 Phi.csv",sep=""),header=F))
  
  
  for(iter in 1:M2){
    set.seed(19+23*iter)  
    rs<-sampler(10000, B=B, P=P, Phi=Phi, seedmultiplier=(19+23*iter),
                Gibbs.n=Gibbs.n,burnin=burnin)
    write.table(rs,file=paste("sample60N", iter, ".txt",sep=""))
  }
  
  
  
  B<-as.matrix(read.csv(file=paste("./Graph/graphGB", "120 B.csv",sep=""),header=F))
  P<-as.matrix(read.csv(file=paste("./Graph/graphGB", "120 P.csv",sep=""),header=F))
  Phi<-as.matrix(read.csv(file=paste("./Graph/graphGB", "120 Phi.csv",sep=""),header=F))
  
  for(iter in 1:M2){
    set.seed(19+23*iter)  
    rs<-sampler(10000, B=B, P=P, Phi=Phi, seedmultiplier=(19+23*iter),
                Gibbs.n=Gibbs.n,burnin=burnin)
    write.table(rs,file=paste("sample120N", iter, ".txt",sep=""))
  }
  
}
