###--------------------------------------------------------------###
### Generating random samples for Fig 1
### Last updated: Mar.13th 2014
###--------------------------------------------------------------###
### Note:
###   Data generation can be very slow depending on the choices of 
###   parameters of Gibbs sampler
###--------------------------------------------------------------###
source("../Sources/MGM_Sampler.r")

M2<-10 # Number of datasets for each graph


B<-as.matrix(read.csv(file=paste("./Graph/graphGB", "30Fix B.csv",sep=""),header=F))
P<-as.matrix(read.csv(file=paste("./Graph/graphGB", "30Fix P.csv",sep=""),header=F))
Phi<-as.matrix(read.csv(file=paste("./Graph/graphGB", "30Fix Phi.csv",sep=""),header=F))


for(iter in 1:M2){
  set.seed(19+23*iter)  
  rs<-sampler(10000, B=B, P=P, Phi=Phi, seedmultiplier=(19+23*iter),
                  Gibbs.n=500,burnin=3000)
  write.table(rs,file=paste("sample30", "N10000", iter, ".txt",sep=""))
}



B<-as.matrix(read.csv(file=paste("./Graph/graphGB", "60Fix B.csv",sep=""),header=F))
P<-as.matrix(read.csv(file=paste("./Graph/graphGB", "60Fix P.csv",sep=""),header=F))
Phi<-as.matrix(read.csv(file=paste("./Graph/graphGB", "60Fix Phi.csv",sep=""),header=F))


for(iter in 1:M2){
  set.seed(19+23*iter)  
  rs<-sampler(10000, B=B, P=P, Phi=Phi, seedmultiplier=(19+23*iter),
              Gibbs.n=500,burnin=3000)
  write.table(rs,file=paste("sample60", "N10000", iter, ".txt",sep=""))
}



B<-as.matrix(read.csv(file=paste("./Graph/graphGB", "120Fix B.csv",sep=""),header=F))
P<-as.matrix(read.csv(file=paste("./Graph/graphGB", "120Fix P.csv",sep=""),header=F))
Phi<-as.matrix(read.csv(file=paste("./Graph/graphGB", "120Fix Phi.csv",sep=""),header=F))

for(iter in 1:M2){
  set.seed(19+23*iter)  
  rs<-sampler(10000, B=B, P=P, Phi=Phi, seedmultiplier=(19+23*iter),
                  Gibbs.n=500,burnin=3000)
  write.table(rs,file=paste("sample120", "N10000", iter, ".txt",sep=""))
}

