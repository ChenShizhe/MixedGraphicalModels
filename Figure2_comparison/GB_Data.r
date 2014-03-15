###--------------------------------------------------------------###
### Generating random samples from existing graphs (Gaussian and binary)
### Last updated: Mar.14th 2014
###--------------------------------------------------------------###
### Note:
###   Data generation can be very slow depending on the choices of 
###   parameters of Gibbs sampler
###--------------------------------------------------------------###

source("../Sources/MGM_Sampler.r")


M1<-5 # Number of graphs
M2<-3 # Number of datasets for each graph


for(iterg in 1:M1){
  
B<-as.matrix(read.csv(file=paste("./Graph/graphGB", "20 B",iterg, ".csv",sep=""),header=F))
P<-as.matrix(read.csv(file=paste("./Graph/graphGB", "20 P",iterg, ".csv",sep=""),header=F))
Phi<-as.matrix(read.csv(file=paste("./Graph/graphGB", "20 Phi",iterg, ".csv",sep=""),header=F))

p<-dim(P)[1]
q<-dim(P)[2]

for(iter in 1:M2){
  set.seed(19+23*iter)  
  ###For simulation in the paper: Gibbs.n=500, burnin=3000
  rs<-sampler(400, B=B, P=P, Phi=Phi, seedmultiplier=(19+23*iter),
              Gibbs.n=20,burnin=200)
  write.table(rs,file=paste("./Data/sample20GB",iterg, "N400", iter, ".txt",sep=""))
}
}
