###--------------------------------------------------------------###
### Generating random samples for Fig 2
### Last updated: Mar.27 2014
###--------------------------------------------------------------###
### Note:
###   Data generation can be very slow depending on the choices of 
###   parameters of Gibbs sampler
###--------------------------------------------------------------###
source("../Sources/MGM_Sampler.r")


## Arguments of the function
# M2: the number of data sets
# Gibbs.n: the number of iterations between two samples
# burnin: the length of samples to ignore
# p: the number of Gaussian/binary nodes
# size: sample size 

P_Data<-function(M2,Gibbs.n=20,burnin=200,p,size=4000){
  # Beta: the matrix of edge potentials for G-G edges
  B<-as.matrix(read.csv(file=paste("./Graph/graphGB",2*p, " B.csv",sep=""),header=F))
  # Rho: the matrix of edge potentials for G-B edges
  P<-as.matrix(read.csv(file=paste("./Graph/graphGB",2*p, " P.csv",sep=""),header=F))
  # Phi: the matrix of edge potentials for B-B edges
  Phi<-as.matrix(read.csv(file=paste("./Graph/graphGB",2*p, " Phi.csv",sep=""),header=F))
  
    for(iter in 1:M2){
    print(paste("Generate Dataset ",iter, " for the graph with ",2*p," nodes", sep=""))
    rs<-sampler(size, B=B, P=P, Phi=Phi, seedmultiplier=(19+23*iter),
                Gibbs.n=Gibbs.n,burnin=burnin)
    write.table(rs,file=paste("./Data/sample",2*p,"N",  iter, ".txt",sep=""))

  }
}

