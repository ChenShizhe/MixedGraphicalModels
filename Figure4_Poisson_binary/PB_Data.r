###--------------------------------------------------------------###
## Generating random samples from existing graphs (Poisson and binary)
## Last updated: Mar.27 2014
###--------------------------------------------------------------###
### Note:
##   Data generation can be very slow depending on the choices of 
##   parameters of Gibbs sampler
###--------------------------------------------------------------###



## Arguments of the function
# M1 is the number of graphs to generate
# M2: the number of data sets
# Gibbs.n: the number of iterations between two samples
# burnin: the length of samples to ignore
# low: \alpha_{1s} for the first p/2 Poisson nodes
# high: \alpha_{1s} for the last p/2 Poisson nodes
# p: the number of Gaussian/binary nodes
# size: sample size 

PB_Data<-function(M1,M2,Gibbs.n=20, burnin=200,low=-3,high=0,size,p){
  # Set the alpha_1 for poisson nodes
  meanlow<- low
  meanhigh <- high
  q=p;
  for(iterg in 1:M1){
    
    # Beta: the matrix of edge potentials for P-P edges
    B<-as.matrix(read.csv(file=paste("./Graph/graphGB",2*p, " B", iterg,".csv",sep=""),header=F))
    # Rho: the matrix of edge potentials for P-B edges
    P<-as.matrix(read.csv(file=paste("./Graph/graphGB",2*p, " P", iterg,".csv",sep=""),header=F))
    # Phi: the matrix of edge potentials for B-B edges
    Phi<-as.matrix(read.csv(file=paste("./Graph/graphGB",2*p, " Phi", iterg,".csv",sep=""),header=F))

    alpha1<-c(rep(meanlow,p/2),rep(meanhigh,p/2))
    
    for(iter in 1:M2){
      print(paste("Generate Dataset ",iter, " for Graph ",iterg, sep=""))
      
        rs<-sampler_poisson_binary(size, B=B, P=P, Phi=Phi,a0=alpha1, seedmultiplier=(19+23*iter),
                                 Gibbs.n=Gibbs.n,burnin=burnin)
      write.table(rs,file=paste("./Data/sample",2*p,"PB",iterg, "N",iter, ".txt",sep=""))
    }
  }
  
}
  