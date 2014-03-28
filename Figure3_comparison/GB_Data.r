###--------------------------------------------------------------###
### Generating random samples from existing graphs (Gaussian and binary)
### Last updated: Mar.19th 2014
###--------------------------------------------------------------###
### Note:
###   Data generation can be very slow depending on the choices of 
###   parameters of Gibbs sampler
###--------------------------------------------------------------###


GB_Data<-function(M1, M2, Gibbs.n=20, burnin=200,size,p){ 
  q<-p
  for(iterg in 1:M1){
    # Beta: the matrix of edge potentials for G-G edges
    B<-as.matrix(read.csv(file=paste("./Graph/graphGB",2*p, " B", iterg,".csv",sep=""),header=F))
    # Rho: the matrix of edge potentials for G-B edges
    P<-as.matrix(read.csv(file=paste("./Graph/graphGB",2*p, " P", iterg,".csv",sep=""),header=F))
    # Phi: the matrix of edge potentials for B-B edges
    Phi<-as.matrix(read.csv(file=paste("./Graph/graphGB",2*p, " Phi", iterg,".csv",sep=""),header=F))
    
    for(iter in 1:M2){
      print(paste("Generate Dataset ",iter, " for Graph ",iterg, sep=""))
      set.seed(19+23*iter)  
      rs<-sampler(size, B=B, P=P, Phi=Phi, seedmultiplier=(19+23*iter),
                  Gibbs.n=Gibbs.n,burnin=burnin)
      write.table(rs,file=paste("./Data/sample",2*p,"GB",iterg,"N", iter, ".txt",sep=""))
    }
  } 
}
