###--------------------------------------------------------------###
### Generating random samples from existing graphs (Gaussian and binary)
### Last updated: Mar.19th 2014
###--------------------------------------------------------------###
### Note:
###   Data generation can be very slow depending on the choices of 
###   parameters of Gibbs sampler
###--------------------------------------------------------------###


GB_Data<-function(M1, M2, Gibbs.n=20, burnin=200,size){ 
  for(iterg in 1:M1){
    
    B<-as.matrix(read.csv(file=paste("./Graph/graphGB", "20 B",iterg, ".csv",sep=""),header=F))
    P<-as.matrix(read.csv(file=paste("./Graph/graphGB", "20 P",iterg, ".csv",sep=""),header=F))
    Phi<-as.matrix(read.csv(file=paste("./Graph/graphGB", "20 Phi",iterg, ".csv",sep=""),header=F))
    
    p<-dim(P)[1]
    q<-dim(P)[2]
    
    for(iter in 1:M2){
      set.seed(19+23*iter)  
      rs<-sampler(size, B=B, P=P, Phi=Phi, seedmultiplier=(19+23*iter),
                  Gibbs.n=Gibbs.n,burnin=burnin)
      write.table(rs,file=paste("./Data/sample20GB",iterg,"N", iter, ".txt",sep=""))
    }
  } 
}
