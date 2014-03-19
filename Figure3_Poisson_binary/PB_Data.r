###--------------------------------------------------------------###
### Generating random samples from existing graphs (Poisson and binary)
### Last updated: Mar.19th 2014
###--------------------------------------------------------------###
### Note:
###   Data generation can be very slow depending on the choices of 
###   parameters of Gibbs sampler
###--------------------------------------------------------------###

source("../Sources/MGM_Sampler.r")


PB_Data<-function(M1,M2,Gibbs.n=20, burnin=200){
  # Set the alpha_1 for poisson nodes
  meanlow<- -3
  meanhigh <- 0
  
  for(iterg in 1:M1){
    
    B<-as.matrix(read.csv(file=paste("./Graph/graphPB", "40 B",iterg, ".csv",sep=""),header=F))
    P<-as.matrix(read.csv(file=paste("./Graph/graphPB", "40 P",iterg, ".csv",sep=""),header=F))
    Phi<-as.matrix(read.csv(file=paste("./Graph/graphPB", "40 Phi",iterg, ".csv",sep=""),header=F))
    
    p<-dim(P)[1]
    q<-dim(P)[2]
    
    alpha1<-c(rep(meanlow,p/2),rep(meanhigh,p/2))
    
    for(iter in 1:M2){
      set.seed(19+23*iter)  
        rs<-sampler_poisson_binary(400, B=B, P=P, Phi=Phi,a0=alpha1, seedmultiplier=(19+23*iter),
                                 Gibbs.n=20,burnin=200)
      write.table(rs,file=paste("./Data/sample40PB",iterg, "N400", iter, ".txt",sep=""))
    }
  }
  
}
  