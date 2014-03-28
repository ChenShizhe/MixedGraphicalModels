###--------------------------------------------------------------###
### Generating random samples from existing graphs (Poisson and binary)
### Last updated: Mar.27 2014
###--------------------------------------------------------------###
### Note:
###   Data generation can be very slow depending on the choices of 
###   parameters of Gibbs sampler
###--------------------------------------------------------------###

source("../Sources/MGM_Sampler.r")


PB_Data<-function(M1,M2,Gibbs.n=20, burnin=200,low=-3,high=0,size,p){
  # Set the alpha_1 for poisson nodes
  meanlow<- low
  meanhigh <- high
  
  for(iterg in 1:M1){
    
    B<-as.matrix(read.csv(file=paste("./Graph/graphPB", 2*p, " B",iterg, ".csv",sep=""),header=F))
    P<-as.matrix(read.csv(file=paste("./Graph/graphPB", 2*p, " P",iterg, ".csv",sep=""),header=F))
    Phi<-as.matrix(read.csv(file=paste("./Graph/graphPB", 2*p, " Phi",iterg, ".csv",sep=""),header=F))
    
    p<-dim(P)[1]
    q<-dim(P)[2]
    
    alpha1<-c(rep(meanlow,p/2),rep(meanhigh,p/2))
    
    for(iter in 1:M2){
      set.seed(19+23*iter)  
      print(paste("Generate Dataset ",iter, " for Graph ",iterg, sep=""))
      
        rs<-sampler_poisson_binary(size, B=B, P=P, Phi=Phi,a0=alpha1, seedmultiplier=(19+23*iter),
                                 Gibbs.n=Gibbs.n,burnin=burnin)
      write.table(rs,file=paste("./Data/sample",2*p,"PB",iterg, "N",iter, ".txt",sep=""))
    }
  }
  
}
  