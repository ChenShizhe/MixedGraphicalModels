###--------------------------------------------------------------###
### Generating graphs for Poisson-Binary network
### Last updated: Mar.13th 2014
###--------------------------------------------------------------###


source("../Sources/MGM_Graph.r")

M1<-5

p<-q<-40

set.seed(12) 

#####

for(i in 1:M1){
  gr<-generator(p=p,q=q, lwb=0.8,upb=1,PB=T)
  B<-gr[[1]]
  P<-gr[[2]]
  Phi<-gr[[3]]    
  write.table(P,file=paste("./Graph/graphPB", "40 P", i,".csv",sep=""),row.names=F,col.names=F,sep=",")  
  write.table(B,file=paste("./Graph/graphPB", "40 B", i,".csv",sep=""),row.names=F,col.names=F,sep=",")  
  write.table(Phi,file=paste("./Graph/graphPB", "40 Phi", i,".csv",sep=""),row.names=F,col.names=F,sep=",")  
}
