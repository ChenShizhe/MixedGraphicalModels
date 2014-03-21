###--------------------------------------------------------------###
### Generating graphs for Fig 1
### Last updated: Mar.13th 2014
###--------------------------------------------------------------###

source("../Sources/MGM_Graph.r")


P_Graph<-function(lwb,upb){
  q<-p<-30
  set.seed(12) 
  gr<-generator(p=p, q=q,lwb=lwb, upb=upb)
  
  
  write.table(gr[[1]],file=paste("./Graph/graphGB", "30 B.csv",sep=""),row.names=F,col.names=F,sep=",")
  write.table(gr[[2]],file=paste("./Graph/graphGB", "30 P.csv",sep=""),row.names=F,col.names=F,sep=",")
  write.table(gr[[3]],file=paste("./Graph/graphGB", "30 Phi.csv",sep=""),row.names=F,col.names=F,sep=",")
  
  
  q<-p<-60
  set.seed(12) 
  gr<-generator(p=p, q=q, lwb=lwb, upb=upb)
  
  
  write.table(gr[[1]],file=paste("./Graph/graphGB", "60 B.csv",sep=""),row.names=F,col.names=F,sep=",")
  write.table(gr[[2]],file=paste("./Graph/graphGB", "60 P.csv",sep=""),row.names=F,col.names=F,sep=",")
  write.table(gr[[3]],file=paste("./Graph/graphGB", "60 Phi.csv",sep=""),row.names=F,col.names=F,sep=",")
  
  q<-p<-120
  set.seed(12) 
  gr<-generator(p=p, q=q, lwb=lwb, upb=upb)
  
  
  write.table(gr[[1]],file=paste("./Graph/graphGB", "120 B.csv",sep=""),row.names=F,col.names=F,sep=",")
  write.table(gr[[2]],file=paste("./Graph/graphGB", "120 P.csv",sep=""),row.names=F,col.names=F,sep=",")
  write.table(gr[[3]],file=paste("./Graph/graphGB", "120 Phi.csv",sep=""),row.names=F,col.names=F,sep=",")
  
}
