###--------------------------------------------------------------###
### Generating graphs for Poisson-Binary network
### Last updated: Mar.19th 2014
###--------------------------------------------------------------###


GB_Graph<-function(M1, lwb, upb){
  p<-q<-20
  
  set.seed(12) 
  
  #####
  for(i in 1:M1){
    gr<-generator(p=p, q=q, lwb=lwb, upb=upb)
    write.table(gr[[1]],file=paste("./Graph/graphGB", "20 B", i,".csv",sep=""),row.names=F,col.names=F,sep=",")
    write.table(gr[[2]],file=paste("./Graph/graphGB", "20 P", i,".csv",sep=""),row.names=F,col.names=F,sep=",")
    write.table(gr[[3]],file=paste("./Graph/graphGB", "20 Phi", i,".csv",sep=""),row.names=F,col.names=F,sep=",")
    
  }
  
}
