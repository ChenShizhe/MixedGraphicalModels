###--------------------------------------------------------------###
### Generating graphs for Fig 2
### Last updated: Mar.27th 2014
###--------------------------------------------------------------###

source("../Sources/MGM_Graph.r")

# the magnitude of an edge potential is drawn from Unif(lwb,upb)

P_Graph<-function(lwb,upb,p){
  q=p;
  #note: p=q in this experiment
  
  set.seed(12) 
  gr<-generator(p=p, q=q,lwb=lwb, upb=upb)
  
  
  write.table(gr[[1]],file=paste("./Graph/graphGB",2*p, " B.csv",sep=""),row.names=F,col.names=F,sep=",")
  write.table(gr[[2]],file=paste("./Graph/graphGB",2*p, " P.csv",sep=""),row.names=F,col.names=F,sep=",")
  write.table(gr[[3]],file=paste("./Graph/graphGB",2*p, " Phi.csv",sep=""),row.names=F,col.names=F,sep=",")
}

  
