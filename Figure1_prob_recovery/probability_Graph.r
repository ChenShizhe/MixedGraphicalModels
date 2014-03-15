###--------------------------------------------------------------###
### Generating graphs for Fig 1
### Last updated: Mar.13th 2014
###--------------------------------------------------------------###

source("../Sources/MGM_Graph.r")

##########################################
#### Fixed graphs
q<-p<-30
set.seed(12) 
gr<-generator(p=p, q=q,lwb=0.3, upb=0.3)


write.table(gr[[1]],file=paste("./Graph/graphGB", "30Fix B.csv",sep=""),row.names=F,col.names=F,sep=",")
write.table(gr[[2]],file=paste("./Graph/graphGB", "30Fix P.csv",sep=""),row.names=F,col.names=F,sep=",")
write.table(gr[[3]],file=paste("./Graph/graphGB", "30Fix Phi.csv",sep=""),row.names=F,col.names=F,sep=",")


q<-p<-60
set.seed(12) 
gr<-generator(p=p, q=q, lwb=0.3, upb=0.3)


write.table(gr[[1]],file=paste("./Graph/graphGB", "60Fix B.csv",sep=""),row.names=F,col.names=F,sep=",")
write.table(gr[[2]],file=paste("./Graph/graphGB", "60Fix P.csv",sep=""),row.names=F,col.names=F,sep=",")
write.table(gr[[3]],file=paste("./Graph/graphGB", "60Fix Phi.csv",sep=""),row.names=F,col.names=F,sep=",")

q<-p<-120
set.seed(12) 
gr<-generator(p=p, q=q, lwb=0.3, upb=0.3)


write.table(gr[[1]],file=paste("./Graph/graphGB", "120Fix B.csv",sep=""),row.names=F,col.names=F,sep=",")
write.table(gr[[2]],file=paste("./Graph/graphGB", "120Fix P.csv",sep=""),row.names=F,col.names=F,sep=",")
write.table(gr[[3]],file=paste("./Graph/graphGB", "120Fix Phi.csv",sep=""),row.names=F,col.names=F,sep=",")

