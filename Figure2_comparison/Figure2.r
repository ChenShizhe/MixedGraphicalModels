###--------------------------------------------------------------###
### Draw Figure 2
### Last updated: Mar.19th 2014
###--------------------------------------------------------------###
### Note: 
###     1)  To change the parameters (p,q, etc.), modify the parameters in 
###         corresponding functions (need to edit)
###--------------------------------------------------------------###

## Required library
library(glmnet)
library(glasso)


# Toy example
M1<-5 # Number of graphs
M2<-3 # Number of datasets for each graph
size<-200

source("../Sources/MGM_Graph.r")
source("./GB_Graph.r")

GB_Graph(M1=M1);


source("../Sources/MGM_Sampler.r")
source("./GB_Data.r")
GB_Data(M1=M1,M2=M2);


#To replicate the experiment in the paper, use:
# M1 <- 100
# M2 <- 20
# size <- 200
# GB_Graph(M1=M1);
# GB_Data(M1=M1,M2=M2,Gibbs.n=500, burnin=3000);




source("../Sources/MGM_Evaluation.r")
source("../Sources/MGM_Combine.r")
source("../Sources/MGM_BIC.r")
source("../Sources/MGM_NSelect.r")
source("../Sources/MGM_misc.r")
source("./GB_Comp.r")

GB_Select(M1=M1,M2=M2,size=size);

GB_Ising(M1=M1,M2=M2,size=size);
GB_Gaussian(M1=M1,M2=M2,size=size);
GB_Glasso(M1=M1,M2=M2,size=size);


##############################################
#### Plotting ####


source("../Sources/MGM_misc.r")


edgeMGM_ratio<-process(st1="./Estimates/MGM/G", st2="MGM_ratio.txt",range=(1:M1) )
edgeGR<-process(st1="./Estimates/Glasso/G", st2="GR.txt",range=(1:M1) )
edgeMB_or<-process(st1="./Estimates/MB/G", st2="MB_or.txt",range=(1:M1) )
edgeI_or<-process(st1="./Estimates/Ising/G", st2="I_or.txt",range=(1:M1) )
BIC_count<-process(st1="./Estimates/BIC/G", st2="BICcount200.txt",range=(1:M1) )
BIC_rate<-process(st1="./Estimates/BIC/G", st2="BICrate200.txt",range=(1:M1) )
#Precision and recall rates
BIC_rate 



pdf(file = paste("./plots/GB-diff.pdf",sep=""), width=8, height=6)
par(mfrow=c(1,1))
par(mar=c(4,3.6,4,0.3)) # Margins
plot(edgeMGM_ratio[,2]~edgeMGM_ratio[,6],type="l",xlim=c(0,150),ylim=c(0,22), lwd=4,
     xlab="Number of total estimated edges", yaxt='n', ylab='', cex.lab=2.2,cex.axis=1.5)
lines(edgeMB_or[,2]~edgeMB_or[,6],col="#FF4040", lty="longdash", lwd=4)
lines(edgeGR[,2]~edgeGR[,6],col="#FF404090", lty="dashed",lwd=4)
lines(edgeI_or[,2]~edgeI_or[,6],col="#FF404070",  lwd=4)
points(BIC_count[2,2]~BIC_count[2,6],pch=c(17),cex=2)
axis(2,cex.axis=1.5)
mtext("Number of true edges", side=2, line=2.2,cex=2.2)
mtext('Edges between Gaussian and binary nodes', outer=T, line=-2.8,cex=2.3)
dev.off()

## plot of edges of the same type
edgeMGMR_same<-rbind(apply(edgeMGM_ratio[,c(1,4)],1,sum), apply(edgeMGM_ratio[,c(5,8)],1,sum))
edgeMB_same<-rbind(apply(edgeMB_or[,c(1,4)],1,sum), apply(edgeMB_or[,c(5,8)],1,sum))
edgeGR_same<-rbind(apply(edgeGR[,c(1,4)],1,sum), apply(edgeGR[,c(5,8)],1,sum))
edgeI_same<-rbind(apply(edgeI_or[,c(1,4)],1,sum), apply(edgeI_or[,c(5,8)],1,sum))
BIC_same<-rbind(apply(BIC_count[,c(1,4)],1,sum), apply(BIC_count[,c(5,8)],1,sum))


pdf(file = paste("./plots/GB-same.pdf",sep=""), width=8, height=6)

par(mar=c(4,3.6,4,0.3)) # Margins
plot(edgeMGMR_same[1,]~edgeMGMR_same[2,],type="l",xlim=c(0,240),ylim=c(0,38), 
     xlab="Number of total estimated edges", col="black", lwd=4, yaxt='n', ylab='', cex.lab=2.2,cex.axis=1.5)
lines(edgeMB_same[1,]~edgeMB_same[2,], col="#FF4040",lty="longdash", lwd=4)
lines(edgeGR_same[1,]~edgeGR_same[2,],col="#FF404090", lty="dashed", lwd=4)
lines(edgeI_same[1,]~edgeI_same[2,],col="#FF404070", lwd=4)
points(BIC_same[1,2]~BIC_same[2,2],pch=c(17),cex=2)
axis(2,cex.axis=1.5)
mtext("Number of true edges", side=2, line=2.2,cex=2.2)
mtext('Edges between same types of nodes', outer=T, line=-2.8,cex=2.3)
dev.off()
