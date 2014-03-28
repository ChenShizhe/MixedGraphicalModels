###--------------------------------------------------------------###
## Draw Figure 3
## Last updated: Mar.27 2014
###--------------------------------------------------------------###
###
## This file contains a toy example for the experiment in Figure3.r.
## Here we generate 3 data sets for each of the 5 graphs.

## Required library
library(glmnet)
library(glasso)


# Toy example
M1<-5 # Number of graphs
M2<-3 # Number of datasets for each graph
size<-200

###---------------------------------------------------###
## Data and Graph generation ####
## Note: Graphs and datasets for this toy example are 
## available in ./Graph and ./Data. 


## Generate graphs
## p specifies the number of Gaussian nodes 
## (note that there are equal numbers of Gaussian and binary nodes)
source("../Sources/MGM_Graph.r")
source("./GB_Graph.r")
GB_Graph(M1=M1,lwb=0.3,upb=0.6,p=20);


## Draw samples from the generated graphs
source("../Sources/MGM_Sampler.r")
source("./GB_Data.r")
GB_Data(M1=M1,M2=M2,Gibbs.n=20, burnin=200,size=2*size,p=20);
###---------------------------------------------------###


###---------------------------------------------------###
## Neighbourhood selection ####
## Note: total is the number of tuning parameters to try 
source("../Sources/MGM_Evaluation.r")
source("../Sources/MGM_Combine.r")
source("../Sources/MGM_BIC.r")
source("../Sources/MGM_NSelect.r")
source("../Sources/MGM_misc.r")
source("./GB_Comp.r")

GB_Select(M1=M1,M2=M2,size=size,p=20,total=100);
GB_Ising(M1=M1,M2=M2,size=size,p=20,total=100);
GB_Gaussian(M1=M1,M2=M2,size=size,p=20,total=100);
GB_Glasso(M1=M1,M2=M2,size=size,p=20,total=100);
###---------------------------------------------------###



###---------------------------------------------------###
#### Plotting ####


source("../Sources/MGM_misc.r")


edgeMGM_ratio<-process(st1="./Estimates/MGM/G", st2="MGM_ratio.txt",range=(1:M1) )
edgeGR<-process(st1="./Estimates/Glasso/G", st2="GR.txt",range=(1:M1) )
edgeMB_or<-process(st1="./Estimates/MB/G", st2="MB_or.txt",range=(1:M1) )
edgeI_or<-process(st1="./Estimates/Ising/G", st2="I_or.txt",range=(1:M1) )
BIC_count<-process(st1="./Estimates/BIC/G", st2="BICcount.txt",range=(1:M1) )





pdf(file = paste("./plots/GB-diff.pdf",sep=""), width=8, height=6)
par(mfrow=c(1,1))
par(mar=c(4.2,4,4,0.3)) # Margins
plot(edgeMGM_ratio[,2]~edgeMGM_ratio[,6],type="l",xlim=c(0,150),ylim=c(0,22), lwd=4,
     xlab="Number of Estimated Edges", yaxt='n', ylab='', cex.lab=2.2,cex.axis=1.5)
lines(edgeMB_or[,2]~edgeMB_or[,6],col="#FF4040", lty="longdash", lwd=4)
lines(edgeGR[,2]~edgeGR[,6],col="#FF404090", lty="dashed",lwd=4)
lines(edgeI_or[,2]~edgeI_or[,6],col="#FF404070",  lwd=4)
points(BIC_count[2,2]~BIC_count[2,6],pch=c(17),cex=2)
axis(2,cex.axis=1.5)
mtext("Num. of Correctly Est. Edges", side=2, line=2.2,cex=2.2)
mtext('Binary-Gaussian Edges', outer=T, line=-2.8,cex=2.3)

dev.off()

## plot of edges of the same type
edgeMGMR_same<-rbind(apply(edgeMGM_ratio[,c(1,4)],1,sum), apply(edgeMGM_ratio[,c(5,8)],1,sum))
edgeMB_same<-rbind(apply(edgeMB_or[,c(1,4)],1,sum), apply(edgeMB_or[,c(5,8)],1,sum))
edgeGR_same<-rbind(apply(edgeGR[,c(1,4)],1,sum), apply(edgeGR[,c(5,8)],1,sum))
edgeI_same<-rbind(apply(edgeI_or[,c(1,4)],1,sum), apply(edgeI_or[,c(5,8)],1,sum))
BIC_same<-rbind(apply(BIC_count[,c(1,4)],1,sum), apply(BIC_count[,c(5,8)],1,sum))


pdf(file = paste("./plots/GB-same.pdf",sep=""), width=8, height=6)
par(mar=c(4.2,4,4,0.3)) # Margins
plot(edgeMGMR_same[1,]~edgeMGMR_same[2,],type="l",xlim=c(0,240),ylim=c(0,38), lwd=4,
     xlab="Number of Estimated Edges", col="black", yaxt='n', ylab='',  cex.lab=2.2,cex.axis=1.5)
lines(edgeMB_same[1,]~edgeMB_same[2,], col="#FF4040",lty="longdash", lwd=4)
lines(edgeGR_same[1,]~edgeGR_same[2,],col="#FF404090", lty="dashed", lwd=4)
lines(edgeI_same[1,]~edgeI_same[2,],col="#FF404070", lwd=4)

points(BIC_same[1,2]~BIC_same[2,2],pch=c(17),cex=2)
axis(2,cex.axis=1.5)
mtext("Num. of Correctly Est. Edges", side=2, line=2.2,cex=2.2)
mtext('Binary-Binary and Gaussian-Gaussian Edges', outer=T, line=-2.8,cex=2.3)
dev.off()



BIC_rate<-process(st1="./Estimates/BIC/G", st2="BICrate.txt",range=(1:M1) )
BIC_rate[2,] # the precision and recall rate of the BIC by node type
