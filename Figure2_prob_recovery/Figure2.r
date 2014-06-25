###--------------------------------------------------------------###
## Draw Figure 2
## Last updated: June.25 2014
###--------------------------------------------------------------###
###
## This file contains the code to reproduce Figure 2 in the paper.
## Note that this experiment is computationally intensive in both data generation
## and fitting the models. The whole program takes about a week to finish. A toy example
## is given in Figure2_quick.r


## Required library
library(glmnet)

M2<-100 # Number of datasets for each graph
size<-6000 #  the maximum sample size
# step: the increment of sample sizes
# must be a factor of size
step<-200 

###---------------------------------------------------###
## Data and Graph generation ####

# Generate graphs
# p specifies the number of Gaussian nodes 
# (note that there are equal numbers of Gaussian and binary nodes)
source("../Sources/MGM_Graph.r")
source("./probability_Graph.r")
P_Graph(lwb=0.3,upb=0.3,p=30);
P_Graph(lwb=0.3,upb=0.3,p=60);
P_Graph(lwb=0.3,upb=0.3,p=120);

# Draw samples from the generated graphs
source("../Sources/MGM_Sampler.r")
source("./probability_Data.r")
P_Data(M2=M2,Gibbs.n=500,burnin=3000,p=30,size=size);
P_Data(M2=M2,Gibbs.n=500,burnin=3000,p=60,size=size);
P_Data(M2=M2,Gibbs.n=500,burnin=3000,p=120,size=size);
###---------------------------------------------------###



###---------------------------------------------------###
## Neighbourhood selection ####
## Note: total is the number of tuning parameters to try
source("../Sources/MGM_Evaluation.r")
source("../Sources/MGM_Combine.r")
source("../Sources/MGM_BIC.r")
source("../Sources/MGM_NSelect.r")
source("../Sources/MGM_misc.r")
source("./probability_recover.r")


P_recover(p=30,total=50,M2=M2,size=size,step=step)
P_recover(p=60,total=50,M2=M2,size=size,step=step)
P_recover(p=120,total=50,M2=M2,size=size,step=step)
###---------------------------------------------------###

###---------------------------------------------------###
#### Plotting ####


samplesize.new<-seq(from=step, to=size, by=step)



prob60<-as.matrix(read.table(file=paste("./Estimates/60/","prob60cc.txt",sep="")))

prob120<-as.matrix(read.table(file=paste("./Estimates/120/","prob120cc.txt",sep="")))

prob240<-as.matrix(read.table(file=paste("./Estimates/240/","prob240cc.txt",sep="")))

pdf(file = "./Plots/prob-cc.pdf",width=4.5, height=4.5)
par(mar=c(4,4,2.5,2.5))
plot(prob60[14,]~ I(samplesize.new/log(60)/3),type="b",main=" ",xlab=" ", cex.main=2,
     cex.lab=1.5,cex.axis=1.5, ylab=" ", lty=1, pch=21, cex=1.5,
     lwd=2.5,xlim=c(0,100),ylim=c(0,1))
lines(prob120[14,]~ I(samplesize.new/log(120)/3), type="b",lwd=2.5, col="black", cex=1.5,lty=5,pch=22)
lines(prob240[14,]~ I(samplesize.new/log(240)/3),type="b",lwd=2.5, col="black", cex=1.5,lty=3,pch=24 )

mtext('(a) Gaussian-Gaussian', outer=T, line=-2,cex=2.3)
mtext('Sample size', outer=T, side=1, line=-1.5,cex=2)
mtext('Empirical Probability', outer=T, side=2, line=-1.5,cex=2)
dev.off()


##############################
## continuous- discrete  #####
##############################


prob60<-as.matrix(read.table(file=paste("./Estimates/60/","prob60cd.txt",sep="")))

prob120<-as.matrix(read.table(file=paste("./Estimates/120/","prob120cd.txt",sep="")))

prob240<-as.matrix(read.table(file=paste("./Estimates/240/","prob240cd.txt",sep="")))

pdf(file = "./Plots/prob-cd.pdf", width=4.5, height=4.5)
par(mar=c(4,4,2.5,2.5))
plot(prob60[14,]~ I(samplesize.new/log(60)/3), type="b",main=" ",xlab=" ", cex.main=2,
     cex.lab=1.5,cex.axis=1.5, ylab=" ", lty=1, pch=21, cex=1.5,
     lwd=2.5,xlim=c(0,100),ylim=c(0,1))
lines(prob120[14,]~ I(samplesize.new/log(120)/3),  type="b",lwd=2.5, col="black", cex=1.5,lty=5,pch=22)
lines(prob240[14,]~ I(samplesize.new/log(240)/3),type="b",lwd=2.5, col="black", cex=1.5,lty=3,pch=24 )

mtext('(b) Gaussian-Bernoulli', outer=T, line=-2,cex=2.3)
mtext('Sample size', outer=T, side=1, line=-1.5,cex=2)
mtext('Empirical Probability', outer=T, side=2, line=-1.5,cex=2)
dev.off()




##############################
## discrete - continuous #####
##############################


prob60<-as.matrix(read.table(file=paste("./Estimates/60/","prob60dc.txt",sep="")))

prob120<-as.matrix(read.table(file=paste("./Estimates/120/","prob120dc.txt",sep="")))

prob240<-as.matrix(read.table(file=paste("./Estimates/240/","prob240dc.txt",sep="")))

pdf(file = "./Plots/prob-dc.pdf",  width=4.5, height=4.5)
par(mar=c(4,4,2.5,2.5))
plot(prob60[14,]~ I(samplesize.new/log(60)/3), type="b",main=" ",xlab=" ", cex.main=2,
     cex.lab=1.5,cex.axis=1.5, ylab=" ", lty=1, pch=21, cex=1.5,
     lwd=2.5,xlim=c(100,400),ylim=c(0,1))
lines(prob120[14,]~ I(samplesize.new/log(120)/3),  type="b",lwd=2.5, col="black", cex=1.5,lty=5,pch=22)
lines(prob240[14,]~ I(samplesize.new/log(240)/3), type="b",lwd=2.5, col="black", cex=1.5,lty=3,pch=24 )

mtext('(c) Bernoulli-Gaussian', outer=T, line=-2,cex=2.3)
mtext('Sample size', outer=T, side=1, line=-1.5,cex=2)
mtext('Empirical Probability', outer=T, side=2, line=-1.5,cex=2)

dev.off()

##############################
##   discrete - discrete #####
##############################


prob60<-as.matrix(read.table(file=paste("./Estimates/60/","prob60dd.txt",sep="")))

prob120<-as.matrix(read.table(file=paste("./Estimates/120/","prob120dd.txt",sep="")))

prob240<-as.matrix(read.table(file=paste("./Estimates/240/","prob240dd.txt",sep="")))

pdf(file = "./Plots/prob-dd.pdf", width=4.5, height=4.5)
par(mar=c(4,4,2.5,2.5))
plot(prob60[14,]~ I(samplesize.new/log(60)/3), type="b",main=" ",xlab=" ", cex.main=2,
     cex.lab=1.5,cex.axis=1.5, ylab=" ", lty=1, pch=21, cex=1.5,
     lwd=2.5,xlim=c(100,400),ylim=c(0,1))
lines(prob120[14,]~ I(samplesize.new/log(120)/3),type="b",lwd=2.5, col="black", cex=1.5,lty=5,pch=22)
lines(prob240[14,]~ I(samplesize.new/log(240)/3),type="b",lwd=2.5, col="black", cex=1.5,lty=3,pch=24 )

mtext('(d) Bernoulli-Bernoulli', outer=T, line=-2,cex=2.3)
mtext('Sample size', outer=T, side=1, line=-1.5,cex=2)
mtext('Empirical Probability', outer=T, side=2, line=-1.5,cex=2)

dev.off()
