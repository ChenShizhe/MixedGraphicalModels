###--------------------------------------------------------------###
### Draw Figure 2
### Last updated: Mar.27th 2014
###--------------------------------------------------------------###

## Required library
library(glmnet)

# Toy example
M2<-100 # Number of datasets for each graph


source("../Sources/MGM_Graph.r")
source("./probability_Graph.r")
P_Graph(lwb=0.3,upb=0.3,p=30);
P_Graph(lwb=0.3,upb=0.3,p=60);
P_Graph(lwb=0.3,upb=0.3,p=120);


source("../Sources/MGM_Sampler.r")
source("./probability_Data.r")
P_Data(M2=M2,Gibbs.n=500,burnin=3000,p=30,size=6000);
P_Data(M2=M2,Gibbs.n=500,burnin=3000,p=60,size=6000);
P_Data(M2=M2,Gibbs.n=500,burnin=3000,p=120,size=6000);



source("../Sources/MGM_Evaluation.r")
source("../Sources/MGM_Combine.r")
source("../Sources/MGM_BIC.r")
source("../Sources/MGM_NSelect.r")
source("../Sources/MGM_misc.r")


source("./probability_recover.r")
# Run the probability recovery
P_recover(p=30,total=50,M2=M2,size=6000,step=200)
P_recover(p=60,total=50,M2=M2,size=6000,step=200)
P_recover(p=120,total=50,M2=M2,size=6000,step=200)


##############################################
#### Plotting ####



prob60<-as.matrix(read.table(file=paste("./Estimates/60/","prob60cc.txt",sep="")))


prob120<-as.matrix(read.table(file=paste("./Estimates/120/","prob120cc.txt",sep="")))


prob240<-as.matrix(read.table(file=paste("./Estimates/240/","prob240cc.txt",sep="")))

samplesize.new<-seq(from=200, to=6000, by=200)

pdf(file = "./Plots/prob-cc.pdf", width=4.5, height=4)
par(mar=c(2,2.2,2.5,2))
plot(prob60[14,]~ I(samplesize.new/log(60)/3), type="b",main=" ",xlab=" ", cex.main=2,
     cex.lab=1.5,cex.axis=1.5, ylab=" ",
     lwd=2,xlim=c(0,100),ylim=c(0,1))
lines(prob120[14,]~ I(samplesize.new/log(120)/3), type="b",lwd=2, col="green")
lines(prob240[14,]~ I(samplesize.new/log(240)/3), type="b",lwd=2, col="blue")
mtext('(a) Gaussian-Gaussian', outer=T, line=-2,cex=2.3)
dev.off()



##############################
## continuous- discrete  #####
##############################


prob60<-as.matrix(read.table(file=paste("./Estimates/60/","prob60cd.txt",sep="")))

prob120<-as.matrix(read.table(file=paste("./Estimates/120/","prob120cd.txt",sep="")))

prob240<-as.matrix(read.table(file=paste("./Estimates/240/","prob240cd.txt",sep="")))

pdf(file = "./Plots/prob-cd.pdf", width=4.5, height=4)
par(mar=c(2,2.2,2.5,2))
plot(prob60[14,]~ I(samplesize.new/log(60)/3), type="b", xlab=" ", 
     main=" ",cex.main=1.7, 
     cex.lab=1.5,cex.axis=1.5, ylab=" ", lwd=2,xlim=c(0,100),ylim=c(0,1))
lines(prob120[14,]~ I(samplesize.new/log(120)/3), type="b",lwd=2, col="green")
lines(prob240[14,]~ I(samplesize.new/log(240)/3), type="b",lwd=2, col="blue")
mtext('(b) Gaussian-Binary', outer=T, line=-2,cex=2.3)
dev.off()




##############################
## discrete - continuous #####
##############################


prob60<-as.matrix(read.table(file=paste("./Estimates/60/","prob60dc.txt",sep="")))

prob120<-as.matrix(read.table(file=paste("./Estimates/120/","prob120dc.txt",sep="")))

prob240<-as.matrix(read.table(file=paste("./Estimates/240/","prob240dc.txt",sep="")))

pdf(file = "./Plots/prob-dc.pdf", width=4.5, height=4)
par(mar=c(2,2.2,2.5,2))
plot(prob60[14,]~ I(samplesize.new/log(60)/3), type="b", xlab=" ",
     main=" ",cex.main=1.7,
     cex.lab=1.5,cex.axis=1.5,ylab=" ", lwd=2,xlim=c(100,400),ylim=c(0,1))
lines(prob120[14,]~ I(samplesize.new/log(120)/3), type="b",lwd=2, col="green")
lines(prob240[14,]~ I(samplesize.new/log(240)/3), type="b",lwd=2, col="blue")
mtext('(c) Binary-Gaussian', outer=T, line=-2,cex=2.3)
dev.off()



##############################
##   discrete - discrete #####
##############################


prob60<-as.matrix(read.table(file=paste("./Estimates/60/","prob60dd.txt",sep="")))

prob120<-as.matrix(read.table(file=paste("./Estimates/120/","prob120dd.txt",sep="")))

prob240<-as.matrix(read.table(file=paste("./Estimates/240/","prob240dd.txt",sep="")))


pdf(file = "./Plots/prob-dd.pdf", width=4.5, height=4)
par(mar=c(2,2.2,2.5,2))
plot(prob60[14,]~ I(samplesize.new/log(60)/3), type="b", xlab="",
     main=" ",cex.main=1.7,
     cex.lab=1.5,cex.axis=1.5,ylab=" ", lwd=2,xlim=c(100,400),ylim=c(0,1))
lines(prob120[14,]~ I(samplesize.new/log(120)/3), type="b",lwd=2, col="green")
lines(prob240[14,]~ I(samplesize.new/log(240)/3), type="b",lwd=2, col="blue")
mtext('(d) Binary-Binary', outer=T, line=-2,cex=2.3)
dev.off()

