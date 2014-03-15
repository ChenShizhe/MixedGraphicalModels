###--------------------------------------------------------------###
### Figure 3
### Last updated: Mar.13th 2014
###--------------------------------------------------------------###
### Note:
###     Estimates from GRaFo (Fellinghauer et al. 2013) are not included
###     Code for GRaFo can be found at 
###     http://www.sciencedirect.com/science/article/pii/S0167947313000789
###--------------------------------------------------------------###
source("../Sources/readplots.r")


M<-3

# Edge estimated using the interction (AND) rule
edge_and<-process(st1="./Estimates/G", st2="PB200_AND.txt",range=(1:M) )
  
# Edges estimated using the union (OR) rule
edge_or<-process(st1="./Estimates//G", st2="PB200_OR.txt",range=(1:M) )

# Edges estimated using the selection rule based on true parameters 
edge_S<-process(st1="./Estimates/G", st2="PB200_S.txt",range=(1:M) )
  
# Edges estimated using the selection rule based on consistent estimators 
edge_bic2<-process(st1="./Estimates/G", st2="PB200_BIC2.txt",range=(1:M) )



pdf(file = paste("./plots/PB-diff.pdf",sep=""), width=8, height=6)
par(mar=c(4,3.6,4,0.3)) # Margins
plot(edge_S[,2]~edge_S[,6],type="l",xlim=c(0,400),ylim=c(0,40), lwd=4,
     xlab="Number of total estimated edges", ylab=" ",col="gray", cex.lab=2.2,cex.axis=1.5)
lines(edge_bic2[,2]~edge_bic2[,6],lwd=4,lty="dashed")
lines(edge_and[,2]~edge_and[,6],col="#FF4040",lwd=4,lty="solid")
lines(edge_or[,2]~edge_or[,6],col="#FF404090",lwd=4,lty="dashed")
mtext("Number of true edges", side=2, line=2.2,cex=2.2)
mtext('Edges between Poisson and binary nodes', outer=T, line=-2.8,cex=2.3)
dev.off()

