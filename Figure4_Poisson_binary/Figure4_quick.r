###--------------------------------------------------------------###
### Draw Figure 4
### Last updated: June.25 2014
###--------------------------------------------------------------###
## This file contains a toy example for the experiment in Figure4.r.
## Here we generate 3 data sets for each of the 5 graphs.


## Required library
library(glmnet)

# Toy example
M1<-5 # Number of graphs
M2<-3 # Number of datasets for each graph
size<-200



###-------------------------------------------------###
## Data and Graph generation ####
## Note: Graphs and datasets for this toy example are 
## available in ./Graph and ./Data.

## Generate graphs
## p specifies the number of Poisson nodes 
## (note that there are equal numbers of Poisson and binary nodes)
source("../Sources/MGM_Graph.r")
source("./PB_Graph.r")

PB_Graph(M1=M1,lwb=0.8,upb=1,p=40);


source("../Sources/MGM_Sampler.r")
source("./PB_Data.r")
PB_Data(M1=M1,M2=M2,Gibbs.n=20, burnin=200,size=2*size,p=40);
###-------------------------------------------------###



###---------------------------------------------------###
## Neighbourhood selection ####
## Note: total is the number of tuning parameters to try 
source("../Sources/MGM_Evaluation.r")
source("../Sources/MGM_Combine.r")
source("../Sources/MGM_BIC.r")
source("../Sources/MGM_NSelect.r")
source("../Sources/MGM_misc.r")
source("./PB_Comp.r")
PB_Select(M1=M1,M2=M2,size=size,p=40);
###---------------------------------------------------###



###---------------------------------------------------###
#### Plotting ####


source("../Sources/MGM_misc.r")

# Edge estimated using the interction (AND) rule
edge_and<-process(st1="./Estimates/G", st2="PB_AND.txt",range=(1:M1) )

# Edges estimated using the union (OR) rule
edge_or<-process(st1="./Estimates//G", st2="PB_OR.txt",range=(1:M1) )

# Edges estimated using the selection rule based on true parameters 
edge_S<-process(st1="./Estimates/G", st2="PB_S.txt",range=(1:M1) )

# Edges estimated using the selection rule based on consistent estimators 
edge_bic2<-process(st1="./Estimates/G", st2="PB_BIC2.txt",range=(1:M1) )


pdf(file = paste("./plots/PB-diff.pdf",sep=""), width=8, height=8)
par(mar=c(4.2,4.2,4,4)) # Margins
plot(edge_S[,2]~edge_S[,6],type="l",xlim=c(0,400),ylim=c(0,40), lwd=4,
     xlab=" ", ylab=" ",col="gray", cex.lab=2.2,cex.axis=1.5)

lines(edge_bic2[,2]~edge_bic2[,6],lwd=4,lty="dashed")
lines(edge_and[,2]~edge_and[,6],col="black",lwd=4,lty="dotted")
lines(edge_or[,2]~edge_or[,6],col="black",lwd=4,lty="twodash")
mtext('Number of Estimated Edges', outer=T, side=1, line=-1.8,cex=2.2)
mtext("Number of Correctly Estimated Edges", side=2, line=2.2,cex=2.2)
mtext('Bernoulli-Poisson Edges', outer=T, line=-2.8,cex=2.3)
dev.off()
