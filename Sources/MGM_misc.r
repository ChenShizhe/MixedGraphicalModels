###--------------------------###
## Additional functions
## Last updated: Mar.14th 2014
###-------------------------------------###
## This file contains the following functions:
##  1) process: reading the results for all graph
##  2) is_same: testing whether the sample is degenerate


###-------------------------------------###
## process ####
## reading the results for all graphs
## Input: 
##       st1: the first half of the string
##       st2: the second half of the string
##       range: the indices of the graphs
## Output: 
##       The average performance
###-------------------------------------###
process<-function(st1, st2, range=(1:100)){
  temp<-read.table( file=paste(st1, 1, st2,sep=""))
  a<-dim(temp)
  res<-array(0,c(a,length(range)))
  for(i in range){
    res[,,i]<-as.matrix(read.table(file=paste(st1, i, st2,sep="")))
  }  
  apply(res,c(1,2),mean,na.rm=T)
}


###-------------------------------------###
## is_same ####
## Testing whether the sample is degenerate.
## Input: 
##       a data set
## Output: 
##       True of False
###-------------------------------------###
is_same<-function(dat){
  std<-apply(dat,2,sd) 
  return( (min(std)==0))  
}
