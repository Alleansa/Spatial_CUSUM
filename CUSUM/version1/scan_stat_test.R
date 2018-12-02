library(boot)
library(spdep)
library(DCluster)

scan<-function(data){
data$Expected=mean(data$Observed)
#K&N's method over the centroids
mle<-calculate.mle(data, model="poisson")
knresults<-opgam(data=data, thegrid=data[,c("x","y")], alpha=.05, 
                 iscluster=kn.iscluster, fractpop=.15, R=, model="poisson", mle=mle)

#Plot all centroids and significant ones in red

#plot(data$x, data$y, main="Kulldorff and Nagarwalla's method")
#points(knresults$x, knresults$y, col="red", pch=19)

#Plot first cluster with the highest likelihood ratio test in green
clusters<-get.knclusters(data, knresults)
idx<-which.max(knresults$statistic)
#points(data$x[clusters[[idx]]], data$y[clusters[[idx]]], col="green", pch=19)
p<-min(knresults$pvalue)
return(p)
}
