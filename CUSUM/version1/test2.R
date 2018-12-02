D<-get_distance(data$x,data$y)
e<-data$Observed
k<-which.max(e)
e<-e[order(D[k,])]
Btest(e)
cu<-rep(0,200)
for (i in 1:200){
  cu[i]<-get_k_CUSUM(e,i)
}
t<-which.max(cu)
my_order<-order(D[k,])
ab<-my_order[1:t]
plot(data$x,data$y)
points(data[ab,1], data[ab,2], col="red", pch=19)

ab<-e>mean(e)
