CUSUM<-function(e){
  e_star<-c(0,e)
  T<-cumsum(e_star)-mean(e)*c(0:length(e))
  return(T^2/length(e))
}
p<-rep(0,100)
for (k in c(1:100)){
count=0
for(j in c(1:100)){
m<-rep(0,1000)
v<-rep(0,1000)
# 
# set<-rnorm(100)
for ( i in c(1:1000)){
  m[i]<-mean(rnorm(10))
  v[i]<-rnorm(1)
}
count<-(mean(CUSUM(v[order(m)]/mad(v)))>0.4627)+count
  # count<-(mean(CUSUM(rnorm(10000)))>0.4627)+count
}
p[k]=count
}
# save(p , file="~/desktop/norm.rda")
