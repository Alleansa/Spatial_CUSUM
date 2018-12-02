thre<-rep(0,100)
for (t in c(1:100)){
  n<-25
  x<-rep(c(1:n),times=n)
  y<-rep(c(1:n),each=n)
  Status<-rep(0,n*n)
  data<-data.frame(x=x,y=y,Status=Status)
  #generating 'L' shape
  data[which(data$x<=0.8*n&data$x>=.4*n&data$y<=.8*n&data$y>=.4*n),3]<-1
  data[which(data$x<=0.8*n&data$x>=.6*n&data$y<=.8*n&data$y>=.6*n),3]<-0
  data$Observed<-rnorm(n*n,0,1)
  data[data$Status==1,]$Observed=data[data$Status==1,]$Observed+15
  pvalues<-pvalue_normal(data$Observed)
  res<-FDRLMethod(pvalues,data$x,data$y,window=5,alpha=0.01)
  thre[t]<-res$threshold
}