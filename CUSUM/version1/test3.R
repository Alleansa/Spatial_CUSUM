my_p=ggplot()
for (k in c(1:20)){
n<-20
x<-rep(c(1:n),times=n)
y<-rep(c(1:n),each=n)
Status<-rep(1,n*n)
data<-data.frame(x=x,y=y,Status=Status)
#generating 'L' shape
data[which(data$x<=0.8*n&data$x>=.4*n&data$y<=.8*n&data$y>=.4*n),3]<-2
data[which(data$x<=0.8*n&data$x>=.6*n&data$y<=.8*n&data$y>=.6*n),3]<-1
data$Observed=runif(n*n)
data[data$Status==2,]$Observed=data[data$Status==2,]$Observed
seq<-generating_tree(data)
sT<-rep(0,400)
for (i in 1:400){
  sT[i]<-get_k_CUSUM(seq,i)
}
newT<-data.frame(id=c(1:400),value=sT)
my_p<-my_p+geom_path(data=newT,aes(x=id,y=value))
}
