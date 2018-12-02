##Apr.3 2017
#generating n random point (x,y) and observing data {e} in [-1,1]*[-1,1] 
n=200
x<-runif(n,-1,1)
y<-runif(n,-1,1)
e<-runif(n,10,20)
data<-data.frame(x=x,y=y,Observed=e,label=rep(1,n),id=seq(1,n))
#suppose the abnormal region is c((1/4,1/4),1/2) and the shift for mean is 2
data$r<-(data$x-1/2)^2+(data$y-1/2)^2
data$label[data$r<1/4]<-2
data$Observed[data$label==2]<-data$Observed[data$label==2]+10

data$r<-(data$x+1/2)^2+(data$y+1/2)^2
data$label[data$r<1/4]<-3
data$Observed[data$label==3]<-data$Observed[data$label==3]-10

#plot the region
library(ggplot2)
circleFun <- function(center = c(0,0),diameter = 1, npoints = 100){
  r = diameter / 2
  tt <- seq(0,2*pi,length.out = npoints)
  xx <- center[1] + r * cos(tt)
  yy <- center[2] + r * sin(tt)
  return(data.frame(xx = xx, yy = yy))
}
dat1 <- circleFun(c(1/2,1/2),1)
dat2 <- circleFun(c(-1/2,-1/2),1)
#geom_path will do open circles, geom_polygon will do filled circles
d<-cbind.data.frame(d,dat)
ggplot(dat,aes(x=xx,y=yy)) + geom_path(group=1)
ggplot(data, aes(x = x, y = y,colour=label)) + geom_point()

ggplot() + geom_point(data=data, aes(x = x, y = y,colour=Observed)) +
  geom_path(aes(x=xx, y=yy), data=dat1)+geom_path(aes(x=xx, y=yy), data=dat2)


ggplot()+geom_path(aes(x=xx, y=yy), data=dat)


##Apr.15 2017
#generating grid
n<-25
x<-rep(c(1:n),times=n)
y<-rep(c(1:n),each=n)
Status<-rep(0,n*n)
data<-data.frame(x=x,y=y,Status=Status)
#generating 'L' shape
data[which(data$x<=0.8*n&data$x>=.4*n&data$y<=.8*n&data$y>=.4*n),3]<-1
data[which(data$x<=0.8*n&data$x>=.6*n&data$y<=.8*n&data$y>=.6*n),3]<-0
#plot the data
library(ggplot2)
ggplot() + geom_point(data=data, aes(x = x, y = y,colour=Status))
data$Observed<-rnorm(n*n,0,1)
data[data$Status==1,]$Observed=data[data$Status==1,]$Observed+3
ggplot() + geom_point(data=data, aes(x = x, y = y,colour=Observed))

p1<-rep(0,10)
p2<-p1

for (i in c(1:10)){
  n<-20
  x<-rep(c(1:n),times=n)
  y<-rep(c(1:n),each=n)
  Status<-rep(1,n*n)
  data<-data.frame(x=x,y=y,Status=Status)
  #generating 'L' shape
  data[which(data$x<=0.8*n&data$x>=.4*n&data$y<=.8*n&data$y>=.4*n),3]<-2
  data[which(data$x<=0.8*n&data$x>=.6*n&data$y<=.8*n&data$y>=.6*n),3]<-1
  data$Observed=runif(n*n)
  data[data$Status==2,]$Observed=data[data$Status==2,]$Observed+0.5
  s<-get_sum_CUSUM(generating_tree(data))
  p1[i]<-sum(CV$c1000>s)/1000
  
}



Status<-read.csv('/Users/Allen/Documents/Study/Ph.D/Code/R/CUSUM/test.csv',header = FALSE)
Status$V1[1]=1
x<-rep(c(1:28),time=52)
y<-rep(c(1:52),each=28)
data<-data.frame(x=x,y=y,Status=Status$V1)
data$Observed=runif(1456,0,1)
data[data$Status==0,]$Observed=data[data$Status==0,]$Observed+1
ggplot() + geom_point(data=data, aes(x = y, y = rev(x),colour=Observed))
