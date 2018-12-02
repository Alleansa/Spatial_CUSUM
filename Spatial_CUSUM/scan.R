##load data
load("/Users/Allen/Documents/Research/DATA/Brain_FMRI/Data.RData")
library(ggplot2)
library(grid)
library(magrittr)
library(dplyr)
yyy <- sqrt(121) * atanh(xxx)
#one_sample<-(yyy[,,13,1])
#one_sample<-(yyy[,,10,1])
#one_sample<-(yyy[,,13,6])
#one_sample<-(yyy[,,5,1])
#one_sample<-(yyy[,,15,6])
one_sample<-(yyy[,,1,1])
N.grid<-128
grid.point <- expand.grid(x=seq(1,N.grid,by=1),y=seq(1,N.grid,by=1))
location<-data.frame(element=rep(1:128^2),lon=grid.point$y,lat=grid.point$x)
write.csv(location,"Desktop/location.csv",row.names = FALSE) 

grid.point$Observed<-as.vector(one_sample)
Obs<-data.frame(element=rep(1:128^2),case=grid.point$Observed)
write.csv(Obs,"Desktop/Observed.csv",row.names = FALSE) 



grid.point$scan=0
grid.point$scan[index]=1
a<-ggplot(grid.point) + geom_tile(aes(x = x, y = y,fill=Observed)) +scale_fill_distiller(palette = "RdBu")+ xlab("")+ylab("")+ theme(legend.position = 'none',panel.grid.major =element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(),axis.text = element_blank(),axis.ticks = element_blank())
a+geom_point(data=grid.point, aes(x = x, y = y,colour=as.logical(grid.point$scan)))+
  scale_colour_manual(values = c(NA,"black"))+ 
  xlab("")+ylab("")+
  theme(legend.position='none',
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank())
