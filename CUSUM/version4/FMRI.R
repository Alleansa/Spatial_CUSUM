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
#one_sample<-(yyy[,,5,6])
#one_sample<-(yyy[,,15,6])
#one_sample<-(yyy[,,1,1])
N.grid<-128
grid.point <- expand.grid(x=seq(1,N.grid,by=1),y=seq(1,N.grid,by=1))
grid.point$Observed<-as.vector(one_sample)
#ggplot()+geom_point(data=grid.point, aes(x = x, y = y,colour=Observed))+scale_colour_continuous(limits=c(-1,1),low="green",high="red")

###cut the circle
for (i in c(1:128^2)){
  loc_x<-grid.point[i,1]
  loc_y<-grid.point[i,2]
  grid.point[i,4]<-sqrt((loc_x-64.5)^2+(loc_y-64.5)^2)
}
grid.point$inside[grid.point$V4<=64.5]=1
grid.point$inside[grid.point$V4>64.5]=0
grid.point$V4=NULL
grid.point$Observed[grid.point$inside==0]=NA
grid.point$inside=NULL
a<-ggplot(grid.point) + geom_tile(aes(x = x, y = y,fill=Observed)) +scale_fill_distiller(palette = "RdBu")+ xlab("")+ylab("")+ theme(panel.grid.major =element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(),axis.text = element_blank(),axis.ticks = element_blank())



###plot fmri
grid.point2 <- grid.point %>% mutate(bound = as.factor(bound))
a<-ggplot(grid.point2) + geom_tile(aes(x = x, y = y,fill=Observed)) +
  scale_fill_distiller(palette = "RdBu")+
  xlab("")+ylab("")+
  theme(panel.grid.major =element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(),axis.text = element_blank(),axis.ticks = element_blank())

#b<-ggplot(grid.point2)+geom_point(aes(x = x, y = y, colour=detection))+
#  scale_color_manual(values = c("white", "black"))+
#  xlab("")+ylab("")+
#  theme(legend.position='none',axis.text = element_blank(),axis.ticks = element_blank())



c<-a+  geom_point(aes(x = x, y = y, colour=bound)) +
  scale_color_manual(values = c(NA, "black"))+
  xlab("")+ylab("")+
  theme(panel.grid.major =element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(),axis.text = element_blank(),axis.ticks = element_blank())


grid.newpage()  ##newpage
pushViewport(viewport(layout = grid.layout(1,2))) ### divide into 1*3
vplayout <- function(x,y){
  viewport(layout.pos.row = x, layout.pos.col = y)
}
print(a, vp = vplayout(1,1))   ### put a in (1,1)
print(c, vp = vplayout(1,2))   ### put b in (1,2)
