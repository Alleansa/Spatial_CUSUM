
load("/Users/Allen/Documents/Study/Ph.D/Code/R/CUSUM/Data.RData")
library(ggplot2)
library(grid)
yyy <- sqrt(121) * atanh(xxx)

vplayout <- function(x,y){
  viewport(layout.pos.row = x, layout.pos.col = y)
}
N.grid<-127
grid.point <- expand.grid(x=seq(0,1,by=1/N.grid),y=seq(0,1,by=1/N.grid))
grid.newpage()  ##newpage
pushViewport(viewport(layout = grid.layout(3,3))) ### divide into 1*3
for (i in c(1:9)){
  one_sample<-yyy[,,i,6]
  grid.point$Observed<-as.vector(one_sample)
  a<-ggplot()+geom_point(data=grid.point, aes(x = x, y = y,colour=Observed))+scale_color_gradientn(colours = rev(rainbow(6)))+theme(legend.position='none')
  print(a, vp = vplayout(floor((i-1)/3)+1,(i-1)%%3+1) )
}
grid.newpage()  ##newpage
pushViewport(viewport(layout = grid.layout(3,3))) ### divide into 1*3
for (i in c(10:18)){
  one_sample<-yyy[,,i,6]
  grid.point$Observed<-as.vector(one_sample)
  a<-ggplot()+geom_point(data=grid.point, aes(x = x, y = y,colour=Observed))+scale_color_gradientn(colours = rev(rainbow(6)))+theme(legend.position='none')
  print(a, vp = vplayout(floor((i-1)/3)-2,(i-1)%%3+1) )
}
grid.newpage()  ##newpage
pushViewport(viewport(layout = grid.layout(3,3))) ### divide into 1*3
for (i in c(19:22)){
  one_sample<-yyy[,,i,6]
  grid.point$Observed<-as.vector(one_sample)
  a<-ggplot()+geom_point(data=grid.point, aes(x = x, y = y,colour=Observed))+scale_color_gradientn(colours = rev(rainbow(n=6)))+theme(legend.position='none')
  print(a, vp = vplayout(floor((i-1)/3)-5,(i-1)%%3+1) )
}



grid.newpage()  ##newpage
pushViewport(viewport(layout = grid.layout(3,3))) ### divide into 1*3
for (i in c(1:9)){
  one_sample<-yyy[,,13,i]
  grid.point$Observed<-as.vector(one_sample)
  a<-ggplot()+geom_point(data=grid.point, aes(x = x, y = y,colour=Observed))+scale_color_gradientn(colours = rev(rainbow(6)))+theme(legend.position='none')
  print(a, vp = vplayout(floor((i-1)/3)+1,(i-1)%%3+1) )
}
grid.newpage()  ##newpage
pushViewport(viewport(layout = grid.layout(3,3))) ### divide into 1*3
for (i in c(10:12)){
  one_sample<-yyy[,,13,i]
  grid.point$Observed<-as.vector(one_sample)
  a<-ggplot()+geom_point(data=grid.point, aes(x = x, y = y,colour=Observed))+scale_color_gradientn(colours = rev(rainbow(6)))+theme(legend.position='none')
  print(a, vp = vplayout(floor((i-1)/3)-2,(i-1)%%3+1) )
}