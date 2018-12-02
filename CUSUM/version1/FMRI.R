load("/Users/Allen/Documents/Study/Ph.D/Code/R/CUSUM/Data.RData")
load("/Users/Allen/Documents/Study/Ph.D/Code/R/CUSUM/B_inte.rda")
source("/Users/Allen/Documents/Study/Ph.D/Code/R/CUSUM/FDRL.R")
source("/Users/Allen/Documents/Study/Ph.D/Code/R/CUSUM/CUSUM.R")

library(ggplot2)
library(grid)
yyy <- sqrt(121) * atanh(xxx)
# one_sample<-yyy[,,5,1]
one_sample<-yyy[,,5,10]
N.grid<-127
grid.point <- expand.grid(x=seq(0,1,by=1/N.grid),y=seq(0,1,by=1/N.grid))
grid.point$Observed<-as.vector(one_sample)
a<-ggplot()+geom_point(data=grid.point, aes(x = x, y = y,colour=Observed))+scale_color_gradientn(colours = rev(rainbow(6)))+theme(legend.position='none')
grid.point$pvalues<-pvalue_normal(grid.point$Observed)

detection_map<-region_detection(grid.point,p_thre=0.01,r=0.02, B_inte)
c<-ggplot()+geom_point(data=detection_map, aes(x = x, y = y,colour=detection))+theme(legend.position='none')
d<-ggplot()+geom_point(data=detection_map, aes(x = x, y = y,colour=p))+scale_color_gradientn(colours = rev(rainbow(6)))+theme(legend.position='none')


res<-FDRLMethod(grid.point$pvalues,grid.point$x,grid.point$y,window = 0.02,alpha=0.5)
grid.point$fdr<-res$ind
b<-ggplot()+geom_point(data=grid.point, aes(x = x, y = y,colour=fdr))+theme(legend.position='none')


grid.newpage()  ##newpage
pushViewport(viewport(layout = grid.layout(2,2))) ### divide into 1*3
vplayout <- function(x,y){
  viewport(layout.pos.row = x, layout.pos.col = y)
}
print(a, vp = vplayout(1,1))   ### put a in (1,1)
print(b, vp = vplayout(1,2))   ### put b in (1,2)
print(c, vp = vplayout(2,1))  ###put c in (1,3)
print(d, vp = vplayout(2,2))  ###put c in (1,3)