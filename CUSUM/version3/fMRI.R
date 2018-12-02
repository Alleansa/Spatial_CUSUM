load("/Users/Allen/Documents/Research/DATA/Brain_FMRI/Data.RData")
source("/Users/Allen/Documents/Research/CODE/R/CUSUM/version3/sub_function.R")

library(ggplot2)
library(grid)
yyy <- sqrt(121) * atanh(xxx)
one_sample<-(yyy[,,13,6])
N.grid<-128
grid.point <- expand.grid(x=seq(1/N.grid,1,by=1/N.grid),y=seq(1/N.grid,1,by=1/N.grid))
grid.point$Observed<-as.vector(one_sample)
a<-ggplot()+geom_point(data=grid.point, aes(x = x, y = y,colour=Observed))+scale_color_gradientn(colours = rev(rainbow(6)))+theme(legend.position='none')
grid.point$pvalues<-pvalue_normal(grid.point$Observed)

res<-FDRLMethod(grid.point$pvalues,grid.point$x,grid.point$y,window = 0.1,alpha=0.3)
grid.point$fdr<-as.logical(res$ind)
c<-ggplot()+geom_point(data=grid.point, aes(x = x, y = y,colour=fdr))+theme(legend.position='none')


detection_map<-spatial_CUSUM(clean_data,r=0.01,k,threshold=0.4627)
d<-ggplot()+geom_point(data=detection_map, aes(x = x, y = y,colour=detection))+theme(legend.position='none')



# detection_map<-multi_region_detection(grid.point,r=0.05,k=10, threshold=1.90,M=0)
# detection_map$p<-(detection_map$detection==max(detection_map$detection))
# c<-ggplot()+geom_point(data=detection_map, aes(x = x, y = y,colour=p))+theme(legend.position='none')
# d<-ggplot(detection_map)+geom_point(aes(x = x, y = y,colour=detection))+scale_color_gradientn(colours = rev(rainbow(6)))+theme(legend.position='none')
# 
# 
# res<-FDRLMethod(grid.point$pvalues,grid.point$x,grid.point$y,window = 0.05,alpha=0.5)
# grid.point$fdr<-as.logical(res$ind)
# b<-ggplot()+geom_point(data=grid.point, aes(x = x, y = y,colour=fdr))+theme(legend.position='none')
# 
# 
# grid.newpage()  ##newpage
# pushViewport(viewport(layout = grid.layout(2,2))) ### divide into 1*3
# vplayout <- function(x,y){
#   viewport(layout.pos.row = x, layout.pos.col = y)
# }
# print(a, vp = vplayout(1,1))   ### put a in (1,1)
# print(b, vp = vplayout(1,2))   ### put b in (1,2)
# print(c, vp = vplayout(2,1))  ###put c in (1,3)
# print(d, vp = vplayout(2,2))  ###put c in (1,3)