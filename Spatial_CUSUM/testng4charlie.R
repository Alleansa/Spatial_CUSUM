library(raster)
library(ggplot2)
x<-brick("xin_test_image.tif")
x <- dropLayer(x, 1)
y <- getValues(x)[,5]
grid.point = expand.grid(x=seq(1,146,by=1),y=seq(1,108,by=1))
grid.point$Observed<-y

CUSUM<-function(e){
  e_star<-c(0,e)
  T<-cumsum(e_star)-mean(e)*c(0:length(e))
  return(T^2/length(e))
}

#locate changepoint
locate.change<-function (x) 
{
  x <- x[!is.na(x)]
  x <- as.matrix(x)
  if (dim(x)[2] == 1) 
    x <- t(x)
  n <- dim(x)[2]
  x<-x/mad(x)
  T <-CUSUM(x)
  changepoint <- which.max(T)
  cusum <- mean(T)
  return(c(changepoint+1,cusum))
}
#sample method and sum method
my_sample<-function(x){
  x<-x[!is.na(x)]
  if(length(x)>0){return(sample(x,1))}
  else{return(NA)}
}
my_mean<-function(x){
  x<-x[!is.na(x)]
  if(length(x)>0){return(mean(x))}
  else{return(NA)}
}

signal_detection<-function(grid.point,k=5,times=10){
  detection<-numeric(dim(grid.point)[1])
  for (xi in c((1-k):0)){
    for (yi in c((1-k):0)){ 
      for (i in c(1:times)){
        grid.point <- within(grid.point, {
          grp_x <- cut(x, seq(xi,146+k,by=k), labels = FALSE)
          grp_y <- cut(y, seq(yi,108+k,by=k), labels = FALSE)
        })
        grid.point$group_label<-grid.point$grp_x+(grid.point$grp_y-1)*max(grid.point$grp_x)
        group_sample <- aggregate(grid.point$Observed,by=list(grid.point$group_label),FUN=function(x){my_sample(x)})
        group_mean<-aggregate(grid.point$Observed,by=list(grid.point$group_label),FUN=my_mean)
        group<-cbind(group_sample,group_mean$x)
        names(group)<-c("group","sample","prob")
        rm(group_sample,group_mean)
        location<-locate.change(group$sample[order(-group$prob)])[1]
        signif_group<-order(-group$prob)[1:location]
        detection[grid.point$group_label %in% signif_group]=1+detection[grid.point$group_label %in% signif_group]
      }
    }
  }
  detection<-detection/(times*k^2)
  return(detection)
}


k=5
times=k^2
detection_prob<-signal_detection(grid.point,k,times)

library(grid)
grid.newpage()  ##newpage
pushViewport(viewport(layout = grid.layout(2,1))) ### divide into 1*3
vplayout <- function(x,y){
  viewport(layout.pos.row = x, layout.pos.col = y)
}
print(a, vp = vplayout(1,1))   ### put a in (1,1)
print(b, vp = vplayout(2,1))   ### put b in (1,2)
