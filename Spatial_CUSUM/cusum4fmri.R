#computing the CUSUM for series {e}
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





#############
signal_detection<-function(grid.point,k=5,times=10,N.grid=100){
  detection<-numeric(dim(grid.point)[1])
  for (xi in c((1-k):0)){
    for (yi in c((1-k):0)){ 
      for (i in c(1:times)){
        grid.point <- within(grid.point, {
          grp_x <- cut(x, seq(xi,N.grid+k,by=k), labels = FALSE)
          grp_y <- cut(y, seq(yi,N.grid+k,by=k), labels = FALSE)
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

k=7
times=k^2
detection_prob<-signal_detection(grid.point,k,times,128)
library(lpdensity)
x_seq<-seq(0,1,length.out = times*k^2)
est <- lpdensity(data = detection_prob,grid=x_seq, bwselect = "imse-rot")$Estimate[,4]
est[est<0]=0
location<-as.numeric(names(which.min(est[x_seq>=0.5])))
h0<-c(est[1:location],seq(est[location],0,length.out=length(x_seq)-location))
ratio<-rev(cumsum(rev(h0))/cumsum(rev(est)))

ret<-detection_prob>=(x_seq[which.max(ratio<=alpha)])
ret[ret==TRUE]='Signal'
ret[ret==FALSE]='Indifference'
ggplot(grid.point) + geom_tile(aes(x = x, y = y,fill=Observed)) +
  scale_fill_distiller(palette = "RdBu")+geom_point(data=grid.point, aes(x = x, y = y,colour=ret))+
  scale_colour_manual(values = c(NA,"black"))+ 
  xlab("")+ylab("")+
  theme(legend.title=element_blank(),
        panel.grid.major =element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(),axis.text = element_blank(),axis.ticks = element_blank())
