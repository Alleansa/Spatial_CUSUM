multi_iter<-function(N.grid=100,shift=0,scale=0,times=10,k=5,alpha=c(0.01,0.05,0.1,0.2,0.5)){
  #generating data
  Status<-rep(0,(N.grid)^2)
  grid.point = expand.grid(x=seq(1,N.grid,by=1),y=seq(1,N.grid,by=1))
  if(scale==0){
    #nondependence
    grid.point$Observed<-rnorm((N.grid)^2)
  }else{
    #dependence
    model <- RMspheric(var=1, scale=scale)
    simu.alti <- RFsimulate(model, x=grid.point$x,y=grid.point$y,grid = FALSE)
    grid.point$Observed = simu.alti@data$variable1
  }

  grid.point$Status<-Status
  grid.point[which(grid.point$x<=0.5*N.grid&grid.point$x>=.2*N.grid&grid.point$y<=.5*N.grid&grid.point$y>=.2*N.grid),4]<-1
  grid.point[which(grid.point$x<=0.5*N.grid&grid.point$x>=.3*N.grid&grid.point$y<=.5*N.grid&grid.point$y>=.3*N.grid),4]<-0

  grid.point[which(grid.point$x<=0.9*N.grid&grid.point$x>=.6*N.grid&grid.point$y<=.9*N.grid&grid.point$y>=.6*N.grid),4]<-1
  grid.point[which(grid.point$x<=0.78*N.grid&grid.point$x>=.73*N.grid&grid.point$y<=.9*N.grid&grid.point$y>=.85*N.grid),4]<-0
  grid.point[which(grid.point$x<=0.78*N.grid&grid.point$x>=.73*N.grid&grid.point$y<=.65*N.grid&grid.point$y>=.6*N.grid),4]<-0
  grid.point[which(grid.point$x<=0.8*N.grid&grid.point$x>=.7*N.grid&grid.point$y<=.8*N.grid&grid.point$y>=.7*N.grid),4]<-0


  grid.point[grid.point$Status==1,]$Observed<-grid.point[grid.point$Status==1,]$Observed+shift

  detection_prob<-signal_detection(grid.point,k=k,times=times,N.grid=N.grid)


  x_seq<-seq(0,1,length.out = times*k^2)
  est <- lpdensity(data = detection_prob,grid=x_seq, bwselect = "imse-rot")$Estimate[,4]
  est[est<0]=0
  location<-as.numeric(names(which.min(est[x_seq>=quantile(detection_prob,0.5)])))
  h0<-c(est[1:location],seq(est[location],0,length.out=length(x_seq)-location))
  ratio<-rev(cumsum(rev(h0))/cumsum(rev(est)))
  
  FP=numeric(length(alpha))
  TP=numeric(length(alpha))
  for (i in c(1:length(alpha))){
  ret<-detection_prob>=(x_seq[which.max(ratio<=alpha[i])])
  res<-grid.point$Status-ret
  FP[i]<-length(res[res==-1])/N.grid^2
  TP[i]<-1-length(res[res==1])/N.grid^2
  }
  result<-c(FP,TP)
  return(result)

}

# # ggplot(grid.point)+geom_point(aes(x=x,y=y,colour=detection))+scale_color_continuous(low="white",high="black",limits=c(0,1))+theme(legend.position='none')
library(doParallel)
cl=makeCluster(13)
registerDoParallel(cl)
aa_05=foreach(m_time=1:50,.combine=rbind,
              .packages=c("RandomFields","foreach","lpdensity"),
              .errorhandling="remove")%dopar%
              {
                multi_iter(N.grid=100,shift=0.5,scale=0,times=50,k=10)
              }
stopCluster(cl)
saveRDS(aa_05,"/vol/data/zhuz/xinzhang/spatial_cusum/ROC_05_10.rds")
#
cl=makeCluster(13)
registerDoParallel(cl)
aa_1=foreach(m_time=1:50,.combine=rbind,
             .packages=c("RandomFields","foreach","lpdensity"),
             .errorhandling="remove")%dopar%
             {
               multi_iter(N.grid=100,shift=1,scale=0,times=50,k=10)
             }
stopCluster(cl)
saveRDS(aa_1,"/vol/data/zhuz/xinzhang/spatial_cusum/ROC_1_10.rds")
#
cl=makeCluster(13)
registerDoParallel(cl)
aa_15=foreach(m_time=1:50,.combine=rbind,
              .packages=c("RandomFields","foreach","lpdensity"),
              .errorhandling="remove")%dopar%
              {
                multi_iter(N.grid=100,shift=1.5,scale=0,times=50,k=10)
              }
stopCluster(cl)
saveRDS(aa_15,"/vol/data/zhuz/xinzhang/spatial_cusum/ROC_15_10.rds")

cl=makeCluster(13)
registerDoParallel(cl)
aa_2=foreach(m_time=1:50,.combine=rbind,
             .packages=c("RandomFields","foreach","lpdensity"),
             .errorhandling="remove")%dopar%
             {
               multi_iter(N.grid=100,shift=2,scale=0,times=50,k=10)
             }
stopCluster(cl)
saveRDS(aa_2,"/vol/data/zhuz/xinzhang/spatial_cusum/ROC_2_10.rds")

library(ggplot2)

TP<-c(0,0.087212,0.087512,0.087830,0.088630,0.092276,1,
      0,0.000112,0.000592,0.001596,0.004308,0.022444,1,
      0,0.000252,0.002764,0.007170,0.017738,0.070542,1,
      0,0.000160,0.000838,0.002040,0.005166,0.040912,1,
      0,0.000088,0.000980,0.002744,0.008038,0.065038,1,
      0,0.000816,0.008242,0.016806,0.033372,0.120576,1,
      0,0.000036,0.000460,0.001304,0.005188,0.065080,1,
      0,0.000102,0.001474,0.004076,0.012702,0.084402,1,
      0,0.001764,0.011062,0.021782,0.041224,0.139148,1,
      0,0.000022,0.000494,0.001824,0.009600,0.076818,1,
      0,0.000158,0.002416,0.007144,0.022312,0.101258,1,
      0,0.002246,0.012440,0.024090,0.044238,0.144798,1)

FP<-c(0,0.887308,0.891550,0.894872,0.900106,0.91189,1,
      0,0.883388,0.900894,0.913800,0.930068,0.959506,1,
      0,0.891346,0.925886,0.944638,0.964448,0.988878,1,
      0,0.923216,0.944148,0.956080,0.969140,0.992016,1,
      0,0.945222,0.965026,0.975270,0.986160,0.999070,1,
      0,0.932564,0.968026,0.981658,0.992842,0.999838,1,
      0,0.961350,0.974792,0.981796,0.991354,0.999780,1,
      0,0.962732,0.979110,0.988274,0.996060,0.999978,1,
      0,0.948480,0.979290,0.990770,0.997702,0.999998,1,
      0,0.972562,0.983800,0.991256,0.997966,0.999996,1,
      0,0.971946,0.987580,0.994906,0.998784,1.000000,1,
      0,0.954568,0.983404,0.993514,0.998674,1.000000,1)
roc_data<-data.frame(x=TP,y=FP,
                     signal=rep(c('signal=0.5','signal=1','signal=1.5','signal=2'),each=21),
                     neighbor=rep(rep(c('k=3','k=5','k=10'),each=7),times=4))

ggplot(roc_data)+geom_line(aes(x=TP,y=FP,colour=neighbor))+
  facet_wrap(~signal, nrow = 1, ncol = 4)+
  xlab('true positive')+ylab('false positive')+
  scale_fill_continuous(breaks = c('k=3','k=5','k=10'),
                        labels = paste("k=", c(3, 5, 10))
                        )

