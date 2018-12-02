multi_iter<-function(N.grid=100,shift=0,scale=0,times=10,k=5,alpha=0.05){
  #generating data

  Status<-rep(0,(N.grid)^2)
  grid.point = expand.grid(x=seq(1,N.grid,by=1),y=seq(1,N.grid,by=1))
  if(scale==0){
    #nondependence
    grid.point$Observed<-rnorm((N.grid)^2)
  }else{
    #dependence
    model <- RMexp(var=1, scale=scale)
    simu.alti <- RFsimulate(model, x=grid.point$x,y=grid.point$y)
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

  #testing whether stationary

  #fdrl_method
  pvalues<-pvalue_normal(grid.point$Observed)
  fdrl<-FDRLMethod(pvalues,grid.point$x,grid.point$y,window=(k-1)/2,alpha)
  fdrl_Res<-grid.point$Status-fdrl$ind
  fdrl_wrong<-length(fdrl_Res[fdrl_Res==-1])
  fdrl_miss<-length(fdrl_Res[fdrl_Res==1])
  fdrl_fdr<-fdrl_wrong/fdrl$numAlt
  if(is.nan(fdrl_fdr)){fdrl_fdr=0}


  detection_prob<-signal_detection(grid.point,k=k,times=times,N.grid=N.grid)


  x_seq<-seq(0,1,length.out = times*k^2)
  est <- lpdensity(data = detection_prob,grid=x_seq, bwselect = "imse-rot")$Estimate[,4]
  est[est<0]=0
  location<-as.numeric(names(which.min(est[x_seq>=quantile(detection_prob,0.5)])))
  h0<-c(est[1:location],seq(est[location],0,length.out=length(x_seq)-location))
  ratio<-rev(cumsum(rev(h0))/cumsum(rev(est)))

  ret<-detection_prob>=(x_seq[which.max(ratio<=alpha)])
  result<-c(fdrl$ind,ret)
  return(result)

}

# ggplot(grid.point)+geom_point(aes(x=x,y=y,colour=detection))+scale_color_continuous(low="white",high="black",limits=c(0,1))+theme(legend.position='none')
library(doParallel)

set.seed(1)
cl=makeCluster(13)
registerDoParallel(cl)
s_2=foreach(m_time=1:100,.combine=rbind,
             .packages=c("RandomFields","foreach","lpdensity"),
             .errorhandling="remove")%dopar%
             {
               multi_iter(N.grid=100,shift=1,scale=0.3,times=50,k=5,alpha=0.05)
             }
stopCluster(cl)
saveRDS(s_2,"/vol/data/zhuz/xinzhang/spatial_cusum/dependence_2.rds")

set.seed(1)
cl=makeCluster(13)
registerDoParallel(cl)
s_1=foreach(m_time=1:100,.combine=rbind,
            .packages=c("RandomFields","foreach","lpdensity"),
            .errorhandling="remove")%dopar%
            {
              multi_iter(N.grid=100,shift=1,scale=0.1,times=50,k=5,alpha=0.05)
            }
stopCluster(cl)
saveRDS(s_1,"/vol/data/zhuz/xinzhang/spatial_cusum/dependence_1.rds")

set.seed(1)
cl=makeCluster(13)
registerDoParallel(cl)
s_3=foreach(m_time=1:100,.combine=rbind,
            .packages=c("RandomFields","foreach","lpdensity"),
            .errorhandling="remove")%dopar%
            {
              multi_iter(N.grid=100,shift=1,scale=0.5,times=50,k=5,alpha=0.05)
            }
stopCluster(cl)
saveRDS(s_3,"/vol/data/zhuz/xinzhang/spatial_cusum/dependence_3.rds")

library(ggplot2)
library(grid)
library(tidyr)

s1<-colMeans(readRDS("/Users/Allen/Documents/Research/CODE/R/Spatial_CUSUM/result2_image/dependence_1.rds"))
fdrl_1<-s1[1:10000]
cusum_1<-s1[10001:20000]
s2<-colMeans(readRDS("/Users/Allen/Documents/Research/CODE/R/Spatial_CUSUM/result2_image/dependence_2.rds"))
fdrl_2<-s2[1:10000]
cusum_2<-s2[10001:20000]
s3<-colMeans(readRDS("/Users/Allen/Documents/Research/CODE/R/Spatial_CUSUM/result2_image/dependence_3.rds"))
fdrl_3<-s3[1:10000]
cusum_3<-s3[10001:20000]



library(dplyr)
library(ggplot2)

plot_data <- grid.point %>% select(x, y) %>%
  mutate(`FDRL with dependence scale=0.1` = fdrl_1,
         `SCUSUM with dependence scale=0.1` = cusum_1,
         `FDRL with dependence scale=0.3` = fdrl_2,
         `SCUSUM with dependence scale=0.3` = cusum_2,
         `FDRL with dependence scale=0.5` = fdrl_3,
         `SCUSUM with dependence scale=0.5` = cusum_3)

all_data <- plot_data %>% gather(title, probability, -x, -y)

ggplot(all_data, aes(x = x, y = y, color = probability)) +
  geom_point() +
  scale_colour_continuous(limits=c(0,1),low="white",high="black") +
  facet_wrap(~title, nrow = 2, ncol = 3)  +
  xlab("")+ylab("")+
  theme(panel.grid.major =element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(),axis.text = element_blank(),axis.ticks = element_blank(),legend.title=element_blank())


s_1<-readRDS("/Users/Allen/Documents/Research/CODE/R/Spatial_CUSUM/result2_image/dependence_1.rds")
fdrl_1<-mean(s_1[,1:10000]%*%(1-grid.point$Status)/rowSums(s_1[,1:10000]))
cusum_1<-mean(s_1[,10001:20000]%*%(1-grid.point$Status)/rowSums(s_1[,10001:20000]))

s_2<-readRDS("/Users/Allen/Documents/Research/CODE/R/Spatial_CUSUM/result2_image/dependence_2.rds")
fdrl_2<-mean(s_2[,1:10000]%*%(1-grid.point$Status)/rowSums(s_2[,1:10000]))
cusum_2<-mean(s_2[,10001:20000]%*%(1-grid.point$Status)/rowSums(s_2[,10001:20000]))

s_3<-readRDS("/Users/Allen/Documents/Research/CODE/R/Spatial_CUSUM/result2_image/dependence_3.rds")
fdrl_3<-mean(s_3[,1:10000]%*%(1-grid.point$Status)/rowSums(s_3[,1:10000]),na.rm=1)
cusum_3<-mean(s_3[,10001:20000]%*%(1-grid.point$Status)/rowSums(s_3[,10001:20000]))

