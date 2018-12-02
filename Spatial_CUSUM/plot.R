# multi_iter<-function(N.grid=100,shift=0,scale=0,times=10,k=5,alpha=0.05){
#   #generating data
#   
#   Status<-rep(0,(N.grid)^2)
#   grid.point = expand.grid(x=seq(1,N.grid,by=1),y=seq(1,N.grid,by=1))
#   if(scale==0){
#     #nondependence
#     grid.point$Observed<-rnorm((N.grid)^2)
#   }else{
#     #dependence
#     model <- RMspheric(var=1, scale=scale)
#     simu.alti <- RFsimulate(model, x=grid.point$x,y=grid.point$y,grid = FALSE)
#     grid.point$Observed = simu.alti@data$variable1
#   }
#   
#   grid.point$Status<-Status
#   grid.point[which(grid.point$x<=0.5*N.grid&grid.point$x>=.2*N.grid&grid.point$y<=.5*N.grid&grid.point$y>=.2*N.grid),4]<-1
#   grid.point[which(grid.point$x<=0.5*N.grid&grid.point$x>=.3*N.grid&grid.point$y<=.5*N.grid&grid.point$y>=.3*N.grid),4]<-0
#   
#   grid.point[which(grid.point$x<=0.9*N.grid&grid.point$x>=.6*N.grid&grid.point$y<=.9*N.grid&grid.point$y>=.6*N.grid),4]<-1
#   grid.point[which(grid.point$x<=0.78*N.grid&grid.point$x>=.73*N.grid&grid.point$y<=.9*N.grid&grid.point$y>=.85*N.grid),4]<-0
#   grid.point[which(grid.point$x<=0.78*N.grid&grid.point$x>=.73*N.grid&grid.point$y<=.65*N.grid&grid.point$y>=.6*N.grid),4]<-0
#   grid.point[which(grid.point$x<=0.8*N.grid&grid.point$x>=.7*N.grid&grid.point$y<=.8*N.grid&grid.point$y>=.7*N.grid),4]<-0
#   
#   
#   grid.point[grid.point$Status==1,]$Observed<-grid.point[grid.point$Status==1,]$Observed+shift
#   
#   #testing whether stationary
#   
#   #fdrl_method
#   pvalues<-pvalue_normal(grid.point$Observed)
#   fdrl<-FDRLMethod(pvalues,grid.point$x,grid.point$y,window=(k-1)/2,alpha)
#   fdrl_Res<-grid.point$Status-fdrl$ind
#   fdrl_wrong<-length(fdrl_Res[fdrl_Res==-1])
#   fdrl_miss<-length(fdrl_Res[fdrl_Res==1])
#   fdrl_fdr<-fdrl_wrong/fdrl$numAlt
#   if(is.nan(fdrl_fdr)){fdrl_fdr=0}
#   
#   
#   detection_prob<-signal_detection(grid.point,k=k,times=times,N.grid=N.grid)
#   
#   
#   x_seq<-seq(0,1,length.out = times*k^2)
#   est <- lpdensity(data = detection_prob,grid=x_seq, bwselect = "imse-rot")$Estimate[,4]
#   est[est<0]=0
#   location<-as.numeric(names(which.min(est[x_seq>=quantile(detection_prob,0.5)])))
#   h0<-c(est[1:location],seq(est[location],0,length.out=length(x_seq)-location))
#   ratio<-rev(cumsum(rev(h0))/cumsum(rev(est)))
#   
#   ret<-detection_prob>=(x_seq[which.max(ratio<=alpha)])
#   res<-grid.point$Status-ret
#   wrong_report<-length(res[res==-1])
#   miss_report<-length(res[res==1])
#   
#   
#   grid.point$detection<-ret
#   
#   result<-c(fdrl$ind,ret)
#   return(result)
#   
# }
# 
# # ggplot(grid.point)+geom_point(aes(x=x,y=y,colour=detection))+scale_color_continuous(low="white",high="black",limits=c(0,1))+theme(legend.position='none')
# library(doParallel)
# cl=makeCluster(13)
# registerDoParallel(cl)
# aa_05=foreach(m_time=1:100,.combine=rbind,
#               .packages=c("RandomFields","foreach","lpdensity"),
#               .errorhandling="remove")%dopar%
#               {
#                 multi_iter(N.grid=100,shift=0.5,scale=0,times=50,k=5,alpha=0.05)
#               }
# stopCluster(cl)
# saveRDS(aa_05,"/vol/data/zhuz/xinzhang/spatial_cusum/ret_05.rds")
# # 
# cl=makeCluster(13)
# registerDoParallel(cl)
# aa_08=foreach(m_time=1:100,.combine=rbind,
#               .packages=c("RandomFields","foreach","lpdensity"),
#               .errorhandling="remove")%dopar%
#               {
#                 multi_iter(N.grid=100,shift=0.8,scale=0,times=50,k=5,alpha=0.05)
#               }
# stopCluster(cl)
# saveRDS(aa_08,"/vol/data/zhuz/xinzhang/spatial_cusum/ret_08.rds")
# # 
# cl=makeCluster(13)
# registerDoParallel(cl)
# aa_1=foreach(m_time=1:100,.combine=rbind,
#              .packages=c("RandomFields","foreach","lpdensity"),
#              .errorhandling="remove")%dopar%
#              {
#                multi_iter(N.grid=100,shift=1,scale=0,times=50,k=5,alpha=0.05)
#              }
# stopCluster(cl)
# saveRDS(aa_1,"/vol/data/zhuz/xinzhang/spatial_cusum/ret_1.rds")
# # 
# cl=makeCluster(13)
# registerDoParallel(cl)
# aa_15=foreach(m_time=1:100,.combine=rbind,
#               .packages=c("RandomFields","foreach","lpdensity"),
#               .errorhandling="remove")%dopar%
#               {
#                 multi_iter(N.grid=100,shift=1.5,scale=0,times=50,k=5,alpha=0.05)
#               }
# stopCluster(cl)
# saveRDS(aa_15,"/vol/data/zhuz/xinzhang/spatial_cusum/ret_15.rds")
# 
# cl=makeCluster(13)
# registerDoParallel(cl)
# aa_2=foreach(m_time=1:100,.combine=rbind,
#              .packages=c("RandomFields","foreach","lpdensity"),
#              .errorhandling="remove")%dopar%
#              {
#                multi_iter(N.grid=100,shift=2,scale=0,times=50,k=5,alpha=0.05)
#              }
# stopCluster(cl)
# saveRDS(aa_2,"/vol/data/zhuz/xinzhang/spatial_cusum/ret_2.rds")

library(ggplot2)
library(grid)

aa_05<-colMeans(readRDS("/Users/Allen/Documents/Research/CODE/R/Spatial_CUSUM/result2_image/ret_05.rds"))
fdrl_1<-aa_05[1:10000]
cusum_1<-aa_05[10001:20000]
#scan_1<-colMeans(readRDS("/Users/Allen/Documents/Research/CODE/R/Spatial_CUSUM/result2_image/scan_05.rds"))

aa_08<-colMeans(readRDS("/Users/Allen/Documents/Research/CODE/R/Spatial_CUSUM/result2_image/ret_08.rds"))
fdrl_2<-aa_08[1:10000]
cusum_2<-aa_08[10001:20000]
#scan_2<-colMeans(readRDS("/Users/Allen/Documents/Research/CODE/R/Spatial_CUSUM/result2_image/scan_08.rds"))

aa_10<-colMeans(readRDS("/Users/Allen/Documents/Research/CODE/R/Spatial_CUSUM/result2_image/ret_1.rds"))
fdrl_3<-aa_10[1:10000]
cusum_3<-aa_10[10001:20000]
#scan_3<-colMeans(readRDS("/Users/Allen/Documents/Research/CODE/R/Spatial_CUSUM/result2_image/scan_1.rds"))

aa_15<-colMeans(readRDS("/Users/Allen/Documents/Research/CODE/R/Spatial_CUSUM/result2_image/ret_15.rds"))
fdrl_4<-aa_15[1:10000]
cusum_4<-aa_15[10001:20000]
#scan_4<-colMeans(readRDS("/Users/Allen/Documents/Research/CODE/R/Spatial_CUSUM/result2_image/scan_15.rds"))

aa_20<-colMeans(readRDS("/Users/Allen/Documents/Research/CODE/R/Spatial_CUSUM/result2_image/ret_2.rds"))
fdrl_5<-aa_20[1:10000]
cusum_5<-aa_20[10001:20000]
#scan_5<-colMeans(readRDS("/Users/Allen/Documents/Research/CODE/R/Spatial_CUSUM/result2_image/scan_2.rds"))


library(tidyr)
library(ggplot2)

# plot_data <- grid.point %>% select(x, y) %>%
#   mutate(`FDRL with signal=0.5` = fdrl_1,
#          `Scan Stat. with signal=0.5` = scan_1,
#          `SCUSUM with signal=0.5` = cusum_1,
#          `FDRL with signal=0.8` = fdrl_2,
#          `Scan Stat. with signal=0.8` = scan_2,
#          `SCUSUM with signal=0.8` = cusum_2,
#          `FDRL with signal=1` = fdrl_3,
#          `Scan Stat. with signal=1` = scan_3,
#          `SCUSUM with signal=1` = cusum_3,
#          `FDRL with signal=1.5` = fdrl_4,
#          `Scan Stat. with signal=1,5` = scan_4,
#          `SCUSUM with signal=1.5` = cusum_4,
#          `FDRL with signal=2` = fdrl_5,
#          `Scan Stat. with signal=2` = scan_5,
#          `SCUSUM with signal=2` = cusum_5)

plot_data <- grid.point %>% select(x, y) %>%
  mutate(`FDRL with signal=0.5` = fdrl_1,
         `SCUSUM with signal=0.5` = cusum_1,
         `FDRL with signal=0.8` = fdrl_2,
         `SCUSUM with signal=0.8` = cusum_2,
         `FDRL with signal=1` = fdrl_3,
         `SCUSUM with signal=1` = cusum_3,
         `FDRL with signal=1.5` = fdrl_4,
         `SCUSUM with signal=1.5` = cusum_4,
         `FDRL with signal=2` = fdrl_5,
         `SCUSUM with signal=2` = cusum_5)

all_data <- plot_data %>% gather(title, probability, -x, -y) 

ggplot(all_data, aes(x = x, y = y, color = probability)) + 
  geom_point() + 
  scale_colour_continuous(limits=c(0,1),low="white",high="black") +
  facet_wrap(~title, nrow = 3, ncol = 5)  +
  xlab("")+ylab("")+
  theme(panel.grid.major =element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(),axis.text = element_blank(),axis.ticks = element_blank(),legend.title=element_blank())
