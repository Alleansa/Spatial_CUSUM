rm(list=ls())
#computing the CUSUM for series {e}
CUSUM<-function(e){
  e_star<-c(0,e)
  T<-cumsum(e_star)-mean(e)*c(0:length(e))
  return(T^2/length(e))
}


#neigh info pick
sample_nozero<-function(I,seq,size){
  valid<-t(I*seq)
  return(sample(valid,size,prob=I/sum(I),replace = 1))
}

neigh_pick<-function(e,D,r,k){
  Idct<-D<=r
  Sample<-apply(Idct,1,sample_nozero,seq=e,size=k)
  return(Sample)
}


#main code
multi_iter<-function(N.grid=20,r=0.06,k=5,shift=0,scale=0,times=100){
  Status<-rep(0,(N.grid+1)^2)
  grid.point = expand.grid(x=seq(0,1,by=1/N.grid),y=seq(0,1,by=1/N.grid))
  
  if(scale==0){
  #nondependence
  grid.point$Observed<-rnorm((N.grid+1)^2)
  }else{
  #dependence
  model <- RMspheric(var=1, scale=scale)
  simu.alti <- RFsimulate(model, x=grid.point$x,y=grid.point$y,grid = FALSE)
  grid.point$Observed = simu.alti@data$variable1
  }
  
  
  #generating 'L' shape
  grid.point$Status<-Status
  grid.point[which(grid.point$x<=0.8&grid.point$x>=.4&grid.point$y<=.8&grid.point$y>=.4),4]<-1
  grid.point[which(grid.point$x<=0.8&grid.point$x>=.6&grid.point$y<=.8&grid.point$y>=.6),4]<-0
  grid.point[grid.point$Status==1,]$Observed<-grid.point[grid.point$Status==1,]$Observed+shift
  
  Dist_matrix <- as.matrix(dist(cbind(grid.point$x,grid.point$y),upper=1,diag=1))
  aa=foreach(m_time=1:times,.combine=rbind,
             .packages="RandomFields")%do%
             {
               one_iter(grid.point,r,k,Dist_matrix)
             }

  ret<-sum(aa)/times
  return(ret)
}


one_iter<-function(grid.point,r,k,Dist_matrix,threshold=0.4627){
  fail_count<-0
  n<-length(grid.point$Observed)
  
  neigh<-neigh_pick(grid.point$Observed,Dist_matrix,r,k+1)
  neigh_mean <-matrix(rep(1/k,k),1,k)%*%neigh[2:(k+1),]
  neigh_sample<-neigh[1,]
  #test whether there is anomaly
  test_seq<-neigh_sample[order(-neigh_mean)]
  if(mean(CUSUM(test_seq/mad(test_seq)))<=threshold) {
    #print('There is no anomaly in this region')
    fail_count<-1
    
    # miss_report<-sum(grid.point$Status)
    # aprox_fdr<-match_wrong(pvalues,grid.point$Status,grid.point$x,grid.point$y,miss_report,window=r)
    # fdr_wrong_report<-aprox_fdr[1]
    # fdr_miss_report<-aprox_fdr[2]
    
    
    # res<-FDRLMethod(pvalues,grid.point$x,grid.point$y,window=r,alpha=p_thre)
    # fdr_wrong_report<-sum(res$ind)
    }
  # }else {
    # location<-locate.change(test_seq)[1]
    # grid.point$detection<-rep(0,n)
    # grid.point$detection[order(-grid.point$nn_mean)[1:location]]<-1
    # # grid.point$boundary_score<-neigh_mean(grid.point$detection,Dist_matrix,r)
    # # grid.point$boundary<-rep(0,n)
    # # grid.point$boundary[abs(grid.point$boundary_score-0.5)<0.45]<-1
    # # grid.point$detection<-boundary_detection(grid.point$detection,grid.point$boundary,grid.point$Observed)
    # 
    # res<-grid.point$Status-grid.point$detection
    # wrong_report<-length(res[res==-1])
    # miss_report<-length(res[res==1])
    # aprox_fdr<-match_wrong(pvalues,grid.point$Status,grid.point$x,grid.point$y,miss_report,window =r)
    # fdr_wrong_report<-aprox_fdr[1]
    # fdr_miss_report<-aprox_fdr[2]
    
    # wrong_report<-sum(grid.point$detection)
    # fdrres<-FDRLMethod(pvalues,grid.point$x,grid.point$y,window=r,alpha=p_thre)
    # fdr_wrong_report<-sum(fdrres$ind)
  # }
  
  # return(c(miss_report,
  #          fail_count,
  #          wrong_report,
  #          fdr_miss_report,
  #          fdr_wrong_report))
  return(fail_count)
}


library(doParallel)
cl=makeCluster(13)
registerDoParallel(cl)
aa1=foreach(m_time=1:100,.combine=rbind,
            .packages=c("RandomFields","foreach"))%dopar%
            {
              multi_iter(N.grid=150,r=0.1,k=10,times=100)
            }
stopCluster(cl)
saveRDS(aa1,"/home/xinzhang/spatial_CUSUM/Result_H0/2/aa1.rds")

cl=makeCluster(13)
registerDoParallel(cl)
aa2=foreach(m_time=1:100,.combine=rbind,
            .packages=c("RandomFields","foreach"))%dopar%
            {
              multi_iter(N.grid=150,r=0.1,k=5,times=100)
            }
stopCluster(cl)
saveRDS(aa2,"/home/xinzhang/spatial_CUSUM/Result_H0/2/H0_2.rds")
cl=makeCluster(13)
registerDoParallel(cl)
aa3=foreach(m_time=1:100,.combine=rbind,
            .packages=c("RandomFields","foreach"))%dopar%
            {
              multi_iter(N.grid=150,r=0.1,k=15,times=100)
            }
stopCluster(cl)
saveRDS(aa3,"/home/xinzhang/spatial_CUSUM/Result_H0/2/H0_3.rds")

cl=makeCluster(13)
registerDoParallel(cl)
aa4=foreach(m_time=1:100,.combine=rbind,
            .packages=c("RandomFields","foreach"))%dopar%
            {
              multi_iter(N.grid=150,r=0.15,k=5,times=100)
            }
stopCluster(cl)
saveRDS(aa4,"/home/xinzhang/spatial_CUSUM/Result_H0/2/H0_4.rds")

cl=makeCluster(13)
registerDoParallel(cl)
aa5=foreach(m_time=1:100,.combine=rbind,
            .packages=c("RandomFields","foreach"))%dopar%
            {
              multi_iter(N.grid=150,r=0.05,k=5,times=100)
            }
stopCluster(cl)
saveRDS(aa5,"/home/xinzhang/spatial_CUSUM/Result_H0/2/H0_5.rds")


cl=makeCluster(13)
registerDoParallel(cl)
aa6=foreach(m_time=1:100,.combine=rbind,
            .packages=c("RandomFields","foreach"))%dopar%
            {
              multi_iter(N.grid=150,r=0.15,k=10,times=100)
            }
stopCluster(cl)
saveRDS(aa6,"/home/xinzhang/spatial_CUSUM/Result_H0/2/H0_6.rds")
