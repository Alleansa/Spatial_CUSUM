rm(list=ls())
#computing the CUSUM for series {e}
CUSUM<-function(e){
  e_star<-c(0,e)
  T<-cumsum(e_star)-mean(e)*c(0:length(e))
  return(T^2/length(e))
}

#get dist matrix
dist_matrix<-function(N.grid){
  grid.point = expand.grid(x=seq(0,1,by=1/N.grid),y=seq(0,1,by=1/N.grid))
  return(as.matrix(dist(cbind(grid.point[,1:2]),upper=1,diag=1)))
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

#locate changepoint
locate.change<-function (x) 
{
  x <- as.matrix(x)
  if (dim(x)[2] == 1) 
    x <- t(x)
  n <- dim(x)[2]
  x<-x/mad(x)
  T <-CUSUM(x)
  changepoint <- which.max(T)
  cusum <- mean(T)
  return(c(changepoint,cusum))
}

#main
multi_iter<-function(N.grid=20,r=0.06,k=1,shift=0,scale=0,times=100,Dist_matrix){
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
  
  aa=foreach(m_time=1:times,.combine=rbind,
             .packages="RandomFields")%do%
             {
               one_iter(grid.point,r,k,Dist_matrix)
             }
  if(times==1){
    aa<-t(aa)
  }
  
  ret<-round(c(mean(aa[,1])/sum(grid.point$Status),
               sum(aa[,2])/times,
               sum(aa[,3])/((times-sum(aa[,2]))*((N.grid+1)^2-sum(grid.point$Status))),
               mean((aa[,1]+aa[,3])/(N.grid+1)^2)),3)
  
  return(ret)
}


one_iter<-function(grid.point,r,k,Dist_matrix,threshold=0.4627){
  fail_count<-0
  miss_report<-0
  wrong_report<-0
  
  n<-length(grid.point$Observed)
  
  neigh<-neigh_pick(grid.point$Observed,Dist_matrix,r,k+1)
  neigh_mean <-matrix(rep(1/k,k),1,k)%*%neigh[2:(k+1),]
  neigh_sample<-neigh[1,]
  #test whether there is anomaly
  test_seq<-neigh_sample[order(-neigh_mean)]
  
  if(mean(CUSUM(test_seq/mad(test_seq)))<=threshold) {
    fail_count<-1
    miss_report<-sum(grid.point$Status)
  }else {
    location<-locate.change(test_seq)[1]
    grid.point$detection<-rep(0,n)
    grid.point$detection[order(-neigh_mean)[1:location]]<-1
    
    res<-grid.point$Status-grid.point$detection
    wrong_report<-length(res[res==-1])
    miss_report<-length(res[res==1])
  }
  
  return(c(miss_report,
           fail_count,
           wrong_report))
}

D<-dist_matrix(N.grid=150)

#N.grid=30
library(doParallel)

cl=makeCluster(13)
registerDoParallel(cl)
aa_0=foreach(m_time=1:100,.combine=rbind,
             .packages=c("RandomFields","foreach"))%dopar%
             {
               multi_iter(N.grid=150,r=0.01,k=5,shift=3,scale=0,times=10,Dist_matrix = D)
             }
stopCluster(cl)
saveRDS(aa_0, "/home/xinzhang/spatial_CUSUM/Result_R/2/aa_0.rds")

cl=makeCluster(13)
registerDoParallel(cl)
aa_1=foreach(m_time=1:100,.combine=rbind,
             .packages=c("RandomFields","foreach"))%dopar%
             {
               multi_iter(N.grid=100,r=0.03,k=5,shift=3,scale=0,times=10,Dist_matrix = D)
             }
stopCluster(cl)
saveRDS(aa_1, "/home/xinzhang/spatial_CUSUM/Result_R/2/aa_1.rds")

cl=makeCluster(13)
registerDoParallel(cl)
aa_2=foreach(m_time=1:100,.combine=rbind,
             .packages=c("RandomFields","foreach"))%dopar%
             {
               multi_iter(N.grid=100,r=0.05,k=5,shift=3,scale=0,times=10,Dist_matrix = D)
             }
stopCluster(cl)
saveRDS(aa_2, "/home/xinzhang/spatial_CUSUM/Result_R/2/aa_2.rds")

cl=makeCluster(13)
registerDoParallel(cl)
aa_3=foreach(m_time=1:100,.combine=rbind,
             .packages=c("RandomFields","foreach"))%dopar%
             {
               multi_iter(N.grid=100,r=0.08,k=5,shift=3,scale=0,times=10,Dist_matrix = D)
             }
stopCluster(cl)
saveRDS(aa_3, "/home/xinzhang/spatial_CUSUM/Result_R/2/aa_3.rds")

cl=makeCluster(13)
registerDoParallel(cl)
aa_4=foreach(m_time=1:100,.combine=rbind,
             .packages=c("RandomFields","foreach"),
             .errorhandling="remove")%dopar%
             {
               multi_iter(N.grid=100,r=0.1,k=5,shift=5,scale=0,times=10,Dist_matrix = D)
             }
stopCluster(cl)
saveRDS(aa_4, "/home/xinzhang/spatial_CUSUM/Result_R/2/aa_4.rds")

cl=makeCluster(13)
registerDoParallel(cl)
aa_5=foreach(m_time=1:100,.combine=rbind,
             .packages=c("RandomFields","foreach"),
             .errorhandling="remove")%dopar%
             {
               multi_iter(N.grid=100,r=0.15,k=5,shift=5,scale=0,times=10,Dist_matrix = D)
             }
stopCluster(cl)
saveRDS(aa_5, "/home/xinzhang/spatial_CUSUM/Result_R/2/aa_5.rds")


# p<-ggplot(data=da)+geom_line(aes(x=r,y=V4))+xlim(0.03,0.15)+ylim(0,0.15)+labs(x = "radius",y="misclassification rate")
