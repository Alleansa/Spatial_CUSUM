multi_iter<-function(N.grid=100,shift=0,scale=0,times=10,k=5,alpha=0.05){
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
  
  #testing whether stationary
  detection_prob<-signal_detection(grid.point,k=k,times=times,N.grid=N.grid)
  
  
  x_seq<-seq(0,1,length.out = times*k^2)
  est <- lpdensity(data = detection_prob,grid=x_seq, bwselect = "imse-rot")$Estimate[,4]
  est[est<0]=0
  location<-as.numeric(names(which.min(est[x_seq>=0.5])))
  h0<-c(est[1:location],seq(est[location],0,length.out=length(x_seq)-location))
  ratio<-rev(cumsum(rev(h0))/cumsum(rev(est)))
  
  grid.point$detection<-detection_prob>=(x_seq[which.max(ratio<alpha)])
  
  res<-grid.point$Status-grid.point$detection
  wrong_report<-length(res[res==-1])
  miss_report<-length(res[res==1])
  
  ret<-c(miss_report/sum(grid.point$Status),wrong_report/sum(1-grid.point$Status),wrong_report/sum(grid.point$detection))
  return(ret)
}


library(doParallel)
cl=makeCluster(13)
registerDoParallel(cl)
k_5_05=foreach(m_time=1:100,.combine=rbind,
              .packages=c("RandomFields","foreach","lpdensity"),
              .errorhandling="remove")%dopar%
              {
                multi_iter(N.grid=100,shift=0.5,scale=0,times=90,k=5,alpha=0.05)
              }
stopCluster(cl)
# 
cl=makeCluster(13)
registerDoParallel(cl)
k_5_1=foreach(m_time=1:100,.combine=rbind,
              .packages=c("RandomFields","foreach","lpdensity"),
              .errorhandling="remove")%dopar%
              {
                multi_iter(N.grid=100,shift=1,scale=0,times=90,k=5,alpha=0.05)
              }
stopCluster(cl)
# 
cl=makeCluster(13)
registerDoParallel(cl)
k_5_15=foreach(m_time=1:100,.combine=rbind,
             .packages=c("RandomFields","foreach","lpdensity"),
             .errorhandling="remove")%dopar%
             {
               multi_iter(N.grid=100,shift=1.5,scale=0,times=90,k=5,alpha=0.05)
             }
stopCluster(cl)
# 


cl=makeCluster(13)
registerDoParallel(cl)
k_5_20=foreach(m_time=1:100,.combine=rbind,
             .packages=c("RandomFields","foreach","lpdensity"),
             .errorhandling="remove")%dopar%
             {
               multi_iter(N.grid=100,shift=2,scale=0,times=90,k=5,alpha=0.05)
             }
stopCluster(cl)



