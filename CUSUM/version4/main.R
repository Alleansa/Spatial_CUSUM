rm(list=ls())

pstar_neigh<-function(e,D,r){
  n<-length(e)
  k_neigh<-rep(0,n)
  for (i in 1:n){
    # k_neigh[i]<-get_k_neigh_med(e[order(D[i,])],k)
    tem<-(D[i,]<=r)*e
    k_neigh[i]<-median(tem[tem!=0])
  }
  return(k_neigh)
}

pvalue_normal<-function(Observed){
  
  pv<-1-pnorm(Observed)
  return(pv)
}

FDRLMethod <- function(pvalues,x ,y, window, alpha, nstep = 200, lambda = 0.1){
  
  if( lambda < 0.0 ) stop("lambda must be [0,1]")
  
  tol <- 1.5e-8
  
  m <- length(pvalues)
  total_num<-m
  Dist_matrix <- as.matrix(dist(cbind(x,y),upper=1,diag=1))
  pstar<-pstar_neigh(pvalues,Dist_matrix,window)
  
  #--------------------------------------------------------------------------#
  # Number of non-rejections with a threshold of lambda                      #
  # p > lambda                                                               #
  #--------------------------------------------------------------------------#
  Wlambda <- sum( (pstar - lambda) > tol )
  
  #--------------------------------------------------------------------------#
  # Denominator of Ghat Eq 3.3                                               #
  # tst1 = p >= 0.5                                                          #
  # tst2 = p >  0.5                                                          #
  # tst1 - tst2 = p == 0.5                                                   #
  # 2(p>0.5) + (p==0.5) = 2*tst2 + (tst1 - tst2) = tst2 + tst1               #
  #--------------------------------------------------------------------------#
  tst1 <- sum((pstar - 0.5) > -tol)
  tst2 <- sum((pstar - 0.5) >  tol)
  denom <- 1.0/(tst2 + tst1)
  
  
  #--------------------------------------------------------------------------#
  # Ghat for threshold lambda                                                #
  #--------------------------------------------------------------------------#
  if( (lambda - 0.5) < tol ) {
    gHatLambda <- sum((pstar - (1.0 - lambda)) > -tol) * denom
  } else {
    gHatLambda <- 1.0 - sum((pstar - lambda) > tol) * denom
  }
  
  tVec <- (1L:nstep)/nstep
  
  thresh <- 0.0
  
  lim <- floor(nstep/2)
  
  for(i in 1L:lim) {
    
    #----------------------------------------------------------------------#
    # Number of rejections with a threshold of t                           #
    # p <= t                                                               #
    #----------------------------------------------------------------------#
    Rt <- max(1.0, sum(pstar <= (tVec[i]+tol)))
    
    #----------------------------------------------------------------------#
    # Ghat for threshold t <= 0.5                                          #
    #----------------------------------------------------------------------#
    gHatt <- sum( (pstar - (1.0-tVec[i])) > - tol) * denom
    
    if( ( Wlambda*gHatt / (Rt*(1.0-gHatLambda)) ) < alpha ) thresh <- i/nstep
    
  }
  
  for(i in (lim+1L):nstep) {
    
    #----------------------------------------------------------------------#
    # Number of rejections with a threshold of t                           #
    # p <= t                                                               #
    #----------------------------------------------------------------------#
    Rt <- max(1.0, sum(pstar <= (tVec[i]+tol)))
    
    #----------------------------------------------------------------------#
    # Ghat for threshold t > 0.5                                           #
    #----------------------------------------------------------------------#
    gHatt <- 1.0 - sum( (pstar - tVec[i]) > tol) * denom
    
    if( ( Wlambda*gHatt / (Rt*(1.0-gHatLambda)) ) < alpha ) 
      thresh <- i/nstep
  }
  
  testresult <- as.numeric(pstar <= thresh)
  
  numAlt <- sum(testresult)
  propAlt <- numAlt/m
  
  return( list("ind" = testresult, 
               "threshold" = thresh, 
               "numAlt" = numAlt, 
               "propAlt" = propAlt) )
  
}
#computing the CUSUM for series {e}
CUSUM<-function(e){
  e_star<-c(0,e)
  T<-cumsum(e_star)-mean(e)*c(0:length(e))
  return(T^2/length(e))
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

#signal detection
signal_detection<-function(grid.point,times=500){
  detection<-numeric(dim(grid.point)[1])
  N.grid=100
  for (i in c(1:times)){
    k<-50
    grid.point <- within(grid.point, {
      grp_x <- cut(x, seq(0,N.grid,by=N.grid/k), labels = FALSE)
      grp_y <- cut(y, seq(0,N.grid,by=N.grid/k), labels = FALSE)
    })
    grid.point$group_label<-grid.point$grp_x+(grid.point$grp_y-1)*k
    group_sample<-aggregate(grid.point$Observed,by=list(grid.point$group_label),FUN=function(x){sample(x,1)})
    group_sum<-aggregate(grid.point$Observed,by=list(grid.point$group_label),FUN=sum)
    group_count<-aggregate(rep(1,dim(grid.point)[1]),by=list(grid.point$group_label),FUN=sum)
    group<-cbind(group_sample,(group_sum$x-group_sample$x)/(group_count$x-1))
    names(group)<-c("group","sample","prob")
    rm(group_sample,group_sum,group_count)
    location<-locate.change(group$sample[order(-group$prob)])[1]
    signif_group<-order(-group$prob)[1:location]
    detection[grid.point$group_label %in% signif_group]=1+detection[grid.point$group_label %in% signif_group]
  }
  detection<-detection/times
  detection[detection<1]=0
  return(detection)
}




multi_iter<-function(N.grid=100,shift=0,scale=0,times=100,r=1,alpha=0.05){
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
  
    grid.point$detection<-signal_detection(grid.point,times=100)
    
    res<-grid.point$Status-grid.point$detection
    wrong_report<-length(res[res==-1])
    miss_report<-length(res[res==1])

  
  #fdrl_method
  pvalues<-pvalue_normal(grid.point$Observed)
  fdrl<-FDRLMethod(pvalues,grid.point$x,grid.point$y,window=r,alpha=alpha)
  fdrl_Res<-grid.point$Status-fdrl$ind
  fdrl_wrong<-length(fdrl_Res[fdrl_Res==-1])
  fdrl_miss<-length(fdrl_Res[fdrl_Res==1])
  
  ret<-c(miss_report/sum(grid.point$Status),wrong_report/sum(1-grid.point$Status),wrong_report/sum(grid.point$detection),
         fdrl_miss/sum(grid.point$Status),fdrl_wrong/sum(1-grid.point$Status),fdrl_wrong/fdrl$numAlt)
  return(ret)
  
}

# ggplot(grid.point)+geom_point(aes(x=x,y=y,colour=detection))+scale_color_continuous(low="white",high="black",limits=c(0,1))
library(doParallel)
cl=makeCluster(13)
registerDoParallel(cl)
aa_05=foreach(m_time=1:100,.combine=rbind,
             .packages=c("RandomFields","foreach"),
             .errorhandling="remove")%dopar%
             {
              multi_iter(N.grid=100,shift=05,scale=0,times=100,r=1.1,alpha=0.01)
             }
stopCluster(cl)
saveRDS(aa_05,"/vol/data/zhuz/xinzhang/spatial_cusum/aa_05.rds")

cl=makeCluster(13)
registerDoParallel(cl)
aa_08=foreach(m_time=1:100,.combine=rbind,
             .packages=c("RandomFields","foreach"),
             .errorhandling="remove")%dopar%
             {
               multi_iter(N.grid=100,shift=0.8,scale=0,times=100,r=1.1,alpha=0.01)
             }
stopCluster(cl)
saveRDS(aa_08,"/vol/data/zhuz/xinzhang/spatial_cusum/aa_08.rds")

cl=makeCluster(13)
registerDoParallel(cl)
aa_1=foreach(m_time=1:100,.combine=rbind,
             .packages=c("RandomFields","foreach"),
             .errorhandling="remove")%dopar%
             {
               multi_iter(N.grid=100,shift=1,scale=0,times=100,r=1.1,alpha=0.01)
             }
stopCluster(cl)
saveRDS(aa_1,"/vol/data/zhuz/xinzhang/spatial_cusum/aa_1.rds")

cl=makeCluster(13)
registerDoParallel(cl)
aa_15=foreach(m_time=1:100,.combine=rbind,
             .packages=c("RandomFields","foreach"),
             .errorhandling="remove")%dopar%
             {
               multi_iter(N.grid=100,shift=1.5,scale=0,times=100,r=1.1,alpha=0.01)
             }
stopCluster(cl)
saveRDS(aa_15,"/vol/data/zhuz/xinzhang/spatial_cusum/aa_15.rds")

cl=makeCluster(13)
registerDoParallel(cl)
aa_2=foreach(m_time=1:100,.combine=rbind,
             .packages=c("RandomFields","foreach"),
             .errorhandling="remove")%dopar%
             {
               multi_iter(N.grid=100,shift=1,scale=0,times=100,r=1.1,alpha=0.01)
             }
stopCluster(cl)
saveRDS(aa_2,"/vol/data/zhuz/xinzhang/spatial_cusum/aa_2.rds")

