#version 4
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

#FDR on boundary
boundary_fdr<-function(signif_group,grid.point,alpha=0.05){
  k<-max(grid.point$group_label)
  group<-expand.grid(x=seq(1,sqrt(k),by=1),y=seq(1,sqrt(k),by=1))
  group$label<-numeric(sqrt(k))
  group$label[signif_group]=1
  for (i in c(1:sqrt(k))){
    for (j in c(1:sqrt(k))){
      group$boundary[group$x==i&group$y==j]=mean(c(group$label[group$x==i&group$y==(j+1)],
                                                   group$label[group$x==i&group$y==(j-1)],
                                                   group$label[group$x==(i+1)&group$y==j],
                                                   group$label[group$x==(i-1)&group$y==j],
                                                   group$label[group$x==i&group$y==j])
                                                   )
    }
  }
  group$boundary[which(group$boundary==1)]=0
  group$boundary=ceiling(group$boundary)
  
  grid.point$boundary<-numeric(dim(grid.point)[1])
  grid.point$boundary[grid.point$group_label %in% which(group$boundary==1)]=1
  mu_normal<-median(grid.point$Observed[grid.point$detection==0])
  var_normal<-mad(grid.point$Observed[grid.point$detection==0])
  pvalue<-1-pnorm((grid.point$Observed[grid.point$boundary==1]-mu_normal)/sqrt(var_normal))
  fdr_p<-p.adjust(pvalue, method = "fdr",n=length(pvalue))
  new_detection<-numeric(dim(grid.point)[1])
  new_detection[grid.point$detection==1&grid.point$boundary==0]=1
  new_detection[which(grid.point$boundary==1)[which(fdr_p<alpha)]]=1
  return(new_detection)
}

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
  normalized_Observed<-(Observed-median(Observed))/sqrt(mad(Observed))
  pv<-1-pnorm(normalized_Observed)
  return(pv)
}

FDRLMethod <- function(pvalues,x ,y, window, alpha, nstep = 200, lambda = 0.1){
  if( !is(window,"integer") ) window <- as.integer(round(window,0))
  
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

match_wrong<-function(pvalues,status,x,y,wrong_report,window=5,n=51){
  fdr_miss_report<-numeric(n)
  fdr_wrong_report<-numeric(n)
  for (i in 1L:n) {
    res<-FDRLMethod(pvalues,x,y,window,alpha=1/(n-1)*(i-1))
    che_res<-status-res$ind
    fdr_miss_report[i]<-length(che_res[che_res==1])
    fdr_wrong_report[i]<-length(che_res[che_res==-1])
  }
  diff<-wrong_report-fdr_wrong_report
  location<-which.min(abs(diff))
  return(c(fdr_wrong_report[location],fdr_miss_report[location],location))
}

#main
multi_iter<-function(N.grid=100,k=20,shift=0,scale=0,times=100){
  Status<-rep(0,(N.grid)^2)
  grid.point = expand.grid(x=seq(1/N.grid,1,by=1/N.grid),y=seq(1/N.grid,1,by=1/N.grid))
  
  
  if(scale==0){
    #nondependence
    grid.point$Observed<-rnorm((N.grid)^2)
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
               one_iter(grid.point,k)
             }
  if(times==1){
    aa<-t(aa)
  }
  
  ret<-round(c(mean(aa[,1])/sum(grid.point$Status),
               sum(aa[,2])/times,
               sum(aa[,3])/((times-sum(aa[,2]))*((N.grid)^2-sum(grid.point$Status))),
               mean(aa[,4])/sum(grid.point$Status),
               sum(aa[,4]==sum(grid.point$Status)&aa[,5]==0)/times,
               sum(aa[,5])/((times-sum(aa[,4]==sum(grid.point$Status)&aa[,5]==0))*((N.grid)^2-sum(grid.point$Status)))),3)
  
  return(ret)
}


one_iter<-function(grid.point,k,threshold=0.4627){
  fail_count<-0
  miss_report<-0
  wrong_report<-0
  fdr_miss_report<-0
  fdr_wrong_report<-0
  
  pvalues<-pvalue_normal(grid.point$Observed)

  grid.point <- within(grid.point, {
    grp_x <- cut(x, (0:k)/k, labels = FALSE)
    grp_y <- cut(y, (0:k)/k, labels = FALSE)
  })
  grid.point$group_label<-grid.point$grp_x+(grid.point$grp_y-1)*k
  group_sample<-aggregate(grid.point$Observed,by=list(grid.point$group_label),FUN=function(x){sample(x,1)})
  group_sum<-aggregate(grid.point$Observed,by=list(grid.point$group_label),FUN=sum)
  group_count<-aggregate(rep(1,dim(grid.point)[1]),by=list(grid.point$group_label),FUN=sum)
  group<-cbind(group_sample,(group_sum$x-group_sample$x)/group_count$x)
  names(group)<-c("group","sample","prob")
  rm(group_sample,group_sum,group_count)
  if(mean(CUSUM(group$sample[order(-group$prob)]/mad(group$sample)))<=threshold) {
    fail_count<-1
    
    wrong_report<-0
    miss_report<-sum(grid.point$Status)
    aprox_fdr<-match_wrong(pvalues,grid.point$Status,grid.point$x,grid.point$y,wrong_report,window=1/k)
    fdr_wrong_report<-aprox_fdr[1]
    fdr_miss_report<-aprox_fdr[2]
    
  }else {
    location<-locate.change(group$sample[order(-group$prob)])[1]
    signif_group<-order(-group$prob)[1:location]
    grid.point$detection<-numeric(dim(grid.point)[1])
    grid.point$detection[grid.point$group_label %in% signif_group]=1
    grid.point$detection<-boundary_fdr(signif_group,grid.point)

    res<-grid.point$Status-grid.point$detection
    wrong_report<-length(res[res==-1])
    miss_report<-length(res[res==1])
    aprox_fdr<-match_wrong(pvalues,grid.point$Status,grid.point$x,grid.point$y,wrong_report,window =1/k)
    fdr_wrong_report<-aprox_fdr[1]
    fdr_miss_report<-aprox_fdr[2]
    
  }
  
  return(c(miss_report,
           fail_count,
           wrong_report,
           fdr_miss_report,
           fdr_wrong_report))
}



#N.grid=30
library(doParallel)
cl=makeCluster(13)
registerDoParallel(cl)
aa_1=foreach(m_time=1:100,.combine=rbind,
             .packages=c("RandomFields","foreach"))%dopar%
             {
               multi_iter(N.grid=100,k=20,shift=1,scale=0,times=1)
             }
stopCluster(cl)
saveRDS(aa_1, "/home/xinzhang/spatial_CUSUM/version1/aa_1.rds")

cl=makeCluster(13)
registerDoParallel(cl)
aa_2=foreach(m_time=1:100,.combine=rbind,
             .packages=c("RandomFields","foreach"))%dopar%
             {
               multi_iter(N.grid=100,k=20,shift=2,scale=0,times=1)
             }
stopCluster(cl)
saveRDS(aa_2, "/home/xinzhang/spatial_CUSUM/version1/aa_2.rds")

cl=makeCluster(13)
registerDoParallel(cl)
aa_3=foreach(m_time=1:100,.combine=rbind,
             .packages=c("RandomFields","foreach"))%dopar%
             {
               multi_iter(N.grid=100,k=20,shift=3,scale=0,times=1)
             }
stopCluster(cl)
saveRDS(aa_3,"/home/xinzhang/spatial_CUSUM/version1/aa_3.rds")

cl=makeCluster(13)
registerDoParallel(cl)
aa_5=foreach(m_time=1:100,.combine=rbind,
             .packages=c("RandomFields","foreach"),
             .errorhandling="remove")%dopar%
             {
               multi_iter(N.grid=100,k=20,shift=5,scale=0,times=1)
             }
stopCluster(cl)
saveRDS(aa_5,"/home/xinzhang/spatial_CUSUM/version1/aa_5.rds")
