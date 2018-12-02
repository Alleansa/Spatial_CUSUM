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

get_k_neigh_med<-function(e,k){
  N<-length(e)
  T<-median(e[1:k])
  return(T)
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
  
  
  #normalized_Observed<-(Observed-0)/1
  pv<-2*(1-pnorm(abs(normalized_Observed)))
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

FDRMethod<-function(pvalues,alpha){
  fdr_p<-p.adjust(pvalues, method = "fdr",n=length(pvalues))
  return(fdr_p<alpha)
}
match_alpha<-function(pvalues,status,x,y,miss_report,window=5,n=51){
  fdrl_miss_report<-rep(0,n)
  fdrl_wrong_report<-rep(0,n)
  fdr_miss_report<-rep(0,n)
  fdr_wrong_report<-rep(0,n)
  for (i in 1L:n) {
    fdrlres<-FDRLMethod(pvalues,x,y,window,alpha=1/(n-1)*(i-1))
    fdrlche_res<-status-fdrlres$ind
    fdrl_miss_report[i]<-length(fdrlche_res[fdrlche_res==1])
    fdrl_wrong_report[i]<-length(fdrlche_res[fdrlche_res==-1])
    fdrres<-FDRMethod(pvalues,alpha=1/(n-1)*(i-1))
    fdrche_res<-status-fdrres
    fdr_miss_report[i]<-length(fdrche_res[fdrche_res==1])
    fdr_wrong_report[i]<-length(fdrche_res[fdrche_res==-1])
  }
  fdrldiff<-fdrl_miss_report-miss_report
  fdrl_location<-which.min(abs(fdrldiff))
  fdrdiff<-fdr_miss_report-miss_report
  fdr_location<-which.min(abs(fdrdiff))
  return(c(fdrl_location,fdr_location))
}

one_iter<-function(N.grid,r,k,Dist_matrix,scale=0, shift =0, threshold=0.4627){
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
  pvalues<-pvalue_normal(grid.point$Observed)
  
  n<-length(grid.point$Observed)
  
  neigh<-neigh_pick(grid.point$Observed,Dist_matrix,r,k+1)
  neigh_mean <-matrix(rep(1/k,k),1,k)%*%neigh[2:(k+1),]
  neigh_sample<-neigh[1,]
  #test whether there is anomaly
  test_seq<-neigh_sample[order(-neigh_mean)]

  if(mean(CUSUM(test_seq/mad(test_seq)))>threshold) {
    location<-locate.change(test_seq)[1]
    grid.point$detection<-rep(0,n)
    grid.point$detection[order(-neigh_mean)[1:location]]<-1
    res<-grid.point$Status-grid.point$detection
    wrong_report<-length(res[res==-1])
    miss_report<-length(res[res==1])
  }else{
    miss_report<-sum(grid.point$Status)
    wrong_report<-0
  }
  aprox_fdr<-match_alpha(pvalues,grid.point$Status,grid.point$x,grid.point$y,miss_report,window =r)
  fdrres_1<-FDRLMethod(pvalues,grid.point$x,grid.point$y,window=r,alpha=(aprox_fdr[1]-1)/50)
  fdrres_2<-FDRMethod(pvalues,alpha=(aprox_fdr[2]-1)/50)
  return(c(grid.point$detection,
           fdrres_1$ind,
           fdrres_2
  ))
}

N.grid<-100
D<-dist_matrix(N.grid)
library(doParallel)
cl=makeCluster(13)
registerDoParallel(cl)
# aa=foreach(i=1:500,.combine=rbind,.packages=c("rootSolve","MASS","RandomFields","geoR"),.export="pp",.errorhandling="remove")%dopar%
aa1=foreach(m_time=1:100,.combine=rbind,
            .packages="RandomFields",
            .errorhandling="remove")%dopar%
            {
              one_iter(N.grid,r=0.05,k=10,Dist_matrix=D,scale=0, shift =1 )
            }
stopCluster(cl)
saveRDS(colMeans(aa1), "/home/xinzhang/spatial_CUSUM/diagram/aa1.rds")

cl=makeCluster(13)
registerDoParallel(cl)
# aa=foreach(i=1:500,.combine=rbind,.packages=c("rootSolve","MASS","RandomFields","geoR"),.export="pp",.errorhandling="remove")%dopar%
aa2=foreach(m_time=1:100,.combine=rbind,
            .packages="RandomFields",
            .errorhandling="remove")%dopar%
            {
              one_iter(N.grid,r=0.05,k=10,Dist_matrix=D,scale=0, shift =2 )
            }
stopCluster(cl)
saveRDS(colMeans(aa2), "/home/xinzhang/spatial_CUSUM/diagram/aa2.rds")

cl=makeCluster(13)
registerDoParallel(cl)
# aa=foreach(i=1:500,.combine=rbind,.packages=c("rootSolve","MASS","RandomFields","geoR"),.export="pp",.errorhandling="remove")%dopar%
aa3=foreach(m_time=1:100,.combine=rbind,
            .packages="RandomFields",
            .errorhandling="remove")%dopar%
            {
              one_iter(N.grid,r=0.05,k=10,Dist_matrix=D,scale=0, shift =3 )
            }
stopCluster(cl)
saveRDS(colMeans(aa3), "/home/xinzhang/spatial_CUSUM/diagram/aa3.rds")

library(ggplot2)
aa1<-readRDS("~/desktop/H1/diagram/aa1.rds")
m1<-(aa1>0.9)
N.grid<-100
grid.point = expand.grid(x=seq(0,1,by=1/N.grid),y=seq(0,1,by=1/N.grid))
grid.point$CUSUM<-m1[1:(N.grid+1)^2]
aa<-ggplot(grid.point)+geom_point(aes(x=x,y=y,colour=CUSUM))+theme(legend.position='none')+scale_color_manual(breaks = c("0", "1"),values=c("white", "black"))
grid.point$FDRl<-m1[((N.grid+1)^2+1):(2*(N.grid+1)^2)]
bb<-ggplot(grid.point)+geom_point(aes(x=x,y=y,colour=FDRl))+theme(legend.position='none')+scale_color_manual(breaks = c("0", "1"),values=c("white", "black"))
grid.point$FDR<-m1[(2*(N.grid+1)^2+1):(3*(N.grid+1)^2)]
cc<-ggplot(grid.point)+geom_point(aes(x=x,y=y,colour=FDR))+theme(legend.position='none')+scale_color_manual(breaks = c("0", "1"),values=c("white", "black"))

aa2<-readRDS("~/desktop/H1/diagram/aa2.rds")
m2<-(aa2>0.85)
N.grid<-100
grid.point = expand.grid(x=seq(0,1,by=1/N.grid),y=seq(0,1,by=1/N.grid))
grid.point$CUSUM<-m2[1:(N.grid+1)^2]
dd<-ggplot(grid.point)+geom_point(aes(x=x,y=y,colour=CUSUM))+theme(legend.position='none')+scale_color_manual(breaks = c("0", "1"),values=c("white", "black"))
grid.point$FDRl<-m2[((N.grid+1)^2+1):(2*(N.grid+1)^2)]
ee<-ggplot(grid.point)+geom_point(aes(x=x,y=y,colour=FDRl))+theme(legend.position='none')+scale_color_manual(breaks = c("0", "1"),values=c("white", "black"))
grid.point$FDR<-m2[(2*(N.grid+1)^2+1):(3*(N.grid+1)^2)]
ff<-ggplot(grid.point)+geom_point(aes(x=x,y=y,colour=FDR))+theme(legend.position='none')+scale_color_manual(breaks = c("0", "1"),values=c("white", "black"))

aa3<-readRDS("~/desktop/H1/diagram/aa3.rds")
m3<-(aa3>0.6)
N.grid<-100
grid.point = expand.grid(x=seq(0,1,by=1/N.grid),y=seq(0,1,by=1/N.grid))
grid.point$CUSUM<-m3[1:(N.grid+1)^2]
gg<-ggplot(grid.point)+geom_point(aes(x=x,y=y,colour=CUSUM))+theme(legend.position='none')+scale_color_manual(breaks = c("0", "1"),values=c("white", "black"))
grid.point$FDRl<-m3[((N.grid+1)^2+1):(2*(N.grid+1)^2)]
hh<-ggplot(grid.point)+geom_point(aes(x=x,y=y,colour=FDRl))+theme(legend.position='none')+scale_color_manual(breaks = c("0", "1"),values=c("white", "black"))
grid.point$FDR<-m3[(2*(N.grid+1)^2+1):(3*(N.grid+1)^2)]
ii<-ggplot(grid.point)+geom_point(aes(x=x,y=y,colour=FDR))+theme(legend.position='none')+scale_color_manual(breaks = c("0", "1"),values=c("white", "black"))

library(grid)
grid.newpage()  ##newpage
pushViewport(viewport(layout = grid.layout(2,2))) ### divide into 1*3
vplayout <- function(x,y){
  viewport(layout.pos.row = x, layout.pos.col = y)
}
print(aa, vp = vplayout(1,1))   ### put a in (1,1)
print(bb, vp = vplayout(1,2))   ### put b in (1,2)
print(cc, vp = vplayout(2,1))  ###put c in (1,3)
print(dd, vp = vplayout(2,2))  ###put c in (1,3)
print(ee, vp = vplayout(2,2))   ### put a in (1,1)
print(ff, vp = vplayout(2,3))   ### put b in (1,2)
print(gg, vp = vplayout(3,1))  ###put c in (1,3)
print(hh, vp = vplayout(3,2))  ###put c in (1,3)
print(ii, vp = vplayout(3,3))  ###put c in (1,3)