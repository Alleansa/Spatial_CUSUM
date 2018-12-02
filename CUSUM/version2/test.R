rm(list=ls())
#******************************************************************************#
# FDRL : public function implementing FDR_L method of Zhang et al.             #
#******************************************************************************#
#                                                                              #
# Input                                                                        #
#                                                                              #
# pvalues  : an object of class numeric.                                       #
#            a vector of pvalues                                               #
#                                                                              #
# window   : an object of class numeric                                        #
#            size of the window defining the neighborhood                      #
#                                                                              #
# alpha    : an object of class numeric                                        #
#            the level of significance for determining the critical value      #
#                                                                              #
# nstep    : an object of class numeric                                        #
#            the number of threshold values to consider                        #
#                                                                              #
# lambda   : an object of class numeric                                        #
#            tuning constant                                                   #
#                                                                              #
# Output                                                                       #
#                                                                              #
#  A list is returned                                                          #
#                                                                              #
#    ind : a vector of indicator functions. (1) rejected (0) not-rejected      #
#                                                                              #
#    threshold : critical value                                                #
#                                                                              #
#    numAlt : number of rejected hypotheses                                    #
#                                                                              #
#    propAlt : proportion of hypotheses rejected.                              #
#                                                                              #
#******************************************************************************#



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
  Dist_matrix <- get_distance(x,y)
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

match_wrong<-function(pvalues,status,x,y,miss_report,window=5,n=51){
  fdr_miss_report<-rep(0,n)
  fdr_wrong_report<-rep(0,n)
  for (i in 1L:n) {
    res<-FDRLMethod(pvalues,x,y,window,alpha=1/(n-1)*(i-1))
    che_res<-status-res$ind
    fdr_miss_report[i]<-length(che_res[che_res==1])
    fdr_wrong_report[i]<-length(che_res[che_res==-1])
  }
  diff<-fdr_miss_report-miss_report
  location<-which.min(abs(diff))
  return(c(fdr_wrong_report[location],fdr_miss_report[location],location))
}

#FDR
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


#constructing the distance matrix
get_distance <-function(x,y){
  x_matrix=matrix(rep(x,time=length(x)),length(x),length(x))-t(matrix(rep(x,time=length(x)),length(x),length(x)))
  y_matrix=matrix(rep(y,time=length(y)),length(y),length(y))-t(matrix(rep(y,time=length(y)),length(y),length(y)))
  D=sqrt(x_matrix^2+y_matrix^2)
  return(D)
}

#computing the CUSUM for series {e}
CUSUM<-function(e){
  e_star<-c(0,e)
  T<-cumsum(e_star)-mean(e)*c(0:length(e))
  return(T^2/length(e))
}

# neigh info pick
sample_nozero<-function(I,seq,size){
  valid<-t(I*seq)
  return(sample(valid,size,prob=I/sum(I),replace = 1))
}

neigh_pick<-function(e,D,r,k){
  Idct<-D<=r&D>0
  Sample<-apply(Idct,1,sample_nozero,seq=e,size=k)
  return(Sample)
}

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

#main 1
single_region_detection <- function(data,r=0.1,k=1, threshold=0.47){
  n<-length(data$Observed)
  Dist_matrix <- get_distance(data$x,data$y)
  neigh_sample<-neigh_pick(data$Observed,Dist_matrix,r,k)
  neigh_info <-matrix(rep(1/k,k),1,k)%*%neigh_sample
  #test whether there is anomaly
  test_seq<-data$Observed[order(-neigh_info)]
  if(mean(CUSUM(test_seq/mad(test_seq)))<=threshold) {
    print('There is no anomaly in this region')
  }else {
    location<-locate.change(test_seq)[1]
    data$detection<-rep(0,n)
    data$detection[order(-neigh_info)[1:location]]<-1
    data$p<-rep(0,n)
    data$p[data$detection==1]<-data$Observed[data$detection==1]
  }
  return(data)
}

multi_region_detection<-function(data,r=0.1,k=1, threshold=1.90,M=0){
  n<-length(data$Observed)
  Dist_matrix <- get_distance(data$x,data$y)
  neigh_sample<-neigh_pick(data$Observed,Dist_matrix,r,k)
  neigh_info <-matrix(rep(1/k,k),1,k)%*%neigh_sample
  #test whether there is anomaly
  test_seq<-data$Observed[order(-neigh_info)]
  rnd1 <- sample(0:(n - 2), M, replace = TRUE)
  rnd2 <- sample(0:(n - 2), M, replace = TRUE)
  window_s <- pmin(rnd1, rnd2)
  window_e <- pmax(rnd1, rnd2) + 2
  BinSeg <- function(x, s, e, depth, parent.val) {
    if (e - s <= 2) 
      return(NULL)
    ind <- (window_s >= s) & (window_e <= e)
    max.val <- -1
    cp <- 0
    for (m in c(0, ((1:M)[ind]))) {
      if (m == 0) {
        s_m <- s
        e_m <- e
      }else {
        s_m <- window_s[m]
        e_m <- window_e[m]
      }
      obj <- locate.change(x[(s_m + 1):e_m])
      if (obj[2] > max.val) {
        max.val <- obj[2]
        cp <- s_m + obj[1]
      }
    }
    ret <- NULL
    ret$location <- cp
    ret$max.cusum <- min(parent.val, max.val)
    ret$depth <- depth
    if (ret$max.cusum < threshold |depth>log(n)) {
      return(NULL)
    }else {
      return(cbind(BinSeg(x, s, cp, depth+1, ret$max.cusum), 
                   ret, BinSeg(x, cp, e, depth+1, ret$max.cusum)))
    }
  }
  ret <- NULL
  ret$x <- test_seq
  ret$changepoints <- BinSeg(test_seq, 0, n, depth = 1, parent.val = .Machine$double.xmax)
  ret$changepoints <- t(matrix(as.numeric(ret$changepoints), 
                               nrow = 3))
  dete_cp<-ret$changepoints[,1]
  if(length(dete_cp)==0){
    print("there is no anomaly region")
  }else{
    dete_cp<-unique(c(0,dete_cp,n))
    Num_cp<-length(dete_cp)
    data$detection<-rep(0,n)
    for (ii in c(1:(Num_cp-1))){
      data$detection[order(-neigh_info)[(dete_cp[ii]+1):dete_cp[ii+1]]]<-Num_cp-ii-1
    }
  }
  return(data)
}


multi_iter<-function(N.grid=20,r=0.06,k=1,shift=0,scale=0,times=100){
  Status<-rep(0,(N.grid+1)^2)
  grid.point = expand.grid(x=seq(0,1,by=1/N.grid),y=seq(0,1,by=1/N.grid))
  
  
  #nondependence
  grid.point$Observed<-rnorm((N.grid+1)^2)
  
  # #dependence
  # model <- RMspheric(var=1, scale=0.0007)
  # simu.alti <- RFsimulate(model, x=grid.point$x,y=grid.point$y,grid = FALSE)
  # grid.point$Observed = simu.alti@data$variable1
  
  
  #generating 'L' shape
  grid.point$Status<-Status
  grid.point[which(grid.point$x<=0.8&grid.point$x>=.4&grid.point$y<=.8&grid.point$y>=.4),4]<-1
  grid.point[which(grid.point$x<=0.8&grid.point$x>=.6&grid.point$y<=.8&grid.point$y>=.6),4]<-0
  grid.point[grid.point$Status==1,]$Observed<-grid.point[grid.point$Status==1,]$Observed+shift
  
  
  aa=foreach(m_time=1:times,.combine=rbind,
             .packages="RandomFields")%do%
             {
               one_iter(grid.point,r,k)
             }
  if(times==1){
    aa<-t(aa)
  }
  
  ret<-round(c(mean(aa[,1])/sum(grid.point$Status),
               sum(aa[,2])/times,
               sum(aa[,3])/((times-sum(aa[,2]))*((N.grid+1)^2-sum(grid.point$Status))),
               mean(aa[,4])/sum(grid.point$Status),
               sum(aa[,4]==sum(grid.point$Status)&aa[,5]==0)/times,
               sum(aa[,5])/((times-sum(aa[,4]==sum(grid.point$Status)&aa[,5]==0))*((N.grid+1)^2-sum(grid.point$Status)))),3)
  
  return(ret)
}


one_iter<-function(grid.point,r,k,threshold=0.4627){
  fail_count<-0
  miss_report<-0
  wrong_report<-0
  fdr_miss_report<-0
  fdr_wrong_report<-0
  
  pvalues<-pvalue_normal(grid.point$Observed)
  
  # n<-length(grid.point$Observed)
  #get the distance matrix
  n<-length(grid.point$Observed)
  Dist_matrix <- get_distance(grid.point$x,grid.point$y)
  neigh_sample<-neigh_pick(grid.point$Observed,Dist_matrix,r,k)
  neigh_info <-matrix(rep(1/k,k),1,k)%*%neigh_sample
  #test whether there is anomaly
  test_seq<-neigh_info[order(-grid.point$Observed)]
  if(mean(CUSUM(test_seq/mad(test_seq)))<=threshold) {
    #print('There is no anomaly in this region')
    fail_count<-1
    
    miss_report<-sum(grid.point$Status)
    aprox_fdr<-match_wrong(pvalues,grid.point$Status,grid.point$x,grid.point$y,miss_report,window=r)
    fdr_wrong_report<-aprox_fdr[1]
    fdr_miss_report<-aprox_fdr[2]
    
    
    # res<-FDRLMethod(pvalues,grid.point$x,grid.point$y,window=r,alpha=p_thre)
    # fdr_wrong_report<-sum(res$ind)
  }else {
    location<-locate.change(test_seq)[1]
    grid.point$detection<-rep(0,n)
    grid.point$detection[order(-grid.point$Observed)[1:location]]<-1
    # grid.point$boundary_score<-neigh_mean(grid.point$detection,Dist_matrix,r)
    # grid.point$boundary<-rep(0,n)
    # grid.point$boundary[abs(grid.point$boundary_score-0.5)<0.45]<-1
    # grid.point$detection<-boundary_detection(grid.point$detection,grid.point$boundary,grid.point$Observed)
    
    res<-grid.point$Status-grid.point$detection
    wrong_report<-length(res[res==-1])
    miss_report<-length(res[res==1])
    aprox_fdr<-match_wrong(pvalues,grid.point$Status,grid.point$x,grid.point$y,miss_report,window =r)
    fdr_wrong_report<-aprox_fdr[1]
    fdr_miss_report<-aprox_fdr[2]
    
    # wrong_report<-sum(grid.point$detection)
    # fdrres<-FDRLMethod(pvalues,grid.point$x,grid.point$y,window=r,alpha=p_thre)
    # fdr_wrong_report<-sum(fdrres$ind)
  }
  
  return(c(miss_report,
           fail_count,
           wrong_report,
           fdr_miss_report,
           fdr_wrong_report))
}



#N.grid=30
library(doParallel)
cl=makeCluster(12)
registerDoParallel(cl)
aa1=foreach(m_time=1:100,.combine=rbind,
            .packages=c("RandomFields","foreach"))%dopar%
            {
              multi_iter(N.grid=20,r=0.06,1,times =100)
            }

stopCluster(cl)

saveRDS(aa1, "/work/STAT/xinzhang/CUSUM/Results/aa1.rds")





