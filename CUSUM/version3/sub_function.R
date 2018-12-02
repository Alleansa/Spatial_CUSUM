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
  pv<-(1-pnorm(normalized_Observed))
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


spatial_CUSUM<-function(grid.point,r,threshold=0.4627){
  n<-length(grid.point$Observed)
  Dist_matrix<-as.matrix(dist(cbind(grid.point[,1:2]),upper=1,diag=1))
  neigh<-neigh_pick(grid.point$Observed,Dist_matrix,r,k)
  neigh_sample<-matrix(rep(1/k,k),1,k)%*%neigh
  test_seq<-neigh_sample[order(-grid.point$Observed)]
  if(mean(CUSUM(test_seq/mad(test_seq)))<=threshold) {
    print('There is no anomaly in this region')
    grid.point$detection<-rep(0,n)
  }else {
    location<-locate.change(test_seq)[1]
    grid.point$detection<-rep(0,n)
    grid.point$detection[order(-grid.point$Observed)[1:location]]<-1
  }
  return(grid.point)
}