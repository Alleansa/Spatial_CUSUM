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


