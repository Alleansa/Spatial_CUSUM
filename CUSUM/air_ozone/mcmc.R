
########################################################:
#################     Inputs             ###############:
########################################################:
#
# y          = n-vector of responses
# s          = n x 2 matrix of spatial locations
# X          = n x p matrix of spatial covariates
# sp         = np x 2 matrix of prediction locations
# Xp         = np x p matrix of spatial covariates
# cutoff     = Ho: signal < cutoff; Ha: signal > cutoff
# mean_nu    = prior mean of the log matern smoothness parameter
# sd_nu      = prior sd of the log matern smoothness parameter
# init_nu    = initial value of the log matern smoothness parameter
# mean_range = prior mean of the log matern spatial range parameter
# sd_range   = prior sd of the log matern spatial range parameter
# init_range = initial value of the log matern spatial range parameter
# a_sill     = the prior for the varaince is InvG(a_sill,b_sill)
# init_sill  = initial value of the sill
# sd_beta    = the regression coefs have N(0,sd_beta) priros
# init_beta  = initial value of the regression coefs
# nugget     = include a nugget?
# oracle     = fix the regression coefs and spatial covariance parameters
#              at the initial values?
# iters      = number of MCMC samples to generate
# burn       = number to discard as burnin
# update     = number of iterations between graphical displays
# thin       = degree of thinning (thin*iters update are made, iters are returned)
#
########################################################:


########################################################:
#################     Outputs            ###############:
########################################################:
#
# pred       = (iters-burn) x np matrix of predicted values at sp
# beta       = iters x p matrix of posterior samples of the regression coefs
# keepers    = posterior samples of spatial covariance parameters
#
########################################################:

SpatialMCMC<-function(y,s,X,
                      sp=NULL,Xp=NULL,cutoff=0,
                      mean_nu=-1,sd_nu=1,init_nu=NULL,
                      mean_range=-1,sd_range=1,init_range=NULL,
                      init_r=NULL,
                      a_sill=.01,b_sill=.01,init_sill=NULL,
                      sd_beta=1000,init_beta=NULL,
                      nugget=TRUE,oracle=FALSE,
                      iters=500,burn=100,update=100,thin=2){
  
  library(emulator)
  
  #Bookkeeping
  n<-length(y)
  p<-dim(X)[2]
  d<-rdist(s,s)
  diag(d)<-0
  predictions<-!is.null(sp) & !is.null(Xp)
  np<-1
  if(predictions){
    np<-nrow(sp)
    d12<-rdist(sp,s)
    d11<-rdist(sp,sp)
    diag(d11)<-0
  }
  tX<-t(X)
  precb<-diag(p)/(sd_beta^2)
  
  #Initial values
  mu<-y*0.9+mean(y)*0.1
  beta<-init_beta
  r<-init_r
  tauinv<-init_sill
  rhos<-rhoe<-init_range
  nus<-nue<-init_nu
  
  lmfit<-lm(y~X-1)
  if(is.null(beta)){beta<-lmfit$coef}
  if(is.null(tauinv)){tauinv<-var(lmfit$res)}
  if(is.null(rhos)){rhos<-rhoe<-quantile(d,0.1)}
  if(is.null(r)){r<-0.8}
  if(is.null(nus)){nus<-nue<-0.5}
  tau<-1/tauinv
  rm(lmfit)
  r<-ifelse(nugget,r,1)
  Es<-inv(d,1,rhos,nus)
  Ee<-inv(d,r,rhoe,nue)
  Xb   <- as.vector(X%*%beta)
  
  if(oracle){
    #Prep to update mu
    VVV_MU <- solve(tau*Es$PREC + Ee$PREC)
    MMM_MU <- tau*(Es$PREC%*%Xb) + Ee$PREC%*%y
    MMM_MU <-VVV_MU%*%MMM_MU
    PPP_MU <- t(chol(VVV_MU))
    
    #Prep to interpolate mu
    S11<-corfx(d11,1,rhos,nus)/tau
    S12<-corfx(d12,1,rhos,nus)/tau
    S22inv<-tau*Es$PREC
    np<-nrow(S12)
    MMM_PROJ<-S12%*%S22inv
    PPP_PROJ<-S11-MMM_PROJ%*%t(S12)
    PPP_PROJ<-t(chol(PPP_PROJ))
  }
  
  yp<-matrix(0,iters,np)
  keep.beta<-matrix(0,iters,p)
  keepers<-matrix(0,iters,6)
  colnames(keepers)<-c("sigma","r","range_s","nu_s","range_e","nu_e")
  
  CCC<-diag(2)*1.5-0.5
  PPPnu<-0.1*t(chol(CCC))
  
  
  for(i in 1:iters){
    for(thinnumber in 1:thin){
      
      ##############################################:
      #####   SPATIAL SIGNAL/RANDOM EFFECTS  #######:
      ##############################################:
      
      if(!oracle){
        VVV_MU <- solve(tau*Es$PREC + Ee$PREC)
        MMM_MU <- tau*(Es$PREC%*%Xb) + Ee$PREC%*%y
        MMM_MU <- VVV_MU%*%MMM_MU
        PPP_MU <- t(chol(VVV_MU))
      }
      mu  <- MMM_MU+PPP_MU%*%rnorm(n)
      mu  <- as.vector(mu)
      
      
      if(!oracle){
        
        ##############################################:
        #####        REGRESSION COEFS          #######:
        ##############################################:
        
        XXX  <- tau*tX%*%Es$PREC
        VVV  <- solve(precb+XXX%*%X)
        MMM  <- rep(0,p)+XXX%*%mu
        beta <- VVV%*%MMM+t(chol(VVV))%*%rnorm(p)
        Xb   <- as.vector(X%*%beta)
        
        
        ##############################################:
        #####     SIGNAL COVAR PARAMETERS      #######:
        ##############################################:
        
        #Sill
        SS  <- quad.form(Es$PREC,mu-Xb)
        tau <- rgamma(1,n/2+a_sill,SS/2+b_sill)
        SS  <- tau*SS
        
        
        #spatial range/smoothness:
        eee<-PPPnu%*%rnorm(2)
        canrhos <- exp(log(rhos)+eee[1])
        cannus <- exp(log(nus)+eee[2])
        canEs<-inv(d,1,canrhos,cannus)
        canSS<-tau*quad.form(canEs$PREC,mu-Xb)
        R<-dnorm(log(cannus),mean_nu,sd_nu,log=T)-
          dnorm(log(nus),mean_nu,sd_nu,log=T)+
          dnorm(log(canrhos),mean_range,sd_range,log=T)-
          dnorm(log(rhos),mean_range,sd_range,log=T)+
          0.5*(canEs$DET-Es$DET)-
          0.5*(canSS-SS)
        if(!is.na(exp(R))){if(runif(1)<exp(R)){
          nus<-cannus;rhos<-canrhos;Es<-canEs
        }}
        
        
        
        ##############################################:
        #####      ERROR COVAR PARAMETERS      #######:
        ##############################################:
        
        #Signal to noise ratio
        SS<-quad.form(Ee$PREC,y-mu)
        if(nugget){
          z     <- log(r/(1-r))
          canz  <- rnorm(1,z,0.5)
          canr  <- exp(canz)/(1+exp(canz))
          canEe <- inv(d,canr,rhoe,nue)
          canSS <- quad.form(canEe$PREC,y-mu)
          R<-dnorm(canz,log=T)-
            dnorm(z,log=T)+
            0.5*(canEe$DET-Ee$DET)-
            0.5*(canSS-SS)
          if(!is.na(exp(R))){if(runif(1)<exp(R)){
            r<-canr;Ee<-canEe;SS<-canSS
          }}
        }
        
        #spatial range/smoothness:
        eee<-PPPnu%*%rnorm(2)
        canrhoe <- exp(log(rhoe)+eee[1])
        cannue <- exp(log(nue)+eee[2])
        canEe<-inv(d,r,canrhoe,cannue)
        canSS<-quad.form(canEe$PREC,y-mu)
        R<-dnorm(log(cannue),mean_nu,sd_nu,log=T)-
          dnorm(log(nue),mean_nu,sd_nu,log=T)+
          dnorm(log(canrhoe),mean_range,sd_range,log=T)-
          dnorm(log(rhoe),mean_range,sd_range,log=T)+
          0.5*(canEe$DET-Ee$DET)-
          0.5*(canSS-SS)
        if(!is.na(exp(R))){if(runif(1)<exp(R)){
          nue<-cannue;rhoe<-canrhoe;Ee<-canEe
        }} 
      }#end if oracle
      
    }#end thin
    
    keep.beta[i,]<-beta
    keepers[i,]<-c(1/sqrt(tau),r,rhos,nus,rhoe,nue)
    
    # plots the samples
    if(i%%update==0){
      these<-1:i
      these<-these[these>i-500]
      par(mfrow=c(2,2))
      plot(these,keepers[these,1],type="l",ylab="Sill (sd)")
      plot(these,keepers[these,2],type="l",ylab="S2N ratio")
      plot(these,keepers[these,3],type="l",ylab="Range")
      plot(these,keepers[these,4],type="l",ylab="Smoothness")
    }
    
    if(i>=burn){if(predictions & !oracle){
      S11<-corfx(d11,1,rhos,nus)/tau
      S12<-corfx(d12,1,rhos,nus)/tau
      S22inv<-tau*Es$PREC
      yp[i,]<-Xp%*%beta+proj(mu-Xb,S12,S11,S22inv)       
    }}
    
    if(i>=burn){if(predictions & oracle){
      yp[i,]<-Xp%*%beta+MMM_PROJ%*%(mu-Xb)+PPP_PROJ%*%rnorm(np)       
    }}
    
  }   
  
  
  list(pred=yp[burn:iters,],beta=keep.beta,keepers=keepers)}




corfx<-function(d,r,rho,nu){
  library(geoR)
  r*matern(d,rho,nu)+(1-r)*ifelse(d==0,1,0)
}

inv<-function(d,r,range=NULL,nu=0.5,thresh=0.0000001){
  library(geoR)
  Q<-corfx(d,r,range,nu)
  eig<-eigen(Q)
  V<-Re(eig$vectors)
  D<-Re(eig$values)
  D<-1/ifelse(D<thresh,thresh,D)
  COR.INV<-sweep(V,2,D,"*")%*%t(V)
  list(PREC=COR.INV,DET=sum(log(D)))}


proj<-function(y,S12,S11,S22inv){
  np   <- nrow(S12)
  XXX  <- S12%*%S22inv
  mn   <- XXX%*%y
  E    <- eigen(S11-XXX%*%t(S12))
  PPP  <- Re(E$vectors)
  D    <- Re(E$values)
  D    <- sqrt(ifelse(D>0,D,0))
  yp   <- mn+PPP%*%rnorm(np,0,D)
  return(yp)}



FDR<-function(theta,alpha=.1,nthresh=100){
  #pick reject so that E(mean(theta[reject]))<alpha
  
  inprob<-apply(theta,2,mean)
  thresh<-seq(0,1,length=nthresh)
  BFDR<-rep(0,nthresh)
  for(j in 1:nthresh){if(sum(inprob>=thresh[j])>0){
    BFDR[j]<-1-mean(inprob[inprob>=thresh[j]])
  }}
  level<-min(thresh[BFDR<alpha])
  reject<-inprob>level
  list(level=level,reject=reject,thresh=thresh,BFDR=BFDR)}


FDX<-function(theta,alpha=.1,beta=.9,nthresh=100){
  #pick reject so that P(mean(theta[reject])<alpha)<beta
  
  reject<-rep(0,nrow(theta))
  level<-0
  
  inprob<-apply(theta,2,mean)
  thresh<-seq(0,max(inprob),length=nthresh)
  BFDX<-rep(0,nthresh)
  for(j in 1:nthresh){
    these<- inprob>=thresh[j]
    if(sum(these)>1){
      BFDX[j]<-mean(apply(1-theta[,these],1,mean)<alpha)
    }
  }
  level<-1
  if(sum(BFDX>beta)>0){
    level<-min(thresh[BFDX>beta])
  }
  reject<-inprob>level
  list(level=level,reject=reject,thresh=thresh,BFDX=BFDX)}

get.theta.clust<-function(sig,mu0,clust){
  nclust<-max(clust)
  theta<-rep(0,nclust)
  for(j in 1:nclust){
    theta[j]<-ifelse(mean(sig[clust==j])>0.2,1,0)
  }
  return(theta==1)}

sample.s<-function(n,tooclose=0.01,maxiters=25){
  s<-cbind(runif(n),runif(n))
  d<-as.matrix(dist(s))
  diag(d)<-tooclose+1
  md<-apply(d,1,min)
  count<-0
  while(min(md)<tooclose & count<maxiters){
    count<-count+1
    s[which(md<tooclose)[1],]<-runif(2)
    d<-as.matrix(dist(s))
    diag(d)<-tooclose+1
    md<-apply(d,1,min)
  }
  return(s)}
