load("/work/STAT/xinzhang/spatial_mcmc/data.rda")
source("/work/STAT/xinzhang/spatial_mcmc/mcmc.R")
library(fields)
library(geoR)
iters  <- 6000                # MCMC iterations
burn   <- 1000                # MCMC burn-in
thin   <- 1                   # Degree of thinning (it generates thin*iters samples)
update <- 10                  # Iterations between graphical displays
set.seed(rep*0820)
Pred.N.grid<-100
s2<-expand.grid(seq(0,1,length=Pred.N.grid),seq(0,1,length=Pred.N.grid))
s1<-data0301[,17:18]
X1<-matrix(1,176,1)
X2<-matrix(1,Pred.N.grid^2,1)
fit<-SpatialMCMC(y=y1,s=s1,X=X1,sp=s2,Xp=X2,burn=burn,iters=iters,update=update,thin=thin)
saveRDS(fit, "/work/STAT/xinzhang/spatial_mcmc/Results/fit.rds")