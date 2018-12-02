rm(list=ls())

load("/Users/Allen/Documents/Study/Ph.D/Code/R/CUSUM/B_inte.rda")

source("/Users/Allen/Documents/Study/Ph.D/Code/R/CUSUM/version2/FDRL.R")
source("/Users/Allen/Documents/Study/Ph.D/Code/R/CUSUM/version2/CUSUM.R")

one_iter<-function(N.grid,B_inte,p_thre,r,k,shift=0){
  Status<-rep(0,(N.grid+1)^2)
  grid.point = expand.grid(x=seq(0,1,by=1/N.grid),y=seq(0,1,by=1/N.grid))
  
  # #nondependence
  # grid.point$Observed<-rnorm((N.grid+1)^2)
  
  #dependence
  model <- RMspheric(var=1, scale=0.05)
  simu.alti <- RFsimulate(model, x=grid.point$x,y=grid.point$y,grid = FALSE)
  grid.point$Observed = simu.alti@data$variable1
  
  
  #generating 'L' shape
  grid.point$Status<-Status
  grid.point[which(grid.point$x<=0.8&grid.point$x>=.4&grid.point$y<=.8&grid.point$y>=.4),4]<-1
  grid.point[which(grid.point$x<=0.8&grid.point$x>=.6&grid.point$y<=.8&grid.point$y>=.6),4]<-0
  grid.point[grid.point$Status==1,]$Observed<-grid.point[grid.point$Status==1,]$Observed+shift
  pvalues<-pvalue_normal(grid.point$Observed)
  
  # n<-length(grid.point$Observed)
  #get the distance matrix
  n<-length(grid.point$Observed)
  Dist_matrix <- get_distance(grid.point$x,grid.point$y)
  neigh_info <-neigh_pick(grid.point$Observed,Dist_matrix,r,k)
  #test whether there is anomaly
  test_seq<-neigh_info[order(-grid.point$Observed)]
  p<-sum((B_inte)>mean(CUSUM(test_seq)))/length(B_inte)
  # p2<-sum(B_max>max(CUSUM(test_seq)))/length(B_max)
  if(p<p_thre) {
    T<-CUSUM(test_seq)
    location<-which.max(T)
    grid.point$detection<-rep(0,n)
    grid.point$detection[order(-grid.point$Observed)[1:location]]<-1
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

library(doParallel)
cl=makeCluster(3)
registerDoParallel(cl)
# aa=foreach(i=1:500,.combine=rbind,.packages=c("rootSolve","MASS","RandomFields","geoR"),.export="pp",.errorhandling="remove")%dopar%
aa1=foreach(m_time=1:100,.combine=rbind,
            .packages="RandomFields")%dopar%
            {
              one_iter(30,B_inte,0.05,0.22,1,1)
            }
stopCluster(cl)

m1<-colMeans(aa1)>0.3
N.grid<-30
grid.point = expand.grid(x=seq(0,1,by=1/N.grid),y=seq(0,1,by=1/N.grid))
grid.point$CUSUM<-m1[1:961]
aa<-ggplot(grid.point)+geom_point(aes(x=x,y=y,colour=CUSUM))
grid.point$FDRl<-m1[962:1922]
bb<-ggplot(grid.point)+geom_point(aes(x=x,y=y,colour=FDRl))
grid.point$FDR<-m1[1923:2883]
cc<-ggplot(grid.point)+geom_point(aes(x=x,y=y,colour=FDR))

library(doParallel)
cl=makeCluster(3)
registerDoParallel(cl)
# aa=foreach(i=1:500,.combine=rbind,.packages=c("rootSolve","MASS","RandomFields","geoR"),.export="pp",.errorhandling="remove")%dopar%
aa2=foreach(m_time=1:100,.combine=rbind,
            .packages="RandomFields")%dopar%
            {
              one_iter(30,B_inte,0.05,0.22,1,2)
            }
stopCluster(cl)

m2<-colMeans(aa2)>0.55
N.grid<-30
grid.point = expand.grid(x=seq(0,1,by=1/N.grid),y=seq(0,1,by=1/N.grid))
grid.point$CUSUM<-m2[1:961]
dd<-ggplot(grid.point)+geom_point(aes(x=x,y=y,colour=CUSUM))
grid.point$FDRl<-m2[962:1922]
ee<-ggplot(grid.point)+geom_point(aes(x=x,y=y,colour=FDRl))
grid.point$FDR<-m2[1923:2883]
ff<-ggplot(grid.point)+geom_point(aes(x=x,y=y,colour=FDR))



library(doParallel)
cl=makeCluster(3)
registerDoParallel(cl)
# aa=foreach(i=1:500,.combine=rbind,.packages=c("rootSolve","MASS","RandomFields","geoR"),.export="pp",.errorhandling="remove")%dopar%
aa3=foreach(m_time=1:100,.combine=rbind,
            .packages="RandomFields")%dopar%
            {
              one_iter(30,B_inte,0.05,0.22,1,4)
            }
stopCluster(cl)

m4<-colMeans(aa3)>0.9
N.grid<-30
grid.point = expand.grid(x=seq(0,1,by=1/N.grid),y=seq(0,1,by=1/N.grid))
grid.point$CUSUM<-m4[1:961]
gg<-ggplot(grid.point)+geom_point(aes(x=x,y=y,colour=CUSUM))
grid.point$FDRl<-m4[962:1922]
hh<-ggplot(grid.point)+geom_point(aes(x=x,y=y,colour=FDRl))
grid.point$FDR<-m4[1923:2883]
ii<-ggplot(grid.point)+geom_point(aes(x=x,y=y,colour=FDR))


library(grid)
grid.newpage()  ##newpage
pushViewport(viewport(layout = grid.layout(3,3))) ### divide into 1*3
vplayout <- function(x,y){
  viewport(layout.pos.row = x, layout.pos.col = y)
}
print(aa, vp = vplayout(1,1))   ### put a in (1,1)
print(bb, vp = vplayout(1,2))   ### put b in (1,2)
print(cc, vp = vplayout(1,3))  ###put c in (1,3)
print(dd, vp = vplayout(2,1))  ###put c in (1,3)
print(ee, vp = vplayout(2,2))   ### put a in (1,1)
print(ff, vp = vplayout(2,3))   ### put b in (1,2)
print(gg, vp = vplayout(3,1))  ###put c in (1,3)
print(hh, vp = vplayout(3,2))  ###put c in (1,3)
print(ii, vp = vplayout(3,3))  ###put c in (1,3)