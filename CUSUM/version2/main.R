rm(list=ls())

source("/Users/Allen/Documents/Study/Ph.D/Code/R/CUSUM/version2/FDRL.R")
source("/Users/Allen/Documents/Study/Ph.D/Code/R/CUSUM/version2/CUSUM.R")


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
  grid.point$nn_mean<-neigh_mean(grid.point$Observed,Dist_matrix,r)
  #test whether there is anomaly
  test_seq<-neigh_info[order(-grid.point$nn_mean)]
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
    grid.point$detection[order(-grid.point$nn_mean)[1:location]]<-1
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
cl=makeCluster(3)
registerDoParallel(cl)
aa1=foreach(m_time=1:2,.combine=rbind,
           .packages=c("RandomFields","foreach"))%dopar%
           {
             multi_iter(N.grid=20,r=0.06,1,times =5)
           }
stopCluster(cl)
