load("/Users/Allen/Desktop/B_inte.rda")

source("/Users/Allen/Documents/Study/Ph.D/Code/R/CUSUM/FDRL.R")
source("/Users/Allen/Documents/Study/Ph.D/Code/R/CUSUM/CUSUM.R")

one_iter<-function(N.grid,B_inte,p_thre,r){
  fail_count<-0
  miss_report<-0
  wrong_report<-0
  fdr_miss_report<-0
  fdr_wrong_report<-0
  Status<-rep(0,(N.grid+1)^2)
  grid.point = expand.grid(x=seq(0,1,by=1/N.grid),y=seq(0,1,by=1/N.grid))
  
  # #nondependence
  # grid.point$Observed<-rnorm((N.grid+1)^2)
  
  #dependence
  model <- RMspheric(var=1, scale=0.1)
  simu.alti <- RFsimulate(model, x=grid.point$x,y=grid.point$y,grid = FALSE)
  grid.point$Observed = simu.alti@data$variable1
  

  #generating 'L' shape
  grid.point$Status<-Status
  grid.point[which(grid.point$x<=0.8&grid.point$x>=.4&grid.point$y<=.8&grid.point$y>=.4),4]<-1
  grid.point[which(grid.point$x<=0.8&grid.point$x>=.6&grid.point$y<=.8&grid.point$y>=.6),4]<-0
  grid.point[grid.point$Status==1,]$Observed<-grid.point[grid.point$Status==1,]$Observed+5
  pvalues<-pvalue_normal(grid.point$Observed)
  
  # n<-length(grid.point$Observed)
  #get the distance matrix
  n<-length(grid.point$Observed)
  Dist_matrix <- get_distance(grid.point$x,grid.point$y)
  grid.point$k_neigh<-neigh_mean(grid.point$Observed,Dist_matrix,r)
  #test whether there is anomaly
  test_seq<-grid.point$Observed[order(-grid.point$k_neigh)]
  p<-Bridge2test_dire(mean(CUSUM(test_seq)),B_inte)
  if(p>p_thre) {
    #print('There is no anomaly in this region')
    fail_count<-1
    
    miss_report<-sum(grid.point$Status)
    aprox_fdr<-match_wrong(pvalues,grid.point$Status,grid.point$x,grid.point$y,miss_report,window=r)
    fdr_wrong_report<-aprox_fdr[1]
    fdr_miss_report<-aprox_fdr[2]
    

    # res<-FDRLMethod(pvalues,grid.point$x,grid.point$y,window=r,alpha=p_thre)
    # fdr_wrong_report<-sum(res$ind)
  }else {
    T<-CUSUM(test_seq)
    location<-which.max(T)
    grid.point$detection<-rep(0,n)
    grid.point$detection[order(-grid.point$k_neigh)[1:location]]<-1
    grid.point$boundary_score<-neigh_mean(grid.point$detection,Dist_matrix,r)
    grid.point$boundary<-rep(0,n)
    grid.point$boundary[abs(grid.point$boundary_score-0.5)<0.45]<-1
    grid.point$detection<-boundary_detection(grid.point$detection,grid.point$boundary,grid.point$Observed)
    
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


cl=makeCluster(3)
registerDoParallel(cl)
# aa=foreach(i=1:500,.combine=rbind,.packages=c("rootSolve","MASS","RandomFields","geoR"),.export="pp",.errorhandling="remove")%dopar%
aa=foreach(m_time=1:100,.combine=rbind,
           .packages="RandomFields",
           .errorhandling="remove")%dopar%
{
  one_iter(40,B_inte,0.01,0.25)
}
stopCluster(cl)


# #N.grid=10
# mean(aa[,1])/16
# sum(aa[,2])/nrow(aa)
# sum(aa[,3])/((121-16)*(nrow(aa)-sum(aa[,2])))
# 
# mean(aa[,4])/16
# sum(aa[,5]==0)/nrow(aa)
# sum(aa[,5])/((121-16)*(nrow(aa)-sum(aa[,5]==0)))
# 
# # sum(aa[,5]==0&aa[,4]==16)/nrow(aa)
# # sum(aa[,5])/((121-16)*(nrow(aa)-sum(aa[,5]==0&aa[,4]==16)))

# #N.grid=25
# mean(aa[,1])/85
# sum(aa[,2])/nrow(aa)
# sum(aa[,3])/((676-85)*(nrow(aa)-sum(aa[,2])))
# 
# mean(aa[,4])/85
# sum(aa[,5]==0)/nrow(aa)
# sum(aa[,5])/((676-85)*(nrow(aa)-sum(aa[,5]==0)))
# 
# # sum(aa[,5]==0&aa[,4]==85)/nrow(aa)
# # sum(aa[,5])/((676-85)*(nrow(aa)-sum(aa[,5]==0&aa[,4]==85)))

#N.grid=40
mean(aa[,1])/208
sum(aa[,2])/nrow(aa)
sum(aa[,3])/((1671-208)*(nrow(aa)-sum(aa[,2])))

mean(aa[,4])/208
# sum(aa[,5]==0)/nrow(aa)
# sum(aa[,5])/((1671-208)*(nrow(aa)-sum(aa[,5]==0)))

sum(aa[,5]==0&aa[,4]==208)/nrow(aa)
sum(aa[,5])/((1671-208)*(nrow(aa)-sum(aa[,5]==0&aa[,4]==208)))