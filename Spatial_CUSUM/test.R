N.grid=100
shift=0.5
scale=0
times=50
k=3
alpha=0.05

Status<-rep(0,(N.grid)^2)
grid.point = expand.grid(x=seq(1,N.grid,by=1),y=seq(1,N.grid,by=1))
if(scale==0){
  #nondependence
  grid.point$Observed<-rnorm((N.grid)^2)
}else{
  #dependence
  model <- RMspheric(var=1, scale=scale)
  simu.alti <- RFsimulate(model, x=grid.point$x,y=grid.point$y,grid = FALSE)
  grid.point$Observed = simu.alti@data$variable1
}

grid.point$Status<-Status
grid.point[which(grid.point$x<=0.5*N.grid&grid.point$x>=.2*N.grid&grid.point$y<=.5*N.grid&grid.point$y>=.2*N.grid),4]<-1
grid.point[which(grid.point$x<=0.5*N.grid&grid.point$x>=.3*N.grid&grid.point$y<=.5*N.grid&grid.point$y>=.3*N.grid),4]<-0

grid.point[which(grid.point$x<=0.9*N.grid&grid.point$x>=.6*N.grid&grid.point$y<=.9*N.grid&grid.point$y>=.6*N.grid),4]<-1
grid.point[which(grid.point$x<=0.78*N.grid&grid.point$x>=.73*N.grid&grid.point$y<=.9*N.grid&grid.point$y>=.85*N.grid),4]<-0
grid.point[which(grid.point$x<=0.78*N.grid&grid.point$x>=.73*N.grid&grid.point$y<=.65*N.grid&grid.point$y>=.6*N.grid),4]<-0
grid.point[which(grid.point$x<=0.8*N.grid&grid.point$x>=.7*N.grid&grid.point$y<=.8*N.grid&grid.point$y>=.7*N.grid),4]<-0


grid.point[grid.point$Status==1,]$Observed<-grid.point[grid.point$Status==1,]$Observed+shift


detection<-numeric(dim(grid.point)[1])
p<-0
for (xi in c((1-k):0)){
  for (yi in c((1-k):0)){ 
    for (i in c(1:times)){
      grid.point <- within(grid.point, {
        grp_x <- cut(x, seq(xi,N.grid+k,by=k), labels = FALSE)
        grp_y <- cut(y, seq(yi,N.grid+k,by=k), labels = FALSE)
      })
      grid.point$group_label<-grid.point$grp_x+(grid.point$grp_y-1)*max(grid.point$grp_x)
      group_sample<-aggregate(grid.point$Observed,by=list(grid.point$group_label),FUN=function(x){sample(x,1)})
      group_sum<-aggregate(grid.point$Observed,by=list(grid.point$group_label),FUN=sum)
      group_count<-aggregate(rep(1,dim(grid.point)[1]),by=list(grid.point$group_label),FUN=sum)
      group<-cbind(group_sample,(group_sum$x-group_sample$x)/(group_count$x))
      names(group)<-c("group","sample","prob")
      rm(group_sample,group_sum,group_count)
      location<-locate.change(group$sample[order(-group$prob)])[1]
      signif_group<-order(-group$prob)[1:location]
      detection[grid.point$group_label %in% signif_group]=1+detection[grid.point$group_label %in% signif_group]
      # p<-c(p,sum(grid.point$Status==0& grid.point$group_label %in% signif_group)/(10000-sum(grid.point$Status)))
    }
  }
}
detection<-detection/(times*k^2)

est <- lpdensity(data = detection,grid=x_seq, bwselect = "imse-rot")$Estimate[,4]
est[est<0]=0
location<-as.numeric(names(which.min(est[x_seq>=quantile(detection,0.5)])))
h0<-c(est[1:location],seq(est[location],0,length.out=length(x_seq)-location))
ratio<-rev(cumsum(rev(h0))/cumsum(rev(est)))
grid.point$detection=detection>=(x_seq[which.max(ratio<0.05)])

res<-grid.point$Status-grid.point$detection
wrong_report<-length(res[res==-1])
miss_report<-length(res[res==1])


#fdrl_method
pvalues<-pvalue_normal(grid.point$Observed)
fdrl<-FDRLMethod(pvalues,grid.point$x,grid.point$y,window=(k)/2,alpha)
fdrl_Res<-grid.point$Status-fdrl$ind
fdrl_wrong<-length(fdrl_Res[fdrl_Res==-1])
fdrl_miss<-length(fdrl_Res[fdrl_Res==1])

c(miss_report/sum(grid.point$Status),wrong_report/sum(1-grid.point$Status),wrong_report/sum(grid.point$detection),
  fdrl_miss/sum(grid.point$Status),fdrl_wrong/sum(1-grid.point$Status),fdrl_wrong/fdrl$numAlt)

