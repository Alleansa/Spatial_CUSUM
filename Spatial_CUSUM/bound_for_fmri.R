grid.point$detection<-ret==1
grid.point$bound=0
for (i in c(1:128^2)){
  if(grid.point[i,4]==1){
    loc_x<-grid.point[i,1]
    loc_y<-grid.point[i,2]
    tem1<-grid.point[grid.point$x==(loc_x-1)&grid.point$y==loc_y,4]
    tem2<-grid.point[grid.point$x==(loc_x)&grid.point$y==(loc_y-1),4]
    tem3<-grid.point[grid.point$x==(loc_x+1)&grid.point$y==loc_y,4]
    tem4<-grid.point[grid.point$x==(loc_x)&grid.point$y==(loc_y+1),4]
    if(length(tem1)==0){tem1=NA}
    if(length(tem2)==0){tem2=NA}
    if(length(tem3)==0){tem3=NA}
    if(length(tem4)==0){tem4=NA}
    tem<-c(tem1,tem2,tem3,tem4)
    tem<-tem[!is.na(tem)]
    score=mean(tem)
    if(score!=1){grid.point$bound[i]=1}
  }
}

