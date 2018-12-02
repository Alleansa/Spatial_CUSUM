#change r
N.grid=50
k=1
grid.point = expand.grid(x=seq(0,1,by=1/N.grid),y=seq(0,1,by=1/N.grid))
grid.point$Observed<-rnorm((N.grid+1)^2)
n<-length(grid.point$Observed)
Dist_matrix <- get_distance(grid.point$x,grid.point$y)
par(mfcol=c(2,2))
p=rep(0,500)
for(i in c(1:500)){
  neigh_info <-neigh_pick(grid.point$Observed,Dist_matrix,0.35,k)
  #test whether there is anomaly
  test_seq<-neigh_info[order(grid.point$Observed)]
  # p[i]<-sum(B_inte>mean(CUSUM(test_seq)))/length(B_inte)
  p[i]<-sum(B_max>max(CUSUM(test_seq)))/length(B_max)
}
hist(p,main="r=0.1")

p=rep(0,500)
for(i in c(1:500)){
  neigh_info <-neigh_pick(grid.point$Observed,Dist_matrix,0.45,k)
  #test whether there is anomaly
  test_seq<-neigh_info[order(grid.point$Observed)]
  # p[i]<-sum(B_inte>mean(CUSUM(test_seq)))/length(B_inte)
  p[i]<-sum(B_max>max(CUSUM(test_seq)))/length(B_max)
}
hist(p,main='r=0.15')

p=rep(0,500)
for(i in c(1:500)){
  neigh_info <-neigh_pick(grid.point$Observed,Dist_matrix,0.55,k)
  #test whether there is anomaly
  test_seq<-neigh_info[order(grid.point$Observed)]
  # p[i]<-sum(B_inte>mean(CUSUM(test_seq)))/length(B_inte)
  p[i]<-sum(B_max>max(CUSUM(test_seq)))/length(B_max)
}
hist(p,main="r=0.2")

p=rep(0,500)
for(i in c(1:500)){
  neigh_info <-neigh_pick(grid.point$Observed,Dist_matrix,1,k)
  #test whether there is anomaly
  test_seq<-neigh_info[order(grid.point$Observed)]
  # p[i]<-sum(B_inte>mean(CUSUM(test_seq)))/length(B_inte)
  p[i]<-sum(B_max>max(CUSUM(test_seq)))/length(B_max)
}
hist(p,main='r=0.25')


#change N.grid
par(mfcol=c(2,2))

p=rep(0,500)
N.grid=10
grid.point = expand.grid(x=seq(0,1,by=1/N.grid),y=seq(0,1,by=1/N.grid))
Dist_matrix <- get_distance(grid.point$x,grid.point$y)
n<-length(grid.point$Observed)
for(i in c(1:500)){
  grid.point$Observed<-rnorm((N.grid+1)^2)
  neigh_info <-neigh_pick(grid.point$Observed,Dist_matrix,0.2,1)
  #test whether there is anomaly
  test_seq<-neigh_info[order(grid.point$Observed)]
  p[i]<-sum(B_inte>mean(CUSUM(test_seq)))/length(B_inte)
}
hist(p,main='N.grid=10')

N.grid=20
grid.point = expand.grid(x=seq(0,1,by=1/N.grid),y=seq(0,1,by=1/N.grid))
Dist_matrix <- get_distance(grid.point$x,grid.point$y)
n<-length(grid.point$Observed)
for(i in c(1:500)){
  grid.point$Observed<-rnorm((N.grid+1)^2)
  neigh_info <-neigh_pick(grid.point$Observed,Dist_matrix,0.2,1)
  #test whether there is anomaly
  test_seq<-neigh_info[order(grid.point$Observed)]
  p[i]<-sum(B_inte>mean(CUSUM(test_seq)))/length(B_inte)
}
hist(p,main='N.grid=20')

N.grid=30
grid.point = expand.grid(x=seq(0,1,by=1/N.grid),y=seq(0,1,by=1/N.grid))
Dist_matrix <- get_distance(grid.point$x,grid.point$y)
n<-length(grid.point$Observed)
for(i in c(1:500)){
  grid.point$Observed<-rnorm((N.grid+1)^2)
  neigh_info <-neigh_pick(grid.point$Observed,Dist_matrix,0.2,1)
  #test whether there is anomaly
  test_seq<-neigh_info[order(grid.point$Observed)]
  p[i]<-sum(B_inte>mean(CUSUM(test_seq)))/length(B_inte)
}
hist(p,main='N.grid=20')

N.grid=40
grid.point = expand.grid(x=seq(0,1,by=1/N.grid),y=seq(0,1,by=1/N.grid))
Dist_matrix <- get_distance(grid.point$x,grid.point$y)
n<-length(grid.point$Observed)
for(i in c(1:500)){
  grid.point$Observed<-rnorm((N.grid+1)^2)
  neigh_info <-neigh_pick(grid.point$Observed,Dist_matrix,0.2,1)
  #test whether there is anomaly
  test_seq<-neigh_info[order(grid.point$Observed)]
  p[i]<-sum(B_inte>mean(CUSUM(test_seq)))/length(B_inte)
}
hist(p,main='N.grid=40')
