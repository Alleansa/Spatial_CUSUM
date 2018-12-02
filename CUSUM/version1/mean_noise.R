get_distance <-function(x,y){
  x_matrix=matrix(rep(x,time=length(x)),length(x),length(x))-t(matrix(rep(x,time=length(x)),length(x),length(x)))
  y_matrix=matrix(rep(y,time=length(y)),length(y),length(y))-t(matrix(rep(y,time=length(y)),length(y),length(y)))
  D=x_matrix^2+y_matrix^2
  return(D)
}

get_k_neigh_mean<-function(e,k){
  N<-length(e)
  T<-mean(e[1:k])
  return(T)
}

noise_neigh<-function(e,D,k){
  n<-length(e)
  k_neigh<-rep(0,n)
  for (i in 1:n){
    k_neigh[i]<-get_k_neigh_mean(e[order(D[i,])],k)
  }
  return(k_neigh)
}

mean_noise<-function(data,prop){
  k<-prop*nrow(data)
  Dist_matrix <- get_distance(data$x,data$y)
  new_noise<-noise_neigh(data$Observed,Dist_matrix,k)
  data$Observed<-new_noise
  return(data)
}

make_noise<-function(siz,u_size=100,prop=0.05){
  x<-rep(c(1:u_size),times=u_size)
  y<-rep(c(1:u_size),each=u_size)
  data<-data.frame(x=x,y=y)
  data$Observed<-rnorm(u_size*u_size,sd=sqrt(80))
  data<-mean_noise(data,prop)
  return(data[which(data$x>=(u_size/2-siz/2+1)&data$x<=(u_size/2+siz/2)
                    &data$y>=(u_size/2-siz/2+1)&data$y<=(u_size/2+siz/2)),3])
}

product_noise<-function(siz,u_size,trial,prop){
  mvalue<-matrix(rep(0,trial*(siz*siz)),siz*siz,trial)
  for (t in c(1:trial)){
    mvalue[,t]<-make_noise(siz,u_size,prop)
  }
  return(mvalue)
}
