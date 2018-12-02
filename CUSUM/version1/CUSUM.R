#constructing the distance matrix
get_distance <-function(x,y){
  x_matrix=matrix(rep(x,time=length(x)),length(x),length(x))-t(matrix(rep(x,time=length(x)),length(x),length(x)))
  y_matrix=matrix(rep(y,time=length(y)),length(y),length(y))-t(matrix(rep(y,time=length(y)),length(y),length(y)))
  D=sqrt(x_matrix^2+y_matrix^2)
  return(D)
}

#computing the CUSUM for series {e}
CUSUM<-function(e){
  e_star<-c(0,e)
  T<-cumsum(e_star)-mean(e)*c(0:length(e))
  return(T^2/length(e))
}

# neigh mean
neigh_mean<-function(e,D,r){
  return((D>0&D<=r)%*%matrix(e,length(e),1)/rowSums((D>0&D<=r)))
}

#Bridge motion^2
Bridge2prod<-function(N=100000,tim=1000){
  library(Bolstad2)
  library(e1071)
  B_inte <- rep(0,N)
  for (i in c(1:N)){
    x<-rbridge(frequency=tim)
    tem <- sintegral(seq(0,1,1/(tim-1)),(as.vector(x))^2)
    B_inte[i] <- tem$int
    rm(tem)
  }
  return(B_inte)
}

Bridge2test_dire<-function(y,B_inte){
  p=sum(B_inte>y)/length(B_inte)
  return(p)
}

#boundary check
boundary_detection<-function(detection,boundary,observed){
  num<-sum(boundary)
  new_boundary_detection<-numeric(num)
  diff<-abs(boundary-detection)
  if(sum(diff)==0){
    noise_mean<-mean(observed[detection==0],trim=0.05)
    noise_var<-mad(observed[detection==0])
    new_boundary_detection[abs((observed[boundary==1]-noise_mean)/sqrt(noise_var))<qt(0.975,length(detection)-num)]=0
    detection[boundary==1]<-new_boundary_detection
  } else {
    signal_mean<-mean(observed[detection==1],trim=0.05)
    signal_var<-mad(observed[detection==1])
    noise_mean<-mean(observed[detection==0],trim=0.05)
    noise_var<-mad(observed[detection==0])
    ratio<-sqrt(signal_var/noise_var)*exp(-(observed[boundary==1]-noise_mean)^2/noise_var+((observed[boundary==1]-signal_mean)^2/signal_var))
    new_boundary_detection[ratio<1]<-1
    new_boundary_detection[ratio>=1]<-0
    detection[boundary==1]<-new_boundary_detection
  }
  return(detection)
}

#get robust variance
rb_var<-function(e,alpha=0.5){
  n<-length(e)
  v<-var(e[round(alpha*n):n])
  return(v)
}

#main 1
region_detection <- function(data,p_thre=0.01,r=0.1, B_inte){
  library(ggplot2)
  n<-length(data$Observed)
  Dist_matrix <- get_distance(data$x,data$y)
  data$k_neigh<-neigh_mean(data$Observed,Dist_matrix,r)
  #test whether there is anomaly
  test_seq<-grid.point$Observed[order(-data$k_neigh)]
  p<-Bridge2test_dire(mean(CUSUM(test_seq)),B_inte*mad(test_seq))
  if(p>p_thre) {
    print('There is no anomaly in this region')
  }else {
    T<-CUSUM(test_seq)
    location<-which.max(T)
    data$detection<-rep(0,n)
    data$detection[order(-data$k_neigh)[1:location]]<-1
    #boundary check
    data$boundary_score<-neigh_mean(data$detection,Dist_matrix,r)
    data$boundary<-rep(0,n)
    data$boundary[abs(data$boundary_score-0.5)<0.45]<-1
    data$detection<-boundary_detection(data$detection,data$boundary,data$Observed)
    data$p<-rep(0,n)
    data$p[data$detection==1]<-data$k_neigh[data$detection==1]
  }
  return(data)
}