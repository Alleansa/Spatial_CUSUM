#Input data frame: data$id sample id
#                  data$x x-axis
#                  data$y y-axis
#                  data$Observed the observing value
                  
##Apr.3 2017
#plot the scatter point
my_plot<-function(data){
  library(ggplot2)
  ggplot(data, aes(x = x, y = y,colour=Observed)) + geom_point()
}

#constructing the distance matrix
get_distance <-function(x,y){
  x_matrix=matrix(rep(x,time=length(x)),length(x),length(x))-t(matrix(rep(x,time=length(x)),length(x),length(x)))
  y_matrix=matrix(rep(y,time=length(y)),length(y),length(y))-t(matrix(rep(y,time=length(y)),length(y),length(y)))
  D=x_matrix^2+y_matrix^2
  return(D)
}


#computing the k_th order CUSUM for series {e}
get_k_CUSUM<-function(e,k){
  N<-length(e)
  T<-(sum(e[1:k])-k*mean(e))^2/N
  return(T)
}

#computing the average of k_th nearest points
get_k_neigh<-function(e,k){
  N<-length(e)
  T<-mean(e[2:k])
  return(T)
}

get_sum_CUSUM <- function(e){
  S=0
  N=length(e)
  for (i in c(1:N)){
    S=S+get_k_CUSUM(e,i)
  }
  return(S/N)
}



#generating square of brownian bridge at time t with variance sigma^2
Bridge2<-function(sigma=1){
  # N is the number of end-points of the grid including T
  N<-99
  T <- 1
  # length of the interval [0, T] in time units
  Delta <- T/N
  # time increment
  W <- numeric(N+1)
  # initialization of the vector W approximating 
  # Wiener process
  t <- seq(0,T, length=N+1)
  W <- c(0, cumsum( sqrt(Delta) * rnorm(N)))
  B <- (W-W[100]*t)^2*sigma^2
  return(B)
}

#computing p-value under H0
Bridge2test <- function(y,sigma=1){
  library(Bolstad2)
  N<-10000
  Bmatrix <- matrix(rep(0,10000*5000),10000,5000)
  B_inte <- rep(0,10000)
  for (i in c(1:10000)){
    Bmatrix[i,] <- Bridge2(sigma)
    tem <- sintegral(seq(0,1,1/4999),Bmatrix[i,])
    B_inte[i] <- tem$int
    rm(tem)
  }
  p=sum(B_inte>y)/10000
  return(p)
}

Bridge2prod<-function(sigma=1){
  library(Bolstad2)
  N<-10000
  Bmatrix <- matrix(rep(0,10000*5000),10000,5000)
  B_inte <- rep(0,10000)
  for (i in c(1:10000)){
    Bmatrix[i,] <- Bridge2(sigma)
    tem <- sintegral(seq(0,1,1/4999),Bmatrix[i,])
    B_inte[i] <- tem$int
    rm(tem)
  }
  return(B_inte)
}

Bridge2test_dire<-function(y,B_inte){
  p=sum(B_inte>y)/10000
  return(p)
}


#reorder the samples
reorder<-function(e,D,k){
  n<-length(e)
  k_cusum<-rep(0,n)
  for (i in 1:n){
    k_cusum[i]<-get_k_CUSUM(e[order(D[i,])],k)
  }
  return(k_cusum)
}

# neigh mean
neigh_mean<-function(e,D,k){
  n<-length(e)
  k_neigh<-rep(0,n)
  for (i in 1:n){
    k_neigh[i]<-get_k_neigh(e[order(D[i,])],k)
  }
  return(k_neigh)
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
region_detection <- function(data,p_thre=0.05,k=5){
  library(ggplot2)
  n<-length(data$Observed)
  #get the distance matrix
  Dist_matrix <- get_distance(data$x,data$y)
  #data$k_cusum<-reorder(data$Observed,Dist_matrix,k)
  data$k_neigh<-neigh_mean(data$Observed,Dist_matrix,k)
  #test whether there is anomaly
  #test_seq<-data$Observed[order(-data$k_cusum)]
  test_seq<-data$Observed[order(-data$k_neigh)]
  p<-Bridge2test(get_sum_CUSUM(test_seq),sigma=mad(test_seq))
  #p<-Bridge2test_dire(get_sum_CUSUM(test_seq),B_inte)
  if(p>p_thre) {
    print('There is no anomaly in this region')
  }else {
    T<-rep(0,n)
    for (i in 1:n){
      T[i]<-get_k_CUSUM(test_seq,i)
    }
    location<-which.max(T)
    data$detection<-rep(0,n)
    data$detection[order(-data$k_neigh)[1:location]]<-1
    #boundary check
    data$boundary_score<-neigh_mean(data$detection,Dist_matrix,k)
    data$boundary<-rep(0,n)
    data$boundary[abs(data$boundary_score-0.5)<0.45]<-1
    data$detection<-boundary_detection(data$detection,data$boundary,data$Observed)
  }
  ggplot()+geom_point(data=data, aes(x = x, y = y,colour=detection))
  return(data)
}

# #Apr.15 2017
# 
# #get neigh.
# get_neigh<-function(data,id){
#   x<-data$x[id]
#   y<-data$y[id]
#   neigh<-rbind(data[which(data$x==x-1&data$y==y),],
#                data[which(data$x==x+1&data$y==y),],
#                data[which(data$x==x&data$y==y-1),],
#                data[which(data$x==x&data$y==y+1),]
#                )
#   return(neigh)
# }
# #remove point
# rm_point<-function(old_neigh,new_neigh,rm_list){
#   dif_list<-setdiff(new_neigh$id,rm_list)
#   neigh<-rbind(old_neigh,new_neigh[new_neigh$id %in% dif_list,])
#   return(neigh)
# }
# #get the closest point
# add_point<-function(neigh,mean){
#   abs_dif<-abs(neigh$Observed-mean)
#   return(neigh$id[which.min(abs_dif)])
# }
# 
# 
# #main2
# tree_detection<-function(data,p_thre=0.05){
#   n<-dim(data)[1]
#   data$id<-seq(1,n)
#   #get the start point with largest S1
#   mean_Observed<-mean(data$Observed)
#   diff<-abs(data$Observed-mean_Observed)
#   rm_list <- which.max(diff)
#   new_order<-rm_list
#   neigh<-get_neigh(data,rm_list)
#   rm_list<-c(rm_list,neigh$id)
#   for (i in 1:n){
#     new_id<-add_point(neigh,mean(data[new_order,]$Observed))
#     new_order<-c(new_order,new_id)
#     new_neigh<-get_neigh(data,new_id)
#     neigh<-rm_point(neigh,new_neigh,rm_list)
#     neigh<-neigh[-which(neigh$id==new_id),]
#     rm_list<-unique(c(rm_list,new_neigh$id))
#   }
#   test_seq<-data$Observed[new_order]
#   p<-tree_test(get_sum_CUSUM(test_seq))
#   if(p>p_thre) {
#     print('There is no anomaly in this region')
#   }else {
#     T<-rep(0,n)
#     for (i in 1:n){
#       T[i]<-get_k_CUSUM(test_seq,i)
#     }
#     location<-which.max(T)
#     data$detection<-rep(1,n)
#     data$detection[new_order[1:location]]<-2
#   }
#   ggplot()+geom_point(data=data, aes(x = x, y = y,colour=detection))
#   return(data)
# }
# 
# 
# #April.20 1017
# generating_tree<-function(data){
#   n<-dim(data)[1]
#   data$id<-seq(1,n)
#   #get the start point with largest S1
#   mean_Observed<-mean(data$Observed)
#   diff<-abs(data$Observed-mean_Observed)
#   rm_list <- which.max(diff)
#   new_order<-rm_list
#   neigh<-get_neigh(data,rm_list)
#   rm_list<-c(rm_list,neigh$id)
#   for (i in 1:n){
#     new_id<-add_point(neigh,mean(data[new_order,]$Observed))
#     new_order<-c(new_order,new_id)
#     new_neigh<-get_neigh(data,new_id)
#     neigh<-rm_point(neigh,new_neigh,rm_list)
#     neigh<-neigh[-which(neigh$id==new_id),]
#     rm_list<-unique(c(rm_list,new_neigh$id))
#   }
#   seq<-data$Observed[new_order]
#   return(seq)
# }
# 
# random_tree_seq<-function(n=10){
#   x<-rep(c(1:n),times=n)
#   y<-rep(c(1:n),each=n)
#   data<-data.frame(x=x,y=y)
#   data$Observed=runif(n*n)
#   y<- get_sum_CUSUM(generating_tree(data))
#   return(y)
# }
# 
# tree_test<-function(y,times=500,sigma=1, n = 10){
#   CV<-rep(0,times)
#   for (i in c(1:times)){
#      CV[i]<-random_tree_seq(n)*sigma^2
#   }
#   p=sum(CV>y)/times
#   return(p)
# }
# 
# 
