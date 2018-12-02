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

# neigh info pick
sample_nozero<-function(I,seq,size){
  valid<-t(I*seq)
  return(sample(valid,size,prob=I/sum(I),replace = 1))
}


neigh_pick<-function(e,D,r,k){

    Idct<-D<=r&D>0
    Sample<-apply(Idct,1,sample_nozero,seq=e,size=k)
    return(Sample)
  
  
}


neigh_mean<-function(e,D,r){
  return((D>0&D<=r)%*%matrix(e,length(e),1)/rowSums((D>0&D<=r)))
}

locate.change<-function (x) 
{
  x <- as.matrix(x)
  if (dim(x)[2] == 1) 
    x <- t(x)
  n <- dim(x)[2]
  x<-x/mad(x)
  T <-CUSUM(x)
  changepoint <- which.max(T)
  cusum <- mean(T)
  return(c(changepoint,cusum))
}


#main 1
single_region_detection <- function(data,r=0.1,k=1, threshold=0.47){
  n<-length(data$Observed)
  Dist_matrix <- get_distance(data$x,data$y)
  neigh_sample<-neigh_pick(data$Observed,Dist_matrix,r,k)
  data$nn_mean<-neigh_mean(data$Observed,Dist_matrix,r)
  neigh_info <-matrix(rep(1/k,k),1,k)%*%neigh_sample
  #test whether there is anomaly
  test_seq<-neigh_info[order(-data$nn_mean)]
  if(mean(CUSUM(test_seq/mad(test_seq)))<=threshold) {
    print('There is no anomaly in this region')
  }else {
    location<-locate.change(test_seq)[1]
    data$detection<-rep(0,n)
    data$detection[order(-neigh_info)[1:location]]<-1
    data$p<-rep(0,n)
    data$p[data$detection==1]<-data$Observed[data$detection==1]
  }
  return(data)
}

multi_region_detection<-function(data,r=0.1,k=1, threshold=1.90,M=0){
  n<-length(data$Observed)
  Dist_matrix <- get_distance(data$x,data$y)
  data$nn_mean<-neigh_mean(data$Observed,Dist_matrix,r)
  neigh_sample<-neigh_pick(data$Observed,Dist_matrix,r,k)
  neigh_info <-matrix(rep(1/k,k),1,k)%*%neigh_sample
  test_seq<-neigh_info[order(-data$nn_mean)]
  rnd1 <- sample(0:(n - 2), M, replace = TRUE)
  rnd2 <- sample(0:(n - 2), M, replace = TRUE)
  window_s <- pmin(rnd1, rnd2)
  window_e <- pmax(rnd1, rnd2) + 2
  BinSeg <- function(x, s, e, depth, parent.val) {
    if (e - s <= 2) 
      return(NULL)
    ind <- (window_s >= s) & (window_e <= e)
    max.val <- -1
    cp <- 0
    for (m in c(0, ((1:M)[ind]))) {
      if (m == 0) {
        s_m <- s
        e_m <- e
      }else {
        s_m <- window_s[m]
        e_m <- window_e[m]
      }
      obj <- locate.change(x[(s_m + 1):e_m])
      if (obj[2] > max.val) {
        max.val <- obj[2]
        cp <- s_m + obj[1]
      }
    }
    ret <- NULL
    ret$location <- cp
    ret$max.cusum <- min(parent.val, max.val)
    ret$depth <- depth
    if (ret$max.cusum < threshold |depth>log(n)) {
      return(NULL)
    }else {
      return(cbind(BinSeg(x, s, cp, depth+1, ret$max.cusum), 
                   ret, BinSeg(x, cp, e, depth+1, ret$max.cusum)))
    }
  }
  ret <- NULL
  ret$x <- test_seq
  ret$changepoints <- BinSeg(test_seq, 0, n, depth = 1, parent.val = .Machine$double.xmax)
  ret$changepoints <- t(matrix(as.numeric(ret$changepoints), 
                               nrow = 3))
  dete_cp<-ret$changepoints[,1]
  if(length(dete_cp)==0){
    print("there is no anomaly region")
  }else{
    dete_cp<-unique(c(0,dete_cp,n))
    Num_cp<-length(dete_cp)
    data$detection<-rep(0,n)
    for (ii in c(1:(Num_cp-1))){
    data$detection[order(-neigh_info)[(dete_cp[ii]+1):dete_cp[ii+1]]]<-Num_cp-ii-1
    }
  }
  return(data)
}