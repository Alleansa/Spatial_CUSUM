# trial=100
# fail_count<-0
# miss_report<-rep(0,trial)
# wrong_report<-rep(0,trial)
# fdr_miss_report<-rep(0,trial)
# fdr_wrong_report<-rep(0,trial)
# p_thre=0.01
# k<-9
# mvalue<-matrix(rep(0,trial*1600),1600,trial)
# 
# 
# for (m_time in c(1:trial)){
# siz<-40
# x<-rep(c(1:siz),times=siz)
# y<-rep(c(1:siz),each=siz)
# Status<-rep(0,siz*siz)
# data<-data.frame(x=x,y=y,Status=Status)
# #generating 'L' shape
# data[which(data$x<=0.8*siz&data$x>=.4*siz&data$y<=.8*siz&data$y>=.4*siz),3]<-1
# data[which(data$x<=0.8*siz&data$x>=.6*siz&data$y<=.8*siz&data$y>=.6*siz),3]<-0
# data$Observed<-rnorm(siz*siz,0,1)
# data[data$Status==1,]$Observed=data[data$Status==1,]$Observed
# pvalues<-pvalue_normal(data$Observed)
# mvalue[,m_time]<-data$Observed
# 
# n<-length(data$Observed)
# #get the distance matrix
# Dist_matrix <- get_distance(data$x,data$y)
# #data$k_cusum<-reorder(data$Observed,Dist_matrix,k)
# data$k_neigh<-neigh_mean(data$Observed,Dist_matrix,k)
# #test whether there is anomaly
# #test_seq<-data$Observed[order(-data$k_cusum)]
# test_seq<-data$Observed[order(-data$k_neigh)]
# p<-Bridge2test(get_sum_CUSUM(test_seq),1)
# if(p>p_thre) {
#   #print('There is no anomaly in this region')
#   fail_count<-fail_count+1
# 
#   # miss_report[m_time]<-sum(data$Status)
#   # aprox_fdr<-match_wrong(pvalues,data$Status,data$x,data$y,miss_report[m_time])
#   # fdr_wrong_report[m_time]<-aprox_fdr[1]
#   # fdr_miss_report[m_time]<-aprox_fdr[2]
# 
# 
#   res<-FDRLMethod(pvalues,data$x,data$y,window=9,alpha=0.01)
#   fdr_wrong_report[m_time]<-sum(res$ind)
# }else {
#   T<-rep(0,n)
#   for (i in 1:n){
#     T[i]<-get_k_CUSUM(test_seq,i)
#   }
#   location<-which.max(T)
#   data$detection<-rep(0,n)
#   data$detection[order(-data$k_neigh)[1:location]]<-1
#   
#   # res<-data$Status-data$detection
#   # wrong_report[m_time]<-length(res[res==-1])
#   # miss_report[m_time]<-length(res[res==1])
#   # aprox_fdr<-match_wrong(pvalues,data$Status,data$x,data$y,miss_report[m_time])
#   # fdr_wrong_report[m_time]<-aprox_fdr[1]
#   # fdr_miss_report[m_time]<-aprox_fdr[2]
# 
#   wrong_report[m_time]<-location
#   fdrres<-FDRLMethod(pvalues,data$x,data$y,window=9,alpha=0.01)
#   fdr_wrong_report[m_time]<-sum(fdrres$ind)
# }
# print(m_time)
# }
# 
B_inte<-Bridge2prod(sigma=1)

siz<-25
x<-rep(c(1:siz),times=siz)
y<-rep(c(1:siz),each=siz)
Status<-rep(0,siz*siz)
data<-data.frame(x=x,y=y,Status=Status)
#generating 'L' shape
data[which(data$x<=0.8*siz&data$x>=.4*siz&data$y<=.8*siz&data$y>=.4*siz),3]<-1
data[which(data$x<=0.8*siz&data$x>=.6*siz&data$y<=.8*siz&data$y>=.6*siz),3]<-0
trial=100
fail_count<-0
miss_report<-rep(0,trial)
wrong_report<-rep(0,trial)
fdr_miss_report<-rep(0,trial)
fdr_wrong_report<-rep(0,trial)
p_thre=0.01
k<-9
mvalue<-matrix(rep(0,trial*(siz*siz)),siz*siz,trial)


for (m_time in c(1:trial)){
  data$Observed<-rnorm(siz*siz,0,1)
  data<-mean_noise(data)
  data[data$Status==1,]$Observed=data[data$Status==1,]$Observed
  pvalues<-pvalue_normal(data$Observed)
  mvalue[,m_time]<-data$Observed
  
  n<-length(data$Observed)
  #get the distance matrix
  Dist_matrix <- get_distance(data$x,data$y)
  #data$k_cusum<-reorder(data$Observed,Dist_matrix,k)
  data$k_neigh<-neigh_mean(data$Observed,Dist_matrix,k)
  #test whether there is anomaly
  #test_seq<-data$Observed[order(-data$k_cusum)]
  test_seq<-data$Observed[order(-data$k_neigh)]
  p<-Bridge2test_dire(get_sum_CUSUM(test_seq),B_inte)
  if(p>p_thre) {
    #print('There is no anomaly in this region')
    fail_count<-fail_count+1
    
    # miss_report[m_time]<-sum(data$Status)
    # aprox_fdr<-match_wrong(pvalues,data$Status,data$x,data$y,miss_report[m_time],window=k)
    # fdr_wrong_report[m_time]<-aprox_fdr[1]
    # fdr_miss_report[m_time]<-aprox_fdr[2]

    
    res<-FDRLMethod(pvalues,data$x,data$y,window=k,alpha=0.01)
    fdr_wrong_report[m_time]<-sum(res$ind)
  }else {
    T<-rep(0,n)
    for (i in 1:n){
      T[i]<-get_k_CUSUM(test_seq,i)
    }
    location<-which.max(T)
    data$detection<-rep(0,n)
    data$detection[order(-data$k_neigh)[1:location]]<-1
    data$boundary_score<-neigh_mean(data$detection,Dist_matrix,k)
    data$boundary<-rep(0,n)
    data$boundary[abs(data$boundary_score-0.5)<0.45]<-1
    data$detection<-boundary_detection(data$detection,data$boundary,data$Observed)
    
    # res<-data$Status-data$detection
    # wrong_report[m_time]<-length(res[res==-1])
    # miss_report[m_time]<-length(res[res==1])
    # aprox_fdr<-match_wrong(pvalues,data$Status,data$x,data$y,miss_report[m_time],window = k)
    # fdr_wrong_report[m_time]<-aprox_fdr[1]
    # fdr_miss_report[m_time]<-aprox_fdr[2]

    
    # res<-data$Status-data$detection
    # wrong_report[m_time]<-length(res[res==-1])
    # miss_report[m_time]<-length(res[res==1])
    # fdrres<-FDRLMethod(pvalues,data$x,data$y,window=k,alpha=0.01)
    # fdr_res<-data$Status-fdrres$ind
    # fdr_wrong_report[m_time]<-length(res[fdr_res==-1])
    # fdr_miss_report[m_time]<-length(res[fdr_res==1])

    wrong_report[m_time]<-sum(data$detection)
    fdrres<-FDRLMethod(pvalues,data$x,data$y,window=k,alpha=0.01)
    fdr_wrong_report[m_time]<-sum(fdrres$ind)
  }
  print(m_time)
}

fdr_fail<-sum(fdr_miss_report==16&fdr_wrong_report==0)


fail_count
fdr_fail


mean(miss_report)/85
sum(wrong_report)/((100-fail_count)*540)

mean(fdr_miss_report)/85
sum(fdr_wrong_report)/((100-fdr_fail)*540)

