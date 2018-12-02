count=0
B_inte<-Bridge2prod(sigma=1)
for (i in 1L:500){
  x<-rnorm(2000)
  p<-Bridge2test_dire(get_sum_CUSUM(x),B_inte)
  if(p<0.01){
    count=count+1
  }
}