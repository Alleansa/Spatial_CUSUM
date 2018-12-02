
detection<-numeric(1000)

for (i in c(1:200)){
  set<-matrix(rnorm(4000),1000,4)
  sample<-numeric(1000)
  sum<-numeric(1000)
  for (j in c(1:1000)){
    sample[j]<-sample(set[j,],1)
    sum[j]<-sum(set[j,])-sample[j]
  }
  seq<-sample[order(sum)]
  location<-locate.change(seq)[1]

  detection[order(sum)[1:location]]<-detection[order(sum)[1:location]]+1
}


detection<-numeric(1000)
p<-numeric((1000))
set<-matrix(c(rnorm(10000)+1,rnorm(90000)),1000,100,byrow = 1)
for (i in c(1:1000)){
  sample<-numeric(1000)
  sum<-numeric(1000)
  for (j in c(1:1000)){
    sample[j]<-sample(set[j,],1)
    sum[j]<-sum(set[j,])-sample[j]
  }
  seq<-sample[order(-sum)]
  location<-locate.change(seq)[1]
  p[i]<-location/1000
  detection[order(-sum)[1:location]]<-detection[order(-sum)[1:location]]+1
}

p<-numeric((5000))
for (i in c(1:5000)){
  seq<-rnorm(1000)
  location<-locate.change(seq)[1]
  p[i]<-location/1000
}

v<-numeric(200)
for( i in c(1:200)){
  set<-c(rnorm(100)+0.5,rnorm(900))
  v[i]<-mean(CUSUM(set))
}