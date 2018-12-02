kClusters <- 2
kMeans <- kmeans(data$Observed, centers = kClusters)
data$kmeans<-kMeans$cluster
ggplot() + geom_point(data=data, aes(x = y, y = rev(x),colour=kmeans))
