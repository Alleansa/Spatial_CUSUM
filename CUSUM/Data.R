
library(AnalyzeFMRI)
xxx <- f.read.analyze.volume( "CNTS_inter_treg+righthand-corr.img")

# And then you will get a four-dimensional array of the above dimensions.

#The files are on correlations. So, we will variance-stabilize this in
#order to bring it to the normal, since that is what keeps life simple).
#We use Fisher's z-transformation on  the dataset. Thus

#yyy <- atanh(xxx)

#is the dataset of dimensions 128x128x22x12 where 12
#is the number of replications, 22 the slices and the x-y plane is 128x128.
# Let us scale it to be

yyy <- sqrt(121) * atanh(xxx)

#so that we are dealing with unit-variance (marginally) normal random variables.

#zzz[, , , 1] is the first dataset which is what we could stary with for now.

### plot
library(akima)
library(fields)

### look at individual:
pdf("fMRI.case1.pdf")
par(mfrow=c(2,3))
for (i in 1:22) plot.surface(list(x=seq(0,1,l=128), y=seq(0,1,l=128), z=yyy[,,i,1]), type="I", zlim=c(-7,15))
dev.off()

pdf("fMRI.case5.pdf")
par(mfrow=c(2,3))
for (i in 1:22) plot.surface(list(x=seq(0,1,l=128), y=seq(0,1,l=128), z=yyy[,,i,5]), type="I", zlim=c(-7,15))
dev.off()

### save the first case
ycase1 <- yyy[,,,1]
save(ycase1, file= "fMRI.case1.RData")

### look at different slice
pdf("fMRI.22.pdf")
par(mfrow=c(2,3))
for (i in 1:12) image(yyy[,,22,i],col = rainbow(n=128) )
dev.off()
pdf("fMRI.22d.pdf")
par(mfrow=c(2,3))
for (i in 1:12) image(yyy[,,22,i]-yyy.av[,,22],col = rainbow(n=128) )
dev.off()

### look at difference to average over 12 cases:
yyy.av <- apply(yyy, c(1,2,3), mean)
yyy.var <- apply(yyy, c(1,2,3), var)
pdf("fMRI.var.pdf")
par(mfrow=c(2,3))
for (i in 1:22) image(yyy.var[,,i],col = rainbow(n=128) )
dev.off()
pdf("fMRI.av.pdf")
par(mfrow=c(2,3))
for (i in 1:22) image(yyy.av[,,i],col = rainbow(n=128) )
dev.off()




### slice 8 low variation, one high spot. 9 high variation
par(mfrow=c(2,2))
for (i in 8:9) image(yyy.av[,,i],col = rainbow(n=128) )
for (i in 8:9) image(yyy.var[,,i],col = rainbow(n=128) )
for (i in 1:4) image(yyy[,,8,i]-yyy.av[,,8],col = rainbow(n=128) )



par(mfrow=c(4,6))
for (i in 1:22) image(yyy[,,i,1])


par(mfrow=c(3,4))
for (i in 1:12) image(yyy[,,1,i],col = rainbow(n=128))

for (i in 1:12) image(yyy[,,2,i],col = rainbow(n=128))
for (i in 1:12) hist(yyy[,,2,i])

for (i in 1:12) image(yyy[,,i,1]-yyy[,,i,2],col = rainbow(n=128))

par(mfrow=c(2,2))
summary(yyy[,,3,1],col = rainbow(n=128))




### look at the 3rd slice, average over 12 cases
y0 <- apply(yyy[,,3,], c(1,2), mean)
image(y0, col = rainbow(n=128))
image(yyy[,,3,1], col = rainbow(n=128))
image(yyy[,,3,2], col = rainbow(n=128))
image(yyy[,,3,3], col = rainbow(n=128))


image(y0, col = rainbow(n=128))
image(y0<1)

y1 <- yyy[,,3,1]
image(y1<2.5)

y2 <- yyy[,,3,2]
image(y2<2.5)

y3 <- yyy[,,3,3]
image(y3<2.5)
image((y3-y0)<2.5)
image(y3, col = rainbow(n=128))
image(y3-y0, col = rainbow(n=128))



y2 <- yyy[,,3,3]
y2[abs(yyy[,,3,3])< 2] <- 0
y2[abs(yyy[,,3,3])> 2] <- 1
image(y2)

y1 <- yyy[,,3,2]
y1[yyy[,,3,2]< 1.5] <- 0
y1[yyy[,,3,2]> 1.5] <- 1
image(y1)


summary(yyy[,,3,1]-yyy[,,3,2],col = rainbow(n=128))
