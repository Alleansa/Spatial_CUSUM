library(BalancedSampling)
library(doParallel)
library(MASS)
library(geoR)
library(rootSolve)
library(RandomFields)
library(parallel)
library(ggplot2)

N.grid = 100
r = 5 # This is the range of the semivariogram model
N.sample.size = 1000
N.simu = 1000
block.size = c(2,4,5,10,20) # It is better that the block size divides N.grid.


grid.point = expand.grid(x=1:N.grid,y=1:N.grid)

true.theta = c(1,10,0) #the range is about half of the longest distance that two points can achieve in the unit square.
model <- RMspheric(var=true.theta[1], scale=true.theta[2])

simu.alti <- RFsimulate(model, x=grid.point$x,y=grid.point$y,grid = FALSE)
# sample.n$value = simu.alti@data$variable1 + beta0 + beta1 * sample.n$x + beta2 * sample.n$y
grid.point$value = simu.alti@data$variable1
