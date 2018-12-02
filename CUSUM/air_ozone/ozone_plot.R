iowa<-read.csv(file="/Users/Allen/Documents/Research/CODE/R/CUSUM/air_ozone/iowa.csv",header=1)
illinois<-read.csv(file="/Users/Allen/Documents/Research/CODE/R/CUSUM/air_ozone/illinois.csv",header=1)
minnesota<-read.csv(file="/Users/Allen/Documents/Research/CODE/R/CUSUM/air_ozone/minnesota.csv",header=1)
wisconsin<-read.csv(file="/Users/Allen/Documents/Research/CODE/R/CUSUM/air_ozone/wisconsin.csv",header=1)
michigan<-read.csv(file="/Users/Allen/Documents/Research/CODE/R/CUSUM/air_ozone/michigan.csv",header=1)
indiana<-read.csv(file="/Users/Allen/Documents/Research/CODE/R/CUSUM/air_ozone/indiana.csv",header=1)
data0301<-rbind(illinois[illinois$Date=="03/01/2017",],
                iowa[iowa$Date=="03/01/2017",],
                wisconsin[wisconsin$Date=="03/01/2017",],
                michigan[michigan$Date=="03/01/2017",],
                minnesota[minnesota$Date=="03/01/2017",],
                indiana[indiana$Date=="03/01/2017",]
                )


library(ggplot2)
library(maps)
#load us map data
all_states <- map_data("state")
states <- subset(all_states, region %in% c("iowa","indiana","minnesota","wisconsin","illinois","michigan"))
p <- ggplot()
p <- p + geom_polygon( data=states, aes(x=long, y=lat, group = group),colour="black", fill="white" )
p <- p + geom_point(data=data0301,aes(x=SITE_LONGITUDE,y=SITE_LATITUDE,colour=Daily.Max.8.hour.Ozone.Concentration))+theme(legend.position='none')
p

all_counties <- map_data("county")
counties <- subset(all_counties, region %in% c("iowa","indiana","minnesota","wisconsin","illinois","michigan"))
p <- ggplot()
p <- p + geom_polygon( data=counties, aes(x=long, y=lat, group = group),colour="black", fill="white" )
p <- p + geom_point(data=data0301,aes(x=SITE_LONGITUDE,y=SITE_LATITUDE,colour=Daily.Max.8.hour.Ozone.Concentration))+theme(legend.position='none')
p