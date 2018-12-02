library(ggplot2)
library(grid)

N.grid=100
shift=0
Status<-rep(0,(N.grid)^2)
grid.point = expand.grid(x=seq(1,N.grid,by=1),y=seq(1,N.grid,by=1))
grid.point$Status<-Status
grid.point[which(grid.point$x<=0.75*N.grid&grid.point$x>=.25*N.grid&grid.point$y<=.75*N.grid&grid.point$y>=.25*N.grid),3]<-1
grid.point$Observed<-rnorm((N.grid)^2)
grid.point[grid.point$Status==1,]$Observed<-grid.point[grid.point$Status==1,]$Observed+shift

library(ggplot2)
a<-ggplot(grid.point) + 
  geom_tile(aes(x = x, y = y,fill=Observed))+
  xlab("")+ylab("")+ 
  scale_fill_distiller(palette = "RdBu")+
  theme(legend.position = 'none',panel.grid.major =element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(),axis.text = element_blank(),axis.ticks = element_blank())

k=5
times=k^2
grid.point$detection_prob<-signal_detection(grid.point,k,times,100)
x_seq<-seq(0,1,length.out = times*k^2)
est <- lpdensity(data = grid.point$detection_prob,grid=x_seq, bwselect = "imse-rot")$Estimate[,4]
est[est<0]=0
bb<-ggplot(grid.point)+geom_histogram(aes(detection_prob,fill=as.factor(Status)),binwidth=0.05,colour="black", alpha=1, position="identity")+
  xlab("Signal Weight")+ylab("Frequency")+theme_bw()+
  theme(legend.title=element_blank(),legend.position = "top", axis.text=element_text(size = 15,  face = "bold"),axis.title=element_text(size = 12,  face = "bold"))+ 
  scale_fill_discrete(breaks = c('1','0'), labels = c('Signal','Indifference'))

d<-data.frame(x=c(1:625)/625,y=as.vector(est))

d$z<-c(d$y[1:(which.min(d$y))],seq(min(d$y),0,length.out=625-which.min(d$y)))
d$w<-d$y-d$z

dd=data.frame(x=rep(d$x,times=3),y=c(d$y,d$z,d$w),z=rep(c("black","indianred2","mediumaquamarine"),each=625))


cc<-ggplot(dd)+ylim(c(0,4))+
  xlab("Signal Weight")+ylab("Density")+theme_bw()+
  theme(legend.title = element_blank(),legend.position = "top",axis.text=element_text(size = 12,  face = "bold"),axis.title=element_text(size = 12,  face = "bold"))
cc<-cc+geom_line(aes(x=x,y=y,colour=z))+scale_color_manual(breaks=c("black","indianred2","mediumaquamarine"),values=c("black","indianred2","mediumaquamarine"), labels = c(expression(f(0)),expression(tilde(f)[H[1]](x)),expression(tilde(f)[H[0]](x))))



grid.newpage()  ##newpage
pushViewport(viewport(layout = grid.layout(1,3))) ### divide into 1*3
vplayout <- function(x,y){
  viewport(layout.pos.row = x, layout.pos.col = y)
}
print(a, vp = vplayout(1,1))   ### put a in (1,1)
print(bb, vp = vplayout(1,2))   ### put b in (1,2)
print(cc, vp = vplayout(1,3))  ###put c in (1,3)
