library(ggplot2)
d<-data.frame(x=c(1:150),y=c(10+rnorm(60),sample(c(c(10+rnorm(15),c(-10+rnorm(15))))),-10+rnorm(60)))
#d<-data.frame(x=c(1:150),y=rnorm(150))
ggplot(d)+geom_point(aes(x=x,y=y))+xlab("")+ylab("")+
  geom_vline(aes(xintercept=60,colour="red"))+
  geom_vline(aes(xintercept=90,colour="red"))+
  theme(legend.position = "none",axis.text=element_text(size = 12,  face = "bold"))+
  scale_y_continuous(limits = c(-20, 20),breaks=c( -10,10), labels = c( expression(mu[0]),expression(mu[1])))+
  scale_x_continuous(limits = c(0, 150),breaks=c(30,75,120), labels = c("signal",  "interim","indifference"))

# ggplot(d)+geom_point(aes(x=x,y=y))+xlab("")+ylab("")+
#   theme(axis.text.y=element_text(size = 12,  face = "bold"),
#         axis.text.x = element_blank())+
#       scale_y_continuous(limits = c(-20, 20),breaks=c(0), labels = c( expression(mu[0])))
