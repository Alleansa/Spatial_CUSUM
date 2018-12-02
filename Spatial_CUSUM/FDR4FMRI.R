pvalue_normal<-function(Observed){
  pv<-1-pnorm(Observed)
  return(pv)
}

p_val<-pvalue_normal(grid.point$Observed)
grid.point$fdr_p<-p.adjust(p_val, method = "fdr")
grid.point$fdr_p[is.na(grid.point$fdr_p)]=1
a+geom_point(data=grid.point, aes(x = x, y = y,colour=fdr_p<=0.0001))+
  scale_colour_manual(values = c(NA,"black"))+ 
  xlab("")+ylab("")+
  theme(legend.position='none',
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank())


                    