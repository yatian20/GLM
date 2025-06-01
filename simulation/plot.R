library(ggplot2)

##Case A
#Poisson
load("P402A.rda")
load("P412A.rda")
load("P422A.rda")
load("P442A.rda")
load("P482A.rda")
PA_n400 <- cbind(P402A[,3:4],P412A[,3:4],P422A[,3:4],P442A[,3:4],P482A[,3:4])
PA_n400 <- data.frame(median=apply(PA_n400,2,median),Low=apply(PA_n400,2,function(x)quantile(x,probs=0.05)),Up=apply(PA_n400,2,function(x)quantile(x,probs=0.95)),method=rep(c('SS-Hunting','SS-Refinement'),5))
PA_n400$T <- as.factor(sort(rep(c(5,10,20,40,80),2)))
fig <- ggplot(data=PA_n400,aes(x=T,y=median,group=method,color=method,linetype=method)) + geom_line(size=0.85) + geom_point(size=1)
fig <- fig + geom_errorbar(aes(ymin = Low,ymax = Up),linetype = 1,show.legend = FALSE,width = 0.1) + scale_linetype_manual(values=1:2)
fig <- fig + scale_y_continuous(trans = 'log',breaks = c(0.05,0.2,0.8),limits = c(0.03,1.5)) + theme(text = element_text(size = 17))
fig <- fig + theme(legend.position = c(0.67,0.935)) + theme(legend.title = element_blank()) + scale_color_manual(values=c("#F8766D","#7CAE00"))
fig <- fig + labs(y=expression(paste('||',widehat(ZZ^T) - Z^'*',Z^'*T','||'[F]^2,'/',n))) 
fig <- fig + theme(panel.grid.minor = element_blank()) + theme(axis.title = element_text(size=16)) + theme(legend.key.width = unit(2.3,"line"))
fig <- fig + theme(legend.margin=margin(rep(2.5,4))) + theme(plot.margin = margin(t = 1.5,r = 1.6,b = 0,l = 1.5))

pdf(file='Simu_PA.pdf',width=4,height=3.8)
fig
dev.off()

#Bernoulli
load("B402A.rda")
load("B412A.rda")
load("B422A.rda")
load("B442A.rda")
load("B482A.rda")
BA_n400 <- cbind(B402A[,3:6],B412A[,3:6],B422A[,3:6],B442A[,3:6],B482A[,3:6])
BA_n400 <- data.frame(median=apply(BA_n400,2,median),Low=apply(BA_n400,2,function(x)quantile(x,probs=0.05)),Up=apply(BA_n400,2,function(x)quantile(x,probs=0.95)),method=rep(c('SS-Hunting','SS-Refinement','MultiNeSS','MultiNeSS+'),5))
BA_n400$T <- as.factor(sort(rep(c(5,10,20,40,80),4)))
BA_n400$method <- factor(BA_n400$method,levels=c('SS-Hunting','SS-Refinement','MultiNeSS','MultiNeSS+'))
fig <- ggplot(data=BA_n400,aes(x=T,y=median,group=method,color=method,linetype=method)) + geom_line(size=0.85) + geom_point(size=1)
fig <- fig + geom_errorbar(aes(ymin = Low,ymax = Up),linetype = 1,show.legend = FALSE,width = 0.1) + scale_linetype_manual(values=1:4)
fig <- fig + scale_y_continuous(trans = 'log',breaks = c(0.25,1,4,16),limits = c(0.22,40)) + theme(text = element_text(size = 17))
fig <- fig + theme(legend.position = c(0.67,0.86)) + theme(legend.title = element_blank())
fig <- fig + labs(y=expression(paste('||',widehat(ZZ^T) - Z^'*',Z^'*T','||'[F]^2,'/',n))) 
fig <- fig + theme(panel.grid.minor = element_blank()) + theme(axis.title = element_text(size=16)) + theme(legend.key.width = unit(2.3,"line"))
fig <- fig + theme(legend.margin=margin(rep(2.5,4))) + theme(plot.margin = margin(t = 1.5,r = 1.6,b = 0,l = 1.5))

pdf(file='Simu_BA.pdf',width=4,height=3.8)
fig
dev.off()

#Gaussian
load("G402A.rda")
load("G412A.rda")
load("G422A.rda")
load("G442A.rda")
load("G482A.rda")
GA_n400 <- cbind(G402A[,3:6],G412A[,3:6],G422A[,3:6],G442A[,3:6],G482A[,3:6])
GA_n400 <- data.frame(median=apply(GA_n400,2,median),Low=apply(GA_n400,2,function(x)quantile(x,probs=0.05)),Up=apply(GA_n400,2,function(x)quantile(x,probs=0.95)),method=rep(c('SS-Hunting','SS-Refinement','MultiNeSS','MultiNeSS+'),5))
GA_n400$T <- as.factor(sort(rep(c(5,10,20,40,80),4)))
GA_n400$method <- factor(GA_n400$method,levels=c('SS-Hunting','SS-Refinement','MultiNeSS','MultiNeSS+'))
fig <- ggplot(data=GA_n400,aes(x=T,y=median,group=method,color=method,linetype=method)) + geom_line(size=0.85) + geom_point(size=1) 
fig <- fig + geom_errorbar(aes(ymin = Low,ymax = Up),linetype = 1,show.legend = FALSE,width = 0.1) + scale_linetype_manual(values=1:4)
fig <- fig + scale_y_continuous(trans = 'log',breaks = c(0.05,0.2,0.8,3.2),limits = c(0.04,6)) + theme(text = element_text(size = 17))
fig <- fig + theme(legend.position = c(0.67,0.86)) + theme(legend.title = element_blank()) 
fig <- fig + labs(y=expression(paste('||',widehat(ZZ^T) - Z^'*',Z^'*T','||'[F]^2,'/',n))) 
fig <- fig + theme(panel.grid.minor = element_blank()) + theme(axis.title = element_text(size=16)) + theme(legend.key.width = unit(2.3,"line"))
fig <- fig + theme(legend.margin=margin(rep(2.5,4))) + theme(plot.margin = margin(t = 1.5,r = 1.6,b = 0,l = 1.5))

pdf(file='Simu_GA.pdf',width=4,height=3.8)
fig
dev.off()

##Case B
#Poisson
load("P402B.rda")
load("P412B.rda")
load("P422B.rda")
load("P442B.rda")
load("P482B.rda")
PB_n400 <- cbind(P402B[,3:4],P412B[,3:4],P422B[,3:4],P442B[,3:4],P482B[,3:4])
PB_n400 <- data.frame(median=apply(PB_n400,2,median),Low=apply(PB_n400,2,function(x)quantile(x,probs=0.05)),Up=apply(PB_n400,2,function(x)quantile(x,probs=0.95)),method=rep(c('SS-Hunting','SS-Refinement'),5))
PB_n400$T <- as.factor(sort(rep(c(5,10,20,40,80),2)))
fig <- ggplot(data=PB_n400,aes(x=T,y=median,group=method,color=method,linetype=method)) + geom_line(size=0.85) + geom_point(size=1)
fig <- fig + geom_errorbar(aes(ymin = Low,ymax = Up),linetype = 1,show.legend = FALSE,width = 0.1) + scale_linetype_manual(values=1:2)
fig <- fig + scale_y_continuous(trans = 'log',breaks = c(0.05,0.2,0.8),limits = c(0.03,1.5)) + theme(text = element_text(size = 17))
fig <- fig + theme(legend.position = c(0.67,0.935)) + theme(legend.title = element_blank()) + scale_color_manual(values=c("#F8766D","#7CAE00"))
fig <- fig + labs(y=expression(paste('||',widehat(ZZ^T) - Z^'*',Z^'*T','||'[F]^2,'/',n))) 
fig <- fig + theme(panel.grid.minor = element_blank()) + theme(axis.title = element_text(size=16)) + theme(legend.key.width = unit(2.3,"line"))
fig <- fig + theme(legend.margin=margin(rep(2.5,4))) + theme(plot.margin = margin(t = 1.5,r = 1.6,b = 0,l = 1.5))

pdf(file='Simu_PB.pdf',width=4,height=3.8)
fig
dev.off()

#Bernoulli
load("B402B.rda")
load("B412B.rda")
load("B422B.rda")
load("B442B.rda")
load("B482B.rda")
BB_n400 <- cbind(B402B[,3:6],B412B[,3:6],B422B[,3:6],B442B[,3:6],B482B[,3:6])
BB_n400 <- data.frame(median=apply(BB_n400,2,median),Low=apply(BB_n400,2,function(x)quantile(x,probs=0.05)),Up=apply(BB_n400,2,function(x)quantile(x,probs=0.95)),method=rep(c('SS-Hunting','SS-Refinement','MultiNeSS','MultiNeSS+'),5))
BB_n400$T <- as.factor(sort(rep(c(5,10,20,40,80),4)))
BB_n400$method <- factor(BB_n400$method,levels=c('SS-Hunting','SS-Refinement','MultiNeSS','MultiNeSS+'))
fig <- ggplot(data=BB_n400,aes(x=T,y=median,group=method,color=method,linetype=method)) + geom_line(size=0.85) + geom_point(size=1)
fig <- fig + geom_errorbar(aes(ymin = Low,ymax = Up),linetype = 1,show.legend = FALSE,width = 0.1) + scale_linetype_manual(values=1:4)
fig <- fig + scale_y_continuous(trans = 'log',breaks = c(0.25,1,4,16),limits = c(0.23,38)) + theme(text = element_text(size = 17))
fig <- fig + theme(legend.position = c(0.67,0.86)) + theme(legend.title = element_blank()) 
fig <- fig + labs(y=expression(paste('||',widehat(ZZ^T) - Z^'*',Z^'*T','||'[F]^2,'/',n))) 
fig <- fig + theme(panel.grid.minor = element_blank()) + theme(axis.title = element_text(size=16)) + theme(legend.key.width = unit(2.3,"line"))
fig <- fig + theme(legend.margin=margin(rep(2.5,4))) + theme(plot.margin = margin(t = 1.5,r = 1.6,b = 0,l = 1.5))

pdf(file='Simu_BB.pdf',width=4,height=3.8)
fig
dev.off()

#Gaussian
load("G402B.rda")
load("G412B.rda")
load("G422B.rda")
load("G442B.rda")
load("G482B.rda")
GB_n400 <- cbind(G402B[,3:6],G412B[,3:6],G422B[,3:6],G442B[,3:6],G482B[,3:6])
GB_n400 <- data.frame(median=apply(GB_n400,2,median),Low=apply(GB_n400,2,function(x)quantile(x,probs=0.05)),Up=apply(GB_n400,2,function(x)quantile(x,probs=0.95)),method=rep(c('SS-Hunting','SS-Refinement','MultiNeSS','MultiNeSS+'),5))
GB_n400$T <- as.factor(sort(rep(c(5,10,20,40,80),4)))
GB_n400$method <- factor(GB_n400$method,levels=c('SS-Hunting','SS-Refinement','MultiNeSS','MultiNeSS+'))
fig <- ggplot(data=GB_n400,aes(x=T,y=median,group=method,color=method,linetype=method)) + geom_line(size=0.85) + geom_point(size=1) 
fig <- fig + geom_errorbar(aes(ymin = Low,ymax = Up),linetype = 1,show.legend = FALSE,width = 0.1) + scale_linetype_manual(values=1:4)
fig <- fig + scale_y_continuous(trans = 'log',breaks = c(0.05,0.2,0.8,3.2),limits = c(0.045,7)) + theme(text = element_text(size = 17))
fig <- fig + theme(legend.position = c(0.67,0.86)) + theme(legend.title = element_blank()) 
fig <- fig + labs(y=expression(paste('||',widehat(ZZ^T) - Z^'*',Z^'*T','||'[F]^2,'/',n))) 
fig <- fig + theme(panel.grid.minor = element_blank()) + theme(axis.title = element_text(size=16)) + theme(legend.key.width = unit(2.3,"line"))
fig <- fig + theme(legend.margin=margin(rep(2.5,4))) + theme(plot.margin = margin(t = 1.5,r = 1.6,b = 0,l = 1.5))

pdf(file='Simu_GB.pdf',width=4,height=3.8)
fig
dev.off()

##Case C
#Poisson
load("P402C.rda")
load("P412C.rda")
load("P422C.rda")
load("P442C.rda")
load("P482C.rda")
PC_n400 <- cbind(P402C[,3:4],P412C[,3:4],P422C[,3:4],P442C[,3:4],P482C[,3:4])
PC_n400 <- data.frame(median=apply(PC_n400,2,median),Low=apply(PC_n400,2,function(x)quantile(x,probs=0.05)),Up=apply(PC_n400,2,function(x)quantile(x,probs=0.95)),method=rep(c('SS-Hunting','SS-Refinement'),5))
PC_n400$T <- as.factor(sort(rep(c(5,10,20,40,80),2)))
fig <- ggplot(data=PC_n400,aes(x=T,y=median,group=method,color=method,linetype=method)) + geom_line(size=0.85) + geom_point(size=1)
fig <- fig + geom_errorbar(aes(ymin = Low,ymax = Up),linetype = 1,show.legend = FALSE,width = 0.1) + scale_linetype_manual(values=1:2)
fig <- fig + scale_y_continuous(trans = 'log',breaks = c(0.05,0.2,0.8),limits = c(0.03,1.5)) + theme(text = element_text(size = 17))
fig <- fig + theme(legend.position = c(0.67,0.935)) + theme(legend.title = element_blank()) + scale_color_manual(values=c("#F8766D","#7CAE00"))
fig <- fig + labs(y=expression(paste('||',widehat(ZZ^T) - Z^'*',Z^'*T','||'[F]^2,'/',n))) 
fig <- fig + theme(panel.grid.minor = element_blank()) + theme(axis.title = element_text(size=16)) + theme(legend.key.width = unit(2.3,"line"))
fig <- fig + theme(legend.margin=margin(rep(2.5,4))) + theme(plot.margin = margin(t = 1.5,r = 1.6,b = 0,l = 1.5))

pdf(file='Simu_PC.pdf',width=4,height=3.8)
fig
dev.off()

#Bernoulli
load("B402C.rda")
load("B412C.rda")
load("B422C.rda")
load("B442C.rda")
load("B482C.rda")
BC_n400 <- cbind(B402C[,3:4],B412C[,3:4],B422C[,3:4],B442C[,3:4],B482C[,3:4])
BC_n400 <- data.frame(median=apply(BC_n400,2,median),Low=apply(BC_n400,2,function(x)quantile(x,probs=0.05)),Up=apply(BC_n400,2,function(x)quantile(x,probs=0.95)),method=rep(c('SS-Hunting','SS-Refinement'),5))
BC_n400$T <- as.factor(sort(rep(c(5,10,20,40,80),2)))
fig <- ggplot(data=BC_n400,aes(x=T,y=median,group=method,color=method,linetype=method)) + geom_line(size=0.85) + geom_point(size=1)
fig <- fig + geom_errorbar(aes(ymin = Low,ymax = Up),linetype = 1,show.legend = FALSE,width = 0.1) + scale_linetype_manual(values=1:2)
fig <- fig + scale_y_continuous(trans = 'log',breaks = c(0.25,1,4,16),limits = c(0.2,10)) + theme(text = element_text(size = 17))
fig <- fig + theme(legend.position = c(0.67,0.935)) + theme(legend.title = element_blank()) + scale_color_manual(values=c("#F8766D","#7CAE00"))
fig <- fig + labs(y=expression(paste('||',widehat(ZZ^T) - Z^'*',Z^'*T','||'[F]^2,'/',n))) 
fig <- fig + theme(panel.grid.minor = element_blank()) + theme(axis.title = element_text(size=16)) + theme(legend.key.width = unit(2.3,"line"))
fig <- fig + theme(legend.margin=margin(rep(2.5,4))) + theme(plot.margin = margin(t = 1.5,r = 1.6,b = 0,l = 1.5))

pdf(file='Simu_BC.pdf',width=4,height=3.8)
fig
dev.off()

#Gaussian
load("G402C.rda")
load("G412C.rda")
load("G422C.rda")
load("G442C.rda")
load("G482C.rda")
GC_n400 <- cbind(G402C[,3:4],G412C[,3:4],G422C[,3:4],G442C[,3:4],G482C[,3:4])
GC_n400 <- data.frame(median=apply(GC_n400,2,median),Low=apply(GC_n400,2,function(x)quantile(x,probs=0.05)),Up=apply(GC_n400,2,function(x)quantile(x,probs=0.95)),method=rep(c('SS-Hunting','SS-Refinement'),5))
GC_n400$T <- as.factor(sort(rep(c(5,10,20,40,80),2)))
fig <- ggplot(data=GC_n400,aes(x=T,y=median,group=method,color=method,linetype=method)) + geom_line(size=0.85) + geom_point(size=1) 
fig <- fig + geom_errorbar(aes(ymin = Low,ymax = Up),linetype = 1,show.legend = FALSE,width = 0.1) + scale_linetype_manual(values=1:2)
fig <- fig + scale_y_continuous(trans = 'log',breaks = c(0.05,0.2,0.8,3.2),limits = c(0.04,2)) + theme(text = element_text(size = 17))
fig <- fig + theme(legend.position = c(0.67,0.935)) + theme(legend.title = element_blank()) + scale_color_manual(values=c("#F8766D","#7CAE00"))
fig <- fig + labs(y=expression(paste('||',widehat(ZZ^T) - Z^'*',Z^'*T','||'[F]^2,'/',n))) 
fig <- fig + theme(panel.grid.minor = element_blank()) + theme(axis.title = element_text(size=16)) + theme(legend.key.width = unit(2.3,"line"))
fig <- fig + theme(legend.margin=margin(rep(2.5,4))) + theme(plot.margin = margin(t = 1.5,r = 1.6,b = 0,l = 1.5))

pdf(file='Simu_GC.pdf',width=4,height=3.8)
fig
dev.off()

#### error of W1 ####

##Case A
#Poisson
load("P402A.rda")
load("P412A.rda")
load("P422A.rda")
load("P442A.rda")
load("P482A.rda")
PA_n400 <- cbind(P402A[,5:6],P412A[,5:6],P422A[,5:6],P442A[,5:6],P482A[,5:6])
PA_n400 <- data.frame(median=apply(PA_n400,2,median),Low=apply(PA_n400,2,function(x)quantile(x,probs=0.05)),Up=apply(PA_n400,2,function(x)quantile(x,probs=0.95)),method=rep(c('SS-Hunting','SS-Refinement'),5))
PA_n400$T <- as.factor(sort(rep(c(5,10,20,40,80),2)))
fig <- ggplot(data=PA_n400,aes(x=T,y=median,group=method,color=method,linetype=method)) + geom_line(size=0.85) + geom_point(size=1)
fig <- fig + geom_errorbar(aes(ymin = Low,ymax = Up),linetype = 1,show.legend = FALSE,width = 0.1) + scale_linetype_manual(values=1:2)
fig <- fig + scale_y_continuous(trans = 'log',breaks = c(1,2,4,8),limits = c(1.6,8)) + theme(text = element_text(size = 17))
fig <- fig + theme(legend.position = c(0.7,0.935)) + theme(legend.title = element_blank()) + scale_color_manual(values=c("#F8766D","#7CAE00"))
fig <- fig + labs(y=expression(paste('||',widehat(paste(W[1],W[1]^T)) - W[1]^'*',W[1]^'*T','||'[F]^2,'/',n))) 
fig <- fig + theme(panel.grid.minor = element_blank()) + theme(axis.title = element_text(size=16)) + theme(legend.key.width = unit(2.3,"line"))
fig <- fig + theme(legend.margin=margin(rep(2.5,4))) + theme(plot.margin = margin(t = 1.5,r = 1.6,b = 0,l = 1.5))

pdf(file='SimuW_PA.pdf',width=4,height=3.8)
fig
dev.off()

#Bernoulli
load("B402A.rda")
load("B412A.rda")
load("B422A.rda")
load("B442A.rda")
load("B482A.rda")
BA_n400 <- cbind(B402A[,7:10],B412A[,7:10],B422A[,7:10],B442A[,7:10],B482A[,7:10])
BA_n400 <- data.frame(median=apply(BA_n400,2,median),Low=apply(BA_n400,2,function(x)quantile(x,probs=0.05)),Up=apply(BA_n400,2,function(x)quantile(x,probs=0.95)),method=rep(c('SS-Hunting','SS-Refinement','MultiNeSS','MultiNeSS+'),5))
BA_n400$T <- as.factor(sort(rep(c(5,10,20,40,80),4)))
BA_n400$method <- factor(BA_n400$method,levels=c('SS-Hunting','SS-Refinement','MultiNeSS','MultiNeSS+'))
fig <- ggplot(data=BA_n400,aes(x=T,y=median,group=method,color=method,linetype=method)) + geom_line(size=0.85) + geom_point(size=1)
fig <- fig + geom_errorbar(aes(ymin = Low,ymax = Up),linetype = 1,show.legend = FALSE,width = 0.1) + scale_linetype_manual(values=1:4)
fig <- fig + scale_y_continuous(trans = 'log',breaks = c(8,16,32,64),limits = c(15,110)) + theme(text = element_text(size = 17))
fig <- fig + theme(legend.position = c(0.7,0.86)) + theme(legend.title = element_blank()) 
fig <- fig + labs(y=expression(paste('||',widehat(paste(W[1],W[1]^T)) - W[1]^'*',W[1]^'*T','||'[F]^2,'/',n))) 
fig <- fig + theme(panel.grid.minor = element_blank()) + theme(axis.title = element_text(size=16)) + theme(legend.key.width = unit(2.3,"line"))
fig <- fig + theme(legend.margin=margin(rep(2.5,4))) + theme(plot.margin = margin(t = 1.5,r = 1.6,b = 0,l = 1.5))

pdf(file='SimuW_BA.pdf',width=4,height=3.8)
fig
dev.off()

#Gaussian
load("G402A.rda")
load("G412A.rda")
load("G422A.rda")
load("G442A.rda")
load("G482A.rda")
GA_n400 <- cbind(G402A[,7:10],G412A[,7:10],G422A[,7:10],G442A[,7:10],G482A[,7:10])
GA_n400 <- data.frame(median=apply(GA_n400,2,median),Low=apply(GA_n400,2,function(x)quantile(x,probs=0.05)),Up=apply(GA_n400,2,function(x)quantile(x,probs=0.95)),method=rep(c('SS-Hunting','SS-Refinement','MultiNeSS','MultiNeSS+'),5))
GA_n400$T <- as.factor(sort(rep(c(5,10,20,40,80),4)))
GA_n400$method <- factor(GA_n400$method,levels=c('SS-Hunting','SS-Refinement','MultiNeSS','MultiNeSS+'))
fig <- ggplot(data=GA_n400,aes(x=T,y=median,group=method,color=method,linetype=method)) + geom_line(size=0.85) + geom_point(size=1) 
fig <- fig + geom_errorbar(aes(ymin = Low,ymax = Up),linetype = 1,show.legend = FALSE,width = 0.1) + scale_linetype_manual(values=1:4)
fig <- fig + scale_y_continuous(trans = 'log',breaks = c(2,4,8,16,32),limits = c(3,40)) + theme(text = element_text(size = 17))
fig <- fig + theme(legend.position = c(0.7,0.86)) + theme(legend.title = element_blank()) 
fig <- fig + labs(y=expression(paste('||',widehat(paste(W[1],W[1]^T)) - W[1]^'*',W[1]^'*T','||'[F]^2,'/',n))) 
fig <- fig + theme(panel.grid.minor = element_blank()) + theme(axis.title = element_text(size=16)) + theme(legend.key.width = unit(2.3,"line"))
fig <- fig + theme(legend.margin=margin(rep(2.5,4))) + theme(plot.margin = margin(t = 1.5,r = 1.6,b = 0,l = 1.5))

pdf(file='SimuW_GA.pdf',width=4,height=3.8)
fig
dev.off()

##Case B
#Poisson
load("P402B.rda")
load("P412B.rda")
load("P422B.rda")
load("P442B.rda")
load("P482B.rda")
PB_n400 <- cbind(P402B[,5:6],P412B[,5:6],P422B[,5:6],P442B[,5:6],P482B[,5:6])
PB_n400 <- data.frame(median=apply(PB_n400,2,median),Low=apply(PB_n400,2,function(x)quantile(x,probs=0.05)),Up=apply(PB_n400,2,function(x)quantile(x,probs=0.95)),method=rep(c('SS-Hunting','SS-Refinement'),5))
PB_n400$T <- as.factor(sort(rep(c(5,10,20,40,80),2)))
fig <- ggplot(data=PB_n400,aes(x=T,y=median,group=method,color=method,linetype=method)) + geom_line(size=0.85) + geom_point(size=1)
fig <- fig + geom_errorbar(aes(ymin = Low,ymax = Up),linetype = 1,show.legend = FALSE,width = 0.1) + scale_linetype_manual(values=1:2)
fig <- fig + scale_y_continuous(trans = 'log',breaks = c(1,2,4,8),limits = c(1.6,8)) + theme(text = element_text(size = 17))
fig <- fig + theme(legend.position = c(0.7,0.935)) + theme(legend.title = element_blank()) + scale_color_manual(values=c("#F8766D","#7CAE00"))
fig <- fig + labs(y=expression(paste('||',widehat(paste(W[1],W[1]^T)) - W[1]^'*',W[1]^'*T','||'[F]^2,'/',n))) 
fig <- fig + theme(panel.grid.minor = element_blank()) + theme(axis.title = element_text(size=16)) + theme(legend.key.width = unit(2.3,"line"))
fig <- fig + theme(legend.margin=margin(rep(2.5,4))) + theme(plot.margin = margin(t = 1.5,r = 1.6,b = 0,l = 1.5))

pdf(file='SimuW_PB.pdf',width=4,height=3.8)
fig
dev.off()

#Bernoulli
load("B402B.rda")
load("B412B.rda")
load("B422B.rda")
load("B442B.rda")
load("B482B.rda")
BB_n400 <- cbind(B402B[,7:10],B412B[,7:10],B422B[,7:10],B442B[,7:10],B482B[,7:10])
BB_n400 <- data.frame(median=apply(BB_n400,2,median),Low=apply(BB_n400,2,function(x)quantile(x,probs=0.05)),Up=apply(BB_n400,2,function(x)quantile(x,probs=0.95)),method=rep(c('SS-Hunting','SS-Refinement','MultiNeSS','MultiNeSS+'),5))
BB_n400$T <- as.factor(sort(rep(c(5,10,20,40,80),4)))
BB_n400$method <- factor(BB_n400$method,levels=c('SS-Hunting','SS-Refinement','MultiNeSS','MultiNeSS+'))
fig <- ggplot(data=BB_n400,aes(x=T,y=median,group=method,color=method,linetype=method)) + geom_line(size=0.85) + geom_point(size=1)
fig <- fig + geom_errorbar(aes(ymin = Low,ymax = Up),linetype = 1,show.legend = FALSE,width = 0.1) + scale_linetype_manual(values=1:4)
fig <- fig + scale_y_continuous(trans = 'log',breaks = c(8,16,32,64),limits = c(15,110)) + theme(text = element_text(size = 17))
fig <- fig + theme(legend.position = c(0.7,0.86)) + theme(legend.title = element_blank()) 
fig <- fig + labs(y=expression(paste('||',widehat(paste(W[1],W[1]^T)) - W[1]^'*',W[1]^'*T','||'[F]^2,'/',n))) 
fig <- fig + theme(panel.grid.minor = element_blank()) + theme(axis.title = element_text(size=16)) + theme(legend.key.width = unit(2.3,"line"))
fig <- fig + theme(legend.margin=margin(rep(2.5,4))) + theme(plot.margin = margin(t = 1.5,r = 1.6,b = 0,l = 1.5))

pdf(file='SimuW_BB.pdf',width=4,height=3.8)
fig
dev.off()

#Gaussian
load("G402B.rda")
load("G412B.rda")
load("G422B.rda")
load("G442B.rda")
load("G482B.rda")
GB_n400 <- cbind(G402B[,7:10],G412B[,7:10],G422B[,7:10],G442B[,7:10],G482B[,7:10])
GB_n400 <- data.frame(median=apply(GB_n400,2,median),Low=apply(GB_n400,2,function(x)quantile(x,probs=0.05)),Up=apply(GB_n400,2,function(x)quantile(x,probs=0.95)),method=rep(c('SS-Hunting','SS-Refinement','MultiNeSS','MultiNeSS+'),5))
GB_n400$T <- as.factor(sort(rep(c(5,10,20,40,80),4)))
GB_n400$method <- factor(GB_n400$method,levels=c('SS-Hunting','SS-Refinement','MultiNeSS','MultiNeSS+'))
fig <- ggplot(data=GB_n400,aes(x=T,y=median,group=method,color=method,linetype=method)) + geom_line(size=0.85) + geom_point(size=1) 
fig <- fig + geom_errorbar(aes(ymin = Low,ymax = Up),linetype = 1,show.legend = FALSE,width = 0.1) + scale_linetype_manual(values=1:4)
fig <- fig + scale_y_continuous(trans = 'log',breaks = c(2,4,8,16,32),limits = c(3,40)) + theme(text = element_text(size = 17))
fig <- fig + theme(legend.position = c(0.7,0.86)) + theme(legend.title = element_blank()) 
fig <- fig + labs(y=expression(paste('||',widehat(paste(W[1],W[1]^T)) - W[1]^'*',W[1]^'*T','||'[F]^2,'/',n))) 
fig <- fig + theme(panel.grid.minor = element_blank()) + theme(axis.title = element_text(size=16)) + theme(legend.key.width = unit(2.3,"line"))
fig <- fig + theme(legend.margin=margin(rep(2.5,4))) + theme(plot.margin = margin(t = 1.5,r = 1.6,b = 0,l = 1.5))

pdf(file='SimuW_GB.pdf',width=4,height=3.8)
fig
dev.off()

##Case C
#Poisson
load("P402C.rda")
load("P412C.rda")
load("P422C.rda")
load("P442C.rda")
load("P482C.rda")
PC_n400 <- cbind(P402C[,5:6],P412C[,5:6],P422C[,5:6],P442C[,5:6],P482C[,5:6])
PC_n400 <- data.frame(median=apply(PC_n400,2,median),Low=apply(PC_n400,2,function(x)quantile(x,probs=0.05)),Up=apply(PC_n400,2,function(x)quantile(x,probs=0.95)),method=rep(c('SS-Hunting','SS-Refinement'),5))
PC_n400$T <- as.factor(sort(rep(c(5,10,20,40,80),2)))
fig <- ggplot(data=PC_n400,aes(x=T,y=median,group=method,color=method,linetype=method)) + geom_line(size=0.85) + geom_point(size=1)
fig <- fig + geom_errorbar(aes(ymin = Low,ymax = Up),linetype = 1,show.legend = FALSE,width = 0.1) + scale_linetype_manual(values=1:2)
fig <- fig + scale_y_continuous(trans = 'log',breaks = c(1,2,4,8),limits = c(1.6,8)) + theme(text = element_text(size = 17))
fig <- fig + theme(legend.position = c(0.7,0.935)) + theme(legend.title = element_blank()) + scale_color_manual(values=c("#F8766D","#7CAE00"))
fig <- fig + labs(y=expression(paste('||',widehat(paste(W[1],W[1]^T)) - W[1]^'*',W[1]^'*T','||'[F]^2,'/',n))) 
fig <- fig + theme(panel.grid.minor = element_blank()) + theme(axis.title = element_text(size=16)) + theme(legend.key.width = unit(2.3,"line"))
fig <- fig + theme(legend.margin=margin(rep(2.5,4))) + theme(plot.margin = margin(t = 1.5,r = 1.6,b = 0,l = 1.5))

pdf(file='SimuW_PC.pdf',width=4,height=3.8)
fig
dev.off()

#Bernoulli
load("B402C.rda")
load("B412C.rda")
load("B422C.rda")
load("B442C.rda")
load("B482C.rda")
BC_n400 <- cbind(B402C[,7:8],B412C[,7:8],B422C[,7:8],B442C[,7:8],B482C[,7:8])
BC_n400 <- data.frame(median=apply(BC_n400,2,median),Low=apply(BC_n400,2,function(x)quantile(x,probs=0.05)),Up=apply(BC_n400,2,function(x)quantile(x,probs=0.95)),method=rep(c('SS-Hunting','SS-Refinement'),5))
BC_n400$T <- as.factor(sort(rep(c(5,10,20,40,80),2)))
fig <- ggplot(data=BC_n400,aes(x=T,y=median,group=method,color=method,linetype=method)) + geom_line(size=0.85) + geom_point(size=1)
fig <- fig + geom_errorbar(aes(ymin = Low,ymax = Up),linetype = 1,show.legend = FALSE,width = 0.1) + scale_linetype_manual(values=1:2)
fig <- fig + scale_y_continuous(trans = 'log',breaks = c(8,16,32,64),limits = c(12,70)) + theme(text = element_text(size = 17))
fig <- fig + theme(legend.position = c(0.7,0.935)) + theme(legend.title = element_blank()) + scale_color_manual(values=c("#F8766D","#7CAE00"))
fig <- fig + labs(y=expression(paste('||',widehat(paste(W[1],W[1]^T)) - W[1]^'*',W[1]^'*T','||'[F]^2,'/',n))) 
fig <- fig + theme(panel.grid.minor = element_blank()) + theme(axis.title = element_text(size=16)) + theme(legend.key.width = unit(2.3,"line"))
fig <- fig + theme(legend.margin=margin(rep(2.5,4))) + theme(plot.margin = margin(t = 1.5,r = 1.6,b = 0,l = 1.5))

pdf(file='SimuW_BC.pdf',width=4,height=3.8)
fig
dev.off()

#Gaussian
load("G402C.rda")
load("G412C.rda")
load("G422C.rda")
load("G442C.rda")
load("G482C.rda")
GC_n400 <- cbind(G402C[,7:8],G412C[,7:8],G422C[,7:8],G442C[,7:8],G482C[,7:8])
GC_n400 <- data.frame(median=apply(GC_n400,2,median),Low=apply(GC_n400,2,function(x)quantile(x,probs=0.05)),Up=apply(GC_n400,2,function(x)quantile(x,probs=0.95)),method=rep(c('SS-Hunting','SS-Refinement'),5))
GC_n400$T <- as.factor(sort(rep(c(5,10,20,40,80),2)))
fig <- ggplot(data=GC_n400,aes(x=T,y=median,group=method,color=method,linetype=method)) + geom_line(size=0.85) + geom_point(size=1) 
fig <- fig + geom_errorbar(aes(ymin = Low,ymax = Up),linetype = 1,show.legend = FALSE,width = 0.1) + scale_linetype_manual(values=1:2)
fig <- fig + scale_y_continuous(trans = 'log',breaks = c(2,4,8,16,32),limits = c(2,12)) + theme(text = element_text(size = 17))
fig <- fig + theme(legend.position = c(0.7,0.935)) + theme(legend.title = element_blank()) + scale_color_manual(values=c("#F8766D","#7CAE00"))
fig <- fig + labs(y=expression(paste('||',widehat(paste(W[1],W[1]^T)) - W[1]^'*',W[1]^'*T','||'[F]^2,'/',n))) 
fig <- fig + theme(panel.grid.minor = element_blank()) + theme(axis.title = element_text(size=16)) + theme(legend.key.width = unit(2.3,"line"))
fig <- fig + theme(legend.margin=margin(rep(2.5,4))) + theme(plot.margin = margin(t = 1.5,r = 1.6,b = 0,l = 1.5))

pdf(file='SimuW_GC.pdf',width=4,height=3.8)
fig
dev.off()
