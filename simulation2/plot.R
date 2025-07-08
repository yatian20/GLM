library(ggplot2)
setwd('GLM/simulation2/results')

##Figure S4
#Poisson
load("P402A.rda")
BA_n400 <- data.frame(error = as.vector(t(P402A[,2:4])))
BA_n400$method <- factor(rep(c('Simulation','Original','Theory'),100),levels=c('Original','Theory','Simulation'))
fig <- ggplot(data=BA_n400,aes(x=method,y=error,color=method)) + geom_boxplot() + theme(legend.position = "none")
fig <- fig + labs(y=expression(paste('||',widehat(ZZ^T) - Z^'*',Z^'*T','||'[F]^2,'/',n)),x=NULL) 
fig <- fig + theme(panel.grid.minor = element_blank()) + theme(axis.title = element_text(size=14),axis.text.x = element_text(size = 14,face = "bold"))

pdf(file='Simu_SPA.pdf',width=4,height=2)
fig
dev.off()

#Bernoulli
load("B402A.rda")
BA_n400 <- data.frame(error = as.vector(t(B402A[,2:4])))
BA_n400$method <- factor(rep(c('Simulation','Original','Theory'),100),levels=c('Original','Theory','Simulation'))
fig <- ggplot(data=BA_n400,aes(x=method,y=error,color=method)) + geom_boxplot() + theme(legend.position = "none")
fig <- fig + labs(y=expression(paste('||',widehat(ZZ^T) - Z^'*',Z^'*T','||'[F]^2,'/',n)),x=NULL) 
fig <- fig + theme(panel.grid.minor = element_blank()) + theme(axis.title = element_text(size=14),axis.text.x = element_text(size = 14,face = "bold"))

pdf(file='Simu_SBA.pdf',width=4,height=2)
fig
dev.off()

#Gaussian
load("G402A.rda")
BA_n400 <- data.frame(error = as.vector(t(G402A[,2:4])))
BA_n400$method <- factor(rep(c('Simulation','Original','Theory'),100),levels=c('Original','Theory','Simulation'))
fig <- ggplot(data=BA_n400,aes(x=method,y=error,color=method)) + geom_boxplot() + theme(legend.position = "none")
fig <- fig + labs(y=expression(paste('||',widehat(ZZ^T) - Z^'*',Z^'*T','||'[F]^2,'/',n)),x=NULL) 
fig <- fig + theme(panel.grid.minor = element_blank()) + theme(axis.title = element_text(size=14),axis.text.x = element_text(size = 14,face = "bold"))

pdf(file='Simu_SGA.pdf',width=4,height=2)
fig
dev.off()

##Figure S5
#Poisson
load("P402B.rda")
BA_n400 <- data.frame(error = as.vector(t(P402A[,2:4])))
BA_n400$method <- factor(rep(c('Simulation','Original','Theory'),100),levels=c('Original','Theory','Simulation'))
fig <- ggplot(data=BA_n400,aes(x=method,y=error,color=method)) + geom_boxplot() + theme(legend.position = "none")
fig <- fig + labs(y=expression(paste('||',widehat(ZZ^T) - Z^'*',Z^'*T','||'[F]^2,'/',n)),x=NULL) 
fig <- fig + theme(panel.grid.minor = element_blank()) + theme(axis.title = element_text(size=14),axis.text.x = element_text(size = 14,face = "bold"))

pdf(file='Simu_SPB.pdf',width=4,height=2)
fig
dev.off()

#Bernoulli
load("B402B.rda")
BA_n400 <- data.frame(error = as.vector(t(B402A[,2:4])))
BA_n400$method <- factor(rep(c('Simulation','Original','Theory'),100),levels=c('Original','Theory','Simulation'))
fig <- ggplot(data=BA_n400,aes(x=method,y=error,color=method)) + geom_boxplot() + theme(legend.position = "none")
fig <- fig + labs(y=expression(paste('||',widehat(ZZ^T) - Z^'*',Z^'*T','||'[F]^2,'/',n)),x=NULL) 
fig <- fig + theme(panel.grid.minor = element_blank()) + theme(axis.title = element_text(size=14),axis.text.x = element_text(size = 14,face = "bold"))

pdf(file='Simu_SBB.pdf',width=4,height=2)
fig
dev.off()

#Gaussian
load("G402B.rda")
BA_n400 <- data.frame(error = as.vector(t(G402A[,2:4])))
BA_n400$method <- factor(rep(c('Simulation','Original','Theory'),100),levels=c('Original','Theory','Simulation'))
fig <- ggplot(data=BA_n400,aes(x=method,y=error,color=method)) + geom_boxplot() + theme(legend.position = "none")
fig <- fig + labs(y=expression(paste('||',widehat(ZZ^T) - Z^'*',Z^'*T','||'[F]^2,'/',n)),x=NULL) 
fig <- fig + theme(panel.grid.minor = element_blank()) + theme(axis.title = element_text(size=14),axis.text.x = element_text(size = 14,face = "bold"))

pdf(file='Simu_SGB.pdf',width=4,height=2)
fig
dev.off()

##Figure S6
#Poisson
load("P402AC.rda")
BA_n400 <- data.frame(error = as.vector(t(P402A[,2:4])))
BA_n400$method <- factor(rep(c('Simulation','Original','Theory'),100),levels=c('Original','Theory','Simulation'))
fig <- ggplot(data=BA_n400,aes(x=method,y=error,color=method)) + geom_boxplot() + theme(legend.position = "none")
fig <- fig + labs(y=expression(paste('||',widehat(ZZ^T) - Z^'*',Z^'*T','||'[F]^2,'/',n)),x=NULL) 
fig <- fig + theme(panel.grid.minor = element_blank()) + theme(axis.title = element_text(size=14),axis.text.x = element_text(size = 14,face = "bold"))

pdf(file='Simu_SPC.pdf',width=4,height=2)
fig
dev.off()

#Bernoulli
load("B402AC.rda")
BA_n400 <- data.frame(error = as.vector(t(B402A[,2:4])))
BA_n400$method <- factor(rep(c('Simulation','Original','Theory'),100),levels=c('Original','Theory','Simulation'))
fig <- ggplot(data=BA_n400,aes(x=method,y=error,color=method)) + geom_boxplot() + theme(legend.position = "none")
fig <- fig + labs(y=expression(paste('||',widehat(ZZ^T) - Z^'*',Z^'*T','||'[F]^2,'/',n)),x=NULL) 
fig <- fig + theme(panel.grid.minor = element_blank()) + theme(axis.title = element_text(size=14),axis.text.x = element_text(size = 14,face = "bold"))

pdf(file='Simu_SBC.pdf',width=4,height=2)
fig
dev.off()

#Gaussian
load("G402AC.rda")
BA_n400 <- data.frame(error = as.vector(t(G402A[,2:4])))
BA_n400$method <- factor(rep(c('Simulation','Original','Theory'),100),levels=c('Original','Theory','Simulation'))
fig <- ggplot(data=BA_n400,aes(x=method,y=error,color=method)) + geom_boxplot() + theme(legend.position = "none")
fig <- fig + labs(y=expression(paste('||',widehat(ZZ^T) - Z^'*',Z^'*T','||'[F]^2,'/',n)),x=NULL) 
fig <- fig + theme(panel.grid.minor = element_blank()) + theme(axis.title = element_text(size=14),axis.text.x = element_text(size = 14,face = "bold"))

pdf(file='Simu_SGC.pdf',width=4,height=2)
fig
dev.off()

##Figure S7
#Bernoulli
load("B402C.rda")
load("B412C.rda")
load("B422C.rda")
load("B442C.rda")
load("B482C.rda")
BA_n400 <- cbind(B402C[,5:6],B412C[,5:6],B422C[,5:6],B442C[,5:6],B482C[,5:6])
BA_n400 <- data.frame(median=apply(BA_n400,2,median),Low=apply(BA_n400,2,function(x)quantile(x,probs=0.05)),Up=apply(BA_n400,2,function(x)quantile(x,probs=0.95)),method=rep(c('MultiNeSS','MultiNeSS+'),5))
BA_n400$T <- as.factor(sort(rep(c(5,10,20,40,80),2)))
BA_n400$method <- factor(BA_n400$method,levels=c('MultiNeSS','MultiNeSS+'))
fig <- ggplot(data=BA_n400,aes(x=T,y=median,group=method,color=method,linetype=method)) + geom_line(size=0.85) + geom_point(size=1)
fig <- fig + geom_errorbar(aes(ymin = Low,ymax = Up),linetype = 1,show.legend = FALSE,width = 0.1) + scale_linetype_manual(values=1:4)
fig <- fig + scale_y_continuous(trans = 'log',breaks = c(5,20,80),limits = c(3.5,200)) + theme(text = element_text(size = 17))
fig <- fig + theme(legend.position = c(0.75,0.94)) + theme(legend.title = element_blank()) 
fig <- fig + labs(y=expression(paste('||',widehat(ZZ^T) - Z^'*',Z^'*T','||'[F]^2,'/',n))) 
fig <- fig + theme(panel.grid.minor = element_blank()) + theme(axis.title = element_text(size=16)) + theme(legend.key.width = unit(2.3,"line"))
fig <- fig + theme(legend.margin=margin(rep(2.5,4))) + theme(plot.margin = margin(t = 1.5,r = 1.6,b = 0,l = 1.5))

pdf(file='Simu_MB1.pdf',width=4,height=3.8)
fig
dev.off()

#Gaussian                                                                                                                     
load("G402C.rda")
load("G412C.rda")
load("G422C.rda")
load("G442C.rda")
load("G482C.rda")
BA_n400 <- cbind(G402C[,5:6],G412C[,5:6],G422C[,5:6],G442C[,5:6],G482C[,5:6])
BA_n400 <- data.frame(median=apply(BA_n400,2,median),Low=apply(BA_n400,2,function(x)quantile(x,probs=0.05)),Up=apply(BA_n400,2,function(x)quantile(x,probs=0.95)),method=rep(c('MultiNeSS','MultiNeSS+'),5))
BA_n400$T <- as.factor(sort(rep(c(5,10,20,40,80),2)))
BA_n400$method <- factor(BA_n400$method,levels=c('MultiNeSS','MultiNeSS+'))
fig <- ggplot(data=BA_n400,aes(x=T,y=median,group=method,color=method,linetype=method)) + geom_line(size=0.85) + geom_point(size=1)
fig <- fig + geom_errorbar(aes(ymin = Low,ymax = Up),linetype = 1,show.legend = FALSE,width = 0.1) + scale_linetype_manual(values=1:4)
fig <- fig + scale_y_continuous(trans = 'log',breaks = c(1,10,100),limits = c(0.5,300)) + theme(text = element_text(size = 17))
fig <- fig + theme(legend.position = c(0.75,0.94)) + theme(legend.title = element_blank()) 
fig <- fig + labs(y=expression(paste('||',widehat(ZZ^T) - Z^'*',Z^'*T','||'[F]^2,'/',n))) 
fig <- fig + theme(panel.grid.minor = element_blank()) + theme(axis.title = element_text(size=16)) + theme(legend.key.width = unit(2.3,"line"))
fig <- fig + theme(legend.margin=margin(rep(2.5,4))) + theme(plot.margin = margin(t = 1.5,r = 1.6,b = 0,l = 1.5))

pdf(file='Simu_MG1.pdf',width=4,height=3.8)
fig
dev.off()

##Figure S8
#Bernoulli
BA_n400 <- cbind(B402C[,7:8],B412C[,7:8],B422C[,7:8],B442C[,7:8],B482C[,7:8])
BA_n400 <- data.frame(median=apply(BA_n400,2,median),Low=apply(BA_n400,2,function(x)quantile(x,probs=0.05)),Up=apply(BA_n400,2,function(x)quantile(x,probs=0.95)),method=rep(c('MultiNeSS','MultiNeSS+'),5))
BA_n400$T <- as.factor(sort(rep(c(5,10,20,40,80),2)))
BA_n400$method <- factor(BA_n400$method,levels=c('MultiNeSS','MultiNeSS+'))
fig <- ggplot(data=BA_n400,aes(x=T,y=median,group=method,color=method,linetype=method)) + geom_line(size=0.85) + geom_point(size=1)
fig <- fig + geom_errorbar(aes(ymin = Low,ymax = Up),linetype = 1,show.legend = FALSE,width = 0.1) + scale_linetype_manual(values=1:4)
fig <- fig + scale_y_continuous(trans = 'log',breaks = c(1,10,100),limits = c(0.45,200)) + theme(text = element_text(size = 17))
fig <- fig + theme(legend.position = c(0.75,0.94)) + theme(legend.title = element_blank()) 
fig <- fig + labs(y=expression(paste('||',widehat(ZZ^T) - Z^'*',Z^'*T' - W[5]^'*',W[5]^'*T','||'[F]^2,'/',n))) 
fig <- fig + theme(panel.grid.minor = element_blank()) + theme(axis.title = element_text(size=16)) + theme(legend.key.width = unit(2.3,"line"))
fig <- fig + theme(legend.margin=margin(rep(2.5,4))) + theme(plot.margin = margin(t = 1.5,r = 1.6,b = 0,l = 1.5))

pdf(file='Simu_MB2.pdf',width=4,height=3.8)
fig
dev.off()

#Gaussian
BA_n400 <- cbind(G402C[,7:8],G412C[,7:8],G422C[,7:8],G442C[,7:8],G482C[,7:8])
BA_n400 <- data.frame(median=apply(BA_n400,2,median),Low=apply(BA_n400,2,function(x)quantile(x,probs=0.05)),Up=apply(BA_n400,2,function(x)quantile(x,probs=0.95)),method=rep(c('MultiNeSS','MultiNeSS+'),5))
BA_n400$T <- as.factor(sort(rep(c(5,10,20,40,80),2)))
BA_n400$method <- factor(BA_n400$method,levels=c('MultiNeSS','MultiNeSS+'))
fig <- ggplot(data=BA_n400,aes(x=T,y=median,group=method,color=method,linetype=method)) + geom_line(size=0.85) + geom_point(size=1)
fig <- fig + geom_errorbar(aes(ymin = Low,ymax = Up),linetype = 1,show.legend = FALSE,width = 0.1) + scale_linetype_manual(values=1:4)
fig <- fig + scale_y_continuous(trans = 'log',breaks = c(1,10,100),limits = c(0.09,200)) + theme(text = element_text(size = 17))
fig <- fig + theme(legend.position = c(0.75,0.94)) + theme(legend.title = element_blank()) 
fig <- fig + labs(y=expression(paste('||',widehat(ZZ^T) - Z^'*',Z^'*T' - W[5]^'*',W[5]^'*T','||'[F]^2,'/',n))) 
fig <- fig + theme(panel.grid.minor = element_blank()) + theme(axis.title = element_text(size=16)) + theme(legend.key.width = unit(2.3,"line"))
fig <- fig + theme(legend.margin=margin(rep(2.5,4))) + theme(plot.margin = margin(t = 1.5,r = 1.6,b = 0,l = 1.5))

pdf(file='Simu_MG2.pdf',width=4,height=3.8)
fig
dev.off()
                                                                                                                  
##Table S2
#Bernoulli
apply(B402C[,c(3,13:17)],2,mean)
apply(B412C[,c(3,13:22)],2,mean)
apply(B422C[,c(3,13:32)],2,mean)
apply(B442C[,c(3,13:52)],2,mean)
apply(B482C[,c(3,13:92)],2,mean)
 
#Gaussian
apply(G402C[,c(3,13:17)],2,mean)
apply(G412C[,c(3,13:22)],2,mean)
apply(G422C[,c(3,13:32)],2,mean)
apply(G442C[,c(3,13:52)],2,mean)
apply(G482C[,c(3,13:92)],2,mean)                                                                                                                     
