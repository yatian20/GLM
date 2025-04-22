#Data is available at https://www.stats.ox.ac.uk/~snijders/siena/Lazega_lawyers_data.htm
Lwork <- read.table('ELwork.dat')
Lfriend <- read.table('ELfriend.dat')
Ladv <- read.table('ELadv.dat')
Lattr <- read.table('ELattr.dat')

#You need to download the R file functions.R in your working path.
source("functions.R")
library(ggplot2)

#obtain networks with n = 71 and T = 3
Law <- array(0,dim = c(71,71,3))
Law[,,1] <- as.matrix(Lwork)
Law[,,2] <- (as.matrix(Lfriend) + t(as.matrix(Lfriend)))/2
Law[,,3] <- (as.matrix(Ladv) + t(as.matrix(Ladv)))/2

#covariates: office, practice, and status
office <- rep(0,71)
office[Lattr$V4 == 1] <- 'Boston'
office[Lattr$V4 == 2] <- 'Hartford'
office[Lattr$V4 == 3] <- 'Providence'

practice <- rep(0,71)
practice[Lattr$V7 == 1] <- 'litigation'
practice[Lattr$V7 == 2] <- 'corporate'

status <- rep(0,71)
status[Lattr$V2 == 1] <- 'partner'
status[Lattr$V2 == 2] <- 'associate'

#Fit Bernoulli model using proposed method
Law_k <- EST_k2(Law,'bernoulli')
fit_Law <- SILR(Law,Law_k$est_k,Law_k$est_kw,'bernoulli')

#compare with MultiNeSS
Law_P <- Law != 0
library(multiness)
fit_Lawp <- multiness_fit(Law_P,model="logistic",self_loops=FALSE,tuning="adaptive",refit=TRUE,optim_opts=list(return_posns=TRUE))

#compare with COSIE
fit_Lawc <- MASE(Law,2)

save(fit_Law,fit_Lawp,fit_Lawc,file = "lawyers.rda")

#plot the estimates of Z colored by office
Z_pca <- fit_Law$Z_hat %*% eigen(t(fit_Law$Z_hat) %*% fit_Law$Z_hat)$vectors
Z_Law <- as.data.frame(Z_pca)
names(Z_Law) <- c("Z1","Z2")
Z_Law$office <- factor(office)

fig <- ggplot() + geom_point(aes(x=Z1,y=-Z2,col=office,shape=office),data = Z_Law,size = 1)
fig <- fig + labs(x="1st component",y="2nd component")
fig <- fig + theme(legend.position = c(0.82,0.89)) + theme(legend.title=element_blank()) + theme(text = element_text(size = 15))
