Lwork <- read.table('GLM/real_data/raw_data/ELwork.dat')
Lfriend <- read.table('GLM/real_data/raw_data/ELfriend.dat')
Ladv <- read.table('GLM/real_data/raw_data/ELadv.dat')
Lattr <- read.table('GLM/real_data/raw_data/ELattr.dat')

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

##Fit Bernoulli model using proposed method
source('GLM/algorithms/latent_vectors_est.R')
source('GLM/algorithms/latent_dimensions_est.R')

#slightly modify the algorithm for estimate k
Law_k <- EST_k2(Law,'bernoulli')
fit_Law <- SILR(Law,Law_k$est_k,Law_k$est_kw,'bernoulli')

#compare with MultiNeSS
Law_P <- Law != 0
library(multiness)
fit_Lawp <- multiness_fit(Law_P,model="logistic",self_loops=FALSE,tuning="adaptive",refit=TRUE,optim_opts=list(return_posns=TRUE))

#compare with COSIE
fit_Lawc <- MASE(Law,2)

save(fit_Law,fit_Lawp,fit_Lawc,office,practice,status,file = "lawyers.rda")
