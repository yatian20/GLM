##Simulation for the estimation of Z/Wt
#fix n = 400, k = kt = 2; vary T = 5/10/20/40/80, Case A/B/C
#Case A: Z and Wt are orthogonal, i.e. phi = rho = 0
#Case B: Z and Wt are non-orthogonal with phi = 0.1, rho = 0.3
#Case C: Z and 5 Wt are are orthogonal, others are same

source('GLM/algorithms/latent_vectors_est.R')

##### Case A #####
#Poisson model
library(foreach)
library(doParallel)
rep <- 100
cores <- 50
cl <- makeCluster(cores)
registerDoParallel(cl, cores=cores)
P402A <- foreach(i = 1:rep,.combine='rbind') %dopar% {
  n <- 400
  T <- 5
  k <- 2
  
  #generate Z
  Z <- matrix(rnorm(10*n*k),10*n,k)
  Z <- Z[apply(Z^2,1,sum) < k,][1:n,]
  
  #generate Wt's
  Y <- Z
  for(t in 1:T){
    Wt <- matrix(rnorm(10*n*k),10*n,k)
    Wt <- Wt[apply(Wt^2,1,sum) < k,][1:n,]
    Y <- cbind(Y,Wt)
  }
  Gram0 <- t(Y) %*% Y
  
  #add dependency
  phi <- 0
  rho <- 0
  
  Gram <- matrix(rho,1+T,1+T)
  Gram[,1] <- phi
  Gram[1,] <- phi
  diag(Gram) <- 1
  Gram <- kronecker(Gram,diag(rep(1,k)))
  Gram <- (n/(2*sqrt(k)))*Gram
  
  EigG <- eigen(Gram)
  EigG0 <- eigen(Gram0)
  Y <- Y %*% (EigG0$vectors %*% diag(EigG0$values^{-1/2}) %*% t(EigG0$vectors)) %*% (EigG$vectors %*% diag(EigG$values^{1/2}) %*% t(EigG$vectors))
  Z <- Y[,1:k]
  W <- list()
  for(t in 1:T)
    W[[t]] <- Y[,(cumsum(rep(k,T+1))[t]+1):cumsum(rep(k,T+1))[t+1]]
  
  #generate A
  P <- array(0,dim = c(n,n,T))
  A <- array(0,dim = c(n,n,T))
  for(t in 1:T){
    P[,,t] <- exp(Z %*% t(Z) + W[[t]] %*% t(W[[t]]))
    A[,,t] <- matrix(rpois(n*n,as.vector(P[,,t])),n,n)
    A[,,t][upper.tri(A[,,t])] <- t(A[,,t])[upper.tri(A[,,t])]
  }
  
  est <- SILR(A,k,rep(k,T),'poisson')
  O <- (svd(t(est$Z_hat) %*% Z)$v) %*% t(svd(t(est$Z_hat) %*% Z)$u)
  Or <- (svd(t(est$Z_init) %*% Z)$v) %*% t(svd(t(est$Z_init) %*% Z)$u)
  
  err <- NULL
  err <- c(err,sum((est$Z_init - Z %*% Or)^2))
  err <- c(err,sum((est$Z_hat - Z %*% O)^2))
  err <- c(err,sum((est$Z_init %*% t(est$Z_init) - Z %*% t(Z))^2)/n)
  err <- c(err,sum((est$Z_hat %*% t(est$Z_hat) - Z %*% t(Z))^2)/n)
  err <- c(err,sum((est$W_init[[1]] %*% t(est$W_init[[1]]) - W[[1]] %*% t(W[[1]]))^2)/n)
  err <- c(err,sum((est$W_hat[[1]] %*% t(est$W_hat[[1]]) - W[[1]] %*% t(W[[1]]))^2)/n)
  return(err)
}
stopCluster(cl)
colMeans(P402A)
save(P402A,file="P402A.rda")

#Bernoulli model
library(foreach)
library(doParallel)
rep <- 100
cores <- 50
cl <- makeCluster(cores)
registerDoParallel(cl, cores=cores)
B402A <- foreach(i = 1:rep,.combine='rbind') %dopar% {
  n <- 400
  T <- 5
  k <- 2
  
  #generate Z
  Z <- matrix(rnorm(10*n*k),10*n,k)
  Z <- Z[apply(Z^2,1,sum) < k,][1:n,]
  
  #generate Wt's
  Y <- Z
  for(t in 1:T){
    Wt <- matrix(rnorm(10*n*k),10*n,k)
    Wt <- Wt[apply(Wt^2,1,sum) < k,][1:n,]
    Y <- cbind(Y,Wt)
  }
  Gram0 <- t(Y) %*% Y
  
  #add dependency
  phi <- 0
  rho <- 0
  
  Gram <- matrix(rho,1+T,1+T)
  Gram[,1] <- phi
  Gram[1,] <- phi
  diag(Gram) <- 1
  Gram <- kronecker(Gram,diag(rep(1,k)))
  Gram <- (n/(2*sqrt(k)))*Gram
  
  EigG <- eigen(Gram)
  EigG0 <- eigen(Gram0)
  Y <- Y %*% (EigG0$vectors %*% diag(EigG0$values^{-1/2}) %*% t(EigG0$vectors)) %*% (EigG$vectors %*% diag(EigG$values^{1/2}) %*% t(EigG$vectors))
  Z <- Y[,1:k]
  W <- list()
  for(t in 1:T)
    W[[t]] <- Y[,(cumsum(rep(k,T+1))[t]+1):cumsum(rep(k,T+1))[t+1]]
  
  #generate A
  P <- array(0,dim = c(n,n,T))
  A <- array(0,dim = c(n,n,T))
  for(t in 1:T){
    P[,,t] <- 1/(1 + exp(-Z %*% t(Z) - W[[t]] %*% t(W[[t]])))
    A[,,t] <- matrix(rbinom(n*n,1,as.vector(P[,,t])),n,n)
    A[,,t][upper.tri(A[,,t])] <- t(A[,,t])[upper.tri(A[,,t])]
  }
  
  est <- SILR(A,k,rep(k,T),'bernoulli')
  O <- (svd(t(est$Z_hat) %*% Z)$v) %*% t(svd(t(est$Z_hat) %*% Z)$u)
  Or <- (svd(t(est$Z_init) %*% Z)$v) %*% t(svd(t(est$Z_init) %*% Z)$u)
  err <- NULL
  err <- c(err,sum((est$Z_init - Z %*% Or)^2))
  err <- c(err,sum((est$Z_hat - Z %*% O)^2))
  
  #compare the errors of G
  library(multiness)
  est_mac <- multiness_fit(A,model="logistic",tuning="adaptive",refit=FALSE)
  est_macp <- multiness_fit(A,model="logistic",tuning="adaptive",refit=TRUE)
  err <- c(err,sum((est$Z_init %*% t(est$Z_init) - Z %*% t(Z))^2)/n)
  err <- c(err,sum((est$Z_hat %*% t(est$Z_hat) - Z %*% t(Z))^2)/n)
  err <- c(err,sum((est_mac$F_hat - Z %*% t(Z))^2)/n)
  err <- c(err,sum((est_macp$F_hat - Z %*% t(Z))^2)/n)
  
  #compare the errors of F1
  err <- c(err,sum((est$W_init[[1]] %*% t(est$W_init[[1]]) - W[[1]] %*% t(W[[1]]))^2)/n)
  err <- c(err,sum((est$W_hat[[1]] %*% t(est$W_hat[[1]]) - W[[1]] %*% t(W[[1]]))^2)/n)
  err <- c(err,sum((est_mac$G_hat[[1]] - W[[1]] %*% t(W[[1]]))^2)/n)
  err <- c(err,sum((est_macp$G_hat[[1]] - W[[1]] %*% t(W[[1]]))^2)/n)
  return(err)
}
stopCluster(cl)
colMeans(B402A)
save(B402A,file="B402A.rda")

#Gaussian model
library(foreach)
library(doParallel)
rep <- 100
cores <- 50
cl <- makeCluster(cores)
registerDoParallel(cl, cores=cores)
G402A <- foreach(i = 1:rep,.combine='rbind') %dopar% {
  n <- 400
  T <- 5
  k <- 2
  
  #generate Z
  Z <- matrix(rnorm(10*n*k),10*n,k)
  Z <- Z[apply(Z^2,1,sum) < k,][1:n,]
  
  #generate Wt's
  Y <- Z
  for(t in 1:T){
    Wt <- matrix(rnorm(10*n*k),10*n,k)
    Wt <- Wt[apply(Wt^2,1,sum) < k,][1:n,]
    Y <- cbind(Y,Wt)
  }
  Gram0 <- t(Y) %*% Y
  
  #add dependency
  phi <- 0
  rho <- 0
  
  Gram <- matrix(rho,1+T,1+T)
  Gram[,1] <- phi
  Gram[1,] <- phi
  diag(Gram) <- 1
  Gram <- kronecker(Gram,diag(rep(1,k)))
  Gram <- (n/(2*sqrt(k)))*Gram
  
  EigG <- eigen(Gram)
  EigG0 <- eigen(Gram0)
  Y <- Y %*% (EigG0$vectors %*% diag(EigG0$values^{-1/2}) %*% t(EigG0$vectors)) %*% (EigG$vectors %*% diag(EigG$values^{1/2}) %*% t(EigG$vectors))
  Z <- Y[,1:k]
  W <- list()
  for(t in 1:T)
    W[[t]] <- Y[,(cumsum(rep(k,T+1))[t]+1):cumsum(rep(k,T+1))[t+1]]
  
  #generate A
  P <- array(0,dim = c(n,n,T))
  A <- array(0,dim = c(n,n,T))
  for(t in 1:T){
    P[,,t] <- Z %*% t(Z) + W[[t]] %*% t(W[[t]])
    A[,,t] <- matrix(rnorm(n*n,as.vector(P[,,t])),n,n)
    A[,,t][upper.tri(A[,,t])] <- t(A[,,t])[upper.tri(A[,,t])]
  }
  
  est <- SILR(A,k,rep(k,T),'gaussian')
  O <- (svd(t(est$Z_hat) %*% Z)$v) %*% t(svd(t(est$Z_hat) %*% Z)$u)
  Or <- (svd(t(est$Z_init) %*% Z)$v) %*% t(svd(t(est$Z_init) %*% Z)$u)
  err <- NULL
  err <- c(err,sum((est$Z_init - Z %*% Or)^2))
  err <- c(err,sum((est$Z_hat - Z %*% O)^2))
  
  #compare the errors of G
  library(multiness)
  est_mac <- multiness_fit(A,model="gaussian",tuning="adaptive",refit=FALSE)
  est_macp <- multiness_fit(A,model="gaussian",tuning="adaptive",refit=TRUE)
  err <- c(err,sum((est$Z_init %*% t(est$Z_init) - Z %*% t(Z))^2)/n)
  err <- c(err,sum((est$Z_hat %*% t(est$Z_hat) - Z %*% t(Z))^2)/n)
  err <- c(err,sum((est_mac$F_hat - Z %*% t(Z))^2)/n)
  err <- c(err,sum((est_macp$F_hat - Z %*% t(Z))^2)/n)
  
  #compare the errors of F1
  err <- c(err,sum((est$W_init[[1]] %*% t(est$W_init[[1]]) - W[[1]] %*% t(W[[1]]))^2)/n)
  err <- c(err,sum((est$W_hat[[1]] %*% t(est$W_hat[[1]]) - W[[1]] %*% t(W[[1]]))^2)/n)
  err <- c(err,sum((est_mac$G_hat[[1]] - W[[1]] %*% t(W[[1]]))^2)/n)
  err <- c(err,sum((est_macp$G_hat[[1]] - W[[1]] %*% t(W[[1]]))^2)/n)
  return(err)
}
stopCluster(cl)
colMeans(G402A)
save(G402A,file="G402A.rda")

##### Case B #####
#Poisson model
library(foreach)
library(doParallel)
rep <- 100
cores <- 50
cl <- makeCluster(cores)
registerDoParallel(cl, cores=cores)
P402B <- foreach(i = 1:rep,.combine='rbind') %dopar% {
  n <- 400
  T <- 5
  k <- 2
  
  #generate Z
  Z <- matrix(rnorm(10*n*k),10*n,k)
  Z <- Z[apply(Z^2,1,sum) < k,][1:n,]
  
  #generate Wt's
  Y <- Z
  for(t in 1:T){
    Wt <- matrix(rnorm(10*n*k),10*n,k)
    Wt <- Wt[apply(Wt^2,1,sum) < k,][1:n,]
    Y <- cbind(Y,Wt)
  }
  Gram0 <- t(Y) %*% Y
  
  #add dependency
  phi <- 0.1
  rho <- 0.3
  
  Gram <- matrix(rho,1+T,1+T)
  Gram[,1] <- phi
  Gram[1,] <- phi
  diag(Gram) <- 1
  Gram <- kronecker(Gram,diag(rep(1,k)))
  Gram <- (n/(2*sqrt(k)))*Gram
  
  EigG <- eigen(Gram)
  EigG0 <- eigen(Gram0)
  Y <- Y %*% (EigG0$vectors %*% diag(EigG0$values^{-1/2}) %*% t(EigG0$vectors)) %*% (EigG$vectors %*% diag(EigG$values^{1/2}) %*% t(EigG$vectors))
  Z <- Y[,1:k]
  W <- list()
  for(t in 1:T)
    W[[t]] <- Y[,(cumsum(rep(k,T+1))[t]+1):cumsum(rep(k,T+1))[t+1]]
  
  #generate A
  P <- array(0,dim = c(n,n,T))
  A <- array(0,dim = c(n,n,T))
  for(t in 1:T){
    P[,,t] <- exp(Z %*% t(Z) + W[[t]] %*% t(W[[t]]))
    A[,,t] <- matrix(rpois(n*n,as.vector(P[,,t])),n,n)
    A[,,t][upper.tri(A[,,t])] <- t(A[,,t])[upper.tri(A[,,t])]
  }
  
  est <- SILR(A,k,rep(k,T),'poisson')
  O <- (svd(t(est$Z_hat) %*% Z)$v) %*% t(svd(t(est$Z_hat) %*% Z)$u)
  Or <- (svd(t(est$Z_init) %*% Z)$v) %*% t(svd(t(est$Z_init) %*% Z)$u)
  
  err <- NULL
  err <- c(err,sum((est$Z_init - Z %*% Or)^2))
  err <- c(err,sum((est$Z_hat - Z %*% O)^2))
  err <- c(err,sum((est$Z_init %*% t(est$Z_init) - Z %*% t(Z))^2)/n)
  err <- c(err,sum((est$Z_hat %*% t(est$Z_hat) - Z %*% t(Z))^2)/n)
  err <- c(err,sum((est$W_init[[1]] %*% t(est$W_init[[1]]) - W[[1]] %*% t(W[[1]]))^2)/n)
  err <- c(err,sum((est$W_hat[[1]] %*% t(est$W_hat[[1]]) - W[[1]] %*% t(W[[1]]))^2)/n)
  return(err)
}
stopCluster(cl)
colMeans(P402B)
save(P402B,file="P402B.rda")

#Bernoulli model
library(foreach)
library(doParallel)
rep <- 100
cores <- 50
cl <- makeCluster(cores)
registerDoParallel(cl, cores=cores)
B402B <- foreach(i = 1:rep,.combine='rbind') %dopar% {
  n <- 400
  T <- 5
  k <- 2
  
  #generate Z
  Z <- matrix(rnorm(10*n*k),10*n,k)
  Z <- Z[apply(Z^2,1,sum) < k,][1:n,]
  
  #generate Wt's
  Y <- Z
  for(t in 1:T){
    Wt <- matrix(rnorm(10*n*k),10*n,k)
    Wt <- Wt[apply(Wt^2,1,sum) < k,][1:n,]
    Y <- cbind(Y,Wt)
  }
  Gram0 <- t(Y) %*% Y
  
  #add dependency
  phi <- 0.1
  rho <- 0.3
  
  Gram <- matrix(rho,1+T,1+T)
  Gram[,1] <- phi
  Gram[1,] <- phi
  diag(Gram) <- 1
  Gram <- kronecker(Gram,diag(rep(1,k)))
  Gram <- (n/(2*sqrt(k)))*Gram
  
  EigG <- eigen(Gram)
  EigG0 <- eigen(Gram0)
  Y <- Y %*% (EigG0$vectors %*% diag(EigG0$values^{-1/2}) %*% t(EigG0$vectors)) %*% (EigG$vectors %*% diag(EigG$values^{1/2}) %*% t(EigG$vectors))
  Z <- Y[,1:k]
  W <- list()
  for(t in 1:T)
    W[[t]] <- Y[,(cumsum(rep(k,T+1))[t]+1):cumsum(rep(k,T+1))[t+1]]
  
  #generate A
  P <- array(0,dim = c(n,n,T))
  A <- array(0,dim = c(n,n,T))
  for(t in 1:T){
    P[,,t] <- 1/(1 + exp(-Z %*% t(Z) - W[[t]] %*% t(W[[t]])))
    A[,,t] <- matrix(rbinom(n*n,1,as.vector(P[,,t])),n,n)
    A[,,t][upper.tri(A[,,t])] <- t(A[,,t])[upper.tri(A[,,t])]
  }
  
  est <- SILR(A,k,rep(k,T),'bernoulli')
  O <- (svd(t(est$Z_hat) %*% Z)$v) %*% t(svd(t(est$Z_hat) %*% Z)$u)
  Or <- (svd(t(est$Z_init) %*% Z)$v) %*% t(svd(t(est$Z_init) %*% Z)$u)
  err <- NULL
  err <- c(err,sum((est$Z_init - Z %*% Or)^2))
  err <- c(err,sum((est$Z_hat - Z %*% O)^2))
  
  #compare the errors of G
  library(multiness)
  est_mac <- multiness_fit(A,model="logistic",tuning="adaptive",refit=FALSE)
  est_macp <- multiness_fit(A,model="logistic",tuning="adaptive",refit=TRUE)
  err <- c(err,sum((est$Z_init %*% t(est$Z_init) - Z %*% t(Z))^2)/n)
  err <- c(err,sum((est$Z_hat %*% t(est$Z_hat) - Z %*% t(Z))^2)/n)
  err <- c(err,sum((est_mac$F_hat - Z %*% t(Z))^2)/n)
  err <- c(err,sum((est_macp$F_hat - Z %*% t(Z))^2)/n)
  
  #compare the errors of F1
  err <- c(err,sum((est$W_init[[1]] %*% t(est$W_init[[1]]) - W[[1]] %*% t(W[[1]]))^2)/n)
  err <- c(err,sum((est$W_hat[[1]] %*% t(est$W_hat[[1]]) - W[[1]] %*% t(W[[1]]))^2)/n)
  err <- c(err,sum((est_mac$G_hat[[1]] - W[[1]] %*% t(W[[1]]))^2)/n)
  err <- c(err,sum((est_macp$G_hat[[1]] - W[[1]] %*% t(W[[1]]))^2)/n)
  return(err)
}
stopCluster(cl)
colMeans(B402B)
save(B402B,file="B402B.rda")

#Gaussian model
library(foreach)
library(doParallel)
rep <- 100
cores <- 50
cl <- makeCluster(cores)
registerDoParallel(cl, cores=cores)
G402B <- foreach(i = 1:rep,.combine='rbind') %dopar% {
  n <- 400
  T <- 5
  k <- 2
  
  #generate Z
  Z <- matrix(rnorm(10*n*k),10*n,k)
  Z <- Z[apply(Z^2,1,sum) < k,][1:n,]
  
  #generate Wt's
  Y <- Z
  for(t in 1:T){
    Wt <- matrix(rnorm(10*n*k),10*n,k)
    Wt <- Wt[apply(Wt^2,1,sum) < k,][1:n,]
    Y <- cbind(Y,Wt)
  }
  Gram0 <- t(Y) %*% Y
  
  #add dependency
  phi <- 0.1
  rho <- 0.3
  
  Gram <- matrix(rho,1+T,1+T)
  Gram[,1] <- phi
  Gram[1,] <- phi
  diag(Gram) <- 1
  Gram <- kronecker(Gram,diag(rep(1,k)))
  Gram <- (n/(2*sqrt(k)))*Gram
  
  EigG <- eigen(Gram)
  EigG0 <- eigen(Gram0)
  Y <- Y %*% (EigG0$vectors %*% diag(EigG0$values^{-1/2}) %*% t(EigG0$vectors)) %*% (EigG$vectors %*% diag(EigG$values^{1/2}) %*% t(EigG$vectors))
  Z <- Y[,1:k]
  W <- list()
  for(t in 1:T)
    W[[t]] <- Y[,(cumsum(rep(k,T+1))[t]+1):cumsum(rep(k,T+1))[t+1]]
  
  #generate A
  P <- array(0,dim = c(n,n,T))
  A <- array(0,dim = c(n,n,T))
  for(t in 1:T){
    P[,,t] <- Z %*% t(Z) + W[[t]] %*% t(W[[t]])
    A[,,t] <- matrix(rnorm(n*n,as.vector(P[,,t])),n,n)
    A[,,t][upper.tri(A[,,t])] <- t(A[,,t])[upper.tri(A[,,t])]
  }
  
  est <- SILR(A,k,rep(k,T),'gaussian')
  O <- (svd(t(est$Z_hat) %*% Z)$v) %*% t(svd(t(est$Z_hat) %*% Z)$u)
  Or <- (svd(t(est$Z_init) %*% Z)$v) %*% t(svd(t(est$Z_init) %*% Z)$u)
  err <- NULL
  err <- c(err,sum((est$Z_init - Z %*% Or)^2))
  err <- c(err,sum((est$Z_hat - Z %*% O)^2))
  
  #compare the errors of G
  library(multiness)
  est_mac <- multiness_fit(A,model="gaussian",tuning="adaptive",refit=FALSE)
  est_macp <- multiness_fit(A,model="gaussian",tuning="adaptive",refit=TRUE)
  err <- c(err,sum((est$Z_init %*% t(est$Z_init) - Z %*% t(Z))^2)/n)
  err <- c(err,sum((est$Z_hat %*% t(est$Z_hat) - Z %*% t(Z))^2)/n)
  err <- c(err,sum((est_mac$F_hat - Z %*% t(Z))^2)/n)
  err <- c(err,sum((est_macp$F_hat - Z %*% t(Z))^2)/n)
  
  #compare the errors of F1
  err <- c(err,sum((est$W_init[[1]] %*% t(est$W_init[[1]]) - W[[1]] %*% t(W[[1]]))^2)/n)
  err <- c(err,sum((est$W_hat[[1]] %*% t(est$W_hat[[1]]) - W[[1]] %*% t(W[[1]]))^2)/n)
  err <- c(err,sum((est_mac$G_hat[[1]] - W[[1]] %*% t(W[[1]]))^2)/n)
  err <- c(err,sum((est_macp$G_hat[[1]] - W[[1]] %*% t(W[[1]]))^2)/n)
  return(err)
}
stopCluster(cl)
colMeans(G402B)
save(G402B,file="G402B.rda")

##### Case C #####
#Poisson model
library(foreach)
library(doParallel)
rep <- 100
cores <- 50
cl <- makeCluster(cores)
registerDoParallel(cl, cores=cores)
P402C <- foreach(i = 1:rep,.combine='rbind') %dopar% {
  n <- 400
  T <- 5
  k <- 2
  
  #generate Z
  Z <- matrix(rnorm(10*n*k),10*n,k)
  Z <- Z[apply(Z^2,1,sum) < k,][1:n,]
  
  #generate Wt's
  Y <- Z
  for(t in 1:5){
    Wt <- matrix(rnorm(10*n*k),10*n,k)
    Wt <- Wt[apply(Wt^2,1,sum) < k,][1:n,]
    Y <- cbind(Y,Wt)
  }
  Gram0 <- t(Y) %*% Y
  
  #add dependency
  phi <- 0
  rho <- 0
  
  Gram <- matrix(rho,6,6)
  Gram[,1] <- phi
  Gram[1,] <- phi
  diag(Gram) <- 1
  Gram <- kronecker(Gram,diag(rep(1,k)))
  Gram <- (n/(2*sqrt(k)))*Gram
  
  EigG <- eigen(Gram)
  EigG0 <- eigen(Gram0)
  Y <- Y %*% (EigG0$vectors %*% diag(EigG0$values^{-1/2}) %*% t(EigG0$vectors)) %*% (EigG$vectors %*% diag(EigG$values^{1/2}) %*% t(EigG$vectors))
  Z <- Y[,1:k]
  W <- list()
  for(t in 1:T){
    if(t <= 5)
      W[[t]] <- Y[,(cumsum(rep(k,T+1))[t]+1):cumsum(rep(k,T+1))[t+1]]
    else
      W[[t]] <- W[[5]]
  }
    
  #generate A
  P <- array(0,dim = c(n,n,T))
  A <- array(0,dim = c(n,n,T))
  for(t in 1:T){
    P[,,t] <- exp(Z %*% t(Z) + W[[t]] %*% t(W[[t]]))
    A[,,t] <- matrix(rpois(n*n,as.vector(P[,,t])),n,n)
    A[,,t][upper.tri(A[,,t])] <- t(A[,,t])[upper.tri(A[,,t])]
  }
  
  est <- SILR(A,k,rep(k,T),'poisson')
  O <- (svd(t(est$Z_hat) %*% Z)$v) %*% t(svd(t(est$Z_hat) %*% Z)$u)
  Or <- (svd(t(est$Z_init) %*% Z)$v) %*% t(svd(t(est$Z_init) %*% Z)$u)
  
  err <- NULL
  err <- c(err,sum((est$Z_init - Z %*% Or)^2))
  err <- c(err,sum((est$Z_hat - Z %*% O)^2))
  err <- c(err,sum((est$Z_init %*% t(est$Z_init) - Z %*% t(Z))^2)/n)
  err <- c(err,sum((est$Z_hat %*% t(est$Z_hat) - Z %*% t(Z))^2)/n)
  err <- c(err,sum((est$W_init[[1]] %*% t(est$W_init[[1]]) - W[[1]] %*% t(W[[1]]))^2)/n)
  err <- c(err,sum((est$W_hat[[1]] %*% t(est$W_hat[[1]]) - W[[1]] %*% t(W[[1]]))^2)/n)
  return(err)
}
stopCluster(cl)
colMeans(P402C)
save(P402C,file="P402C.rda")

#Bernoulli model
library(foreach)
library(doParallel)
rep <- 100
cores <- 50
cl <- makeCluster(cores)
registerDoParallel(cl, cores=cores)
B402C <- foreach(i = 1:rep,.combine='rbind') %dopar% {
  n <- 400
  T <- 5
  k <- 2
  
  #generate Z
  Z <- matrix(rnorm(10*n*k),10*n,k)
  Z <- Z[apply(Z^2,1,sum) < k,][1:n,]
  
  #generate Wt's
  Y <- Z
  for(t in 1:5){
    Wt <- matrix(rnorm(10*n*k),10*n,k)
    Wt <- Wt[apply(Wt^2,1,sum) < k,][1:n,]
    Y <- cbind(Y,Wt)
  }
  Gram0 <- t(Y) %*% Y
  
  #add dependency
  phi <- 0
  rho <- 0
  
  Gram <- matrix(rho,6,6)
  Gram[,1] <- phi
  Gram[1,] <- phi
  diag(Gram) <- 1
  Gram <- kronecker(Gram,diag(rep(1,k)))
  Gram <- (n/(2*sqrt(k)))*Gram
  
  EigG <- eigen(Gram)
  EigG0 <- eigen(Gram0)
  Y <- Y %*% (EigG0$vectors %*% diag(EigG0$values^{-1/2}) %*% t(EigG0$vectors)) %*% (EigG$vectors %*% diag(EigG$values^{1/2}) %*% t(EigG$vectors))
  Z <- Y[,1:k]
  W <- list()
  for(t in 1:T){
    if(t <= 5)
      W[[t]] <- Y[,(cumsum(rep(k,T+1))[t]+1):cumsum(rep(k,T+1))[t+1]]
    else
      W[[t]] <- W[[5]]
  }
  
  #generate A
  P <- array(0,dim = c(n,n,T))
  A <- array(0,dim = c(n,n,T))
  for(t in 1:T){
    P[,,t] <- 1/(1 + exp(-Z %*% t(Z) - W[[t]] %*% t(W[[t]])))
    A[,,t] <- matrix(rbinom(n*n,1,as.vector(P[,,t])),n,n)
    A[,,t][upper.tri(A[,,t])] <- t(A[,,t])[upper.tri(A[,,t])]
  }
  
  est <- SILR(A,k,rep(k,T),'bernoulli')
  O <- (svd(t(est$Z_hat) %*% Z)$v) %*% t(svd(t(est$Z_hat) %*% Z)$u)
  Or <- (svd(t(est$Z_init) %*% Z)$v) %*% t(svd(t(est$Z_init) %*% Z)$u)
  err <- NULL
  err <- c(err,sum((est$Z_init - Z %*% Or)^2))
  err <- c(err,sum((est$Z_hat - Z %*% O)^2))
  
  #compare the errors of G
  library(multiness)
  est_mac <- multiness_fit(A,model="logistic",tuning="adaptive",refit=FALSE)
  est_macp <- multiness_fit(A,model="logistic",tuning="adaptive",refit=TRUE)
  err <- c(err,sum((est$Z_init %*% t(est$Z_init) - Z %*% t(Z))^2)/n)
  err <- c(err,sum((est$Z_hat %*% t(est$Z_hat) - Z %*% t(Z))^2)/n)
  err <- c(err,sum((est_mac$F_hat - Z %*% t(Z))^2)/n)
  err <- c(err,sum((est_macp$F_hat - Z %*% t(Z))^2)/n)
  
  #compare the errors of F1
  err <- c(err,sum((est$W_init[[1]] %*% t(est$W_init[[1]]) - W[[1]] %*% t(W[[1]]))^2)/n)
  err <- c(err,sum((est$W_hat[[1]] %*% t(est$W_hat[[1]]) - W[[1]] %*% t(W[[1]]))^2)/n)
  err <- c(err,sum((est_mac$G_hat[[1]] - W[[1]] %*% t(W[[1]]))^2)/n)
  err <- c(err,sum((est_macp$G_hat[[1]] - W[[1]] %*% t(W[[1]]))^2)/n)
  return(err)
}
stopCluster(cl)
colMeans(B402C)
save(B402C,file="B402C.rda")

#Gaussian model
library(foreach)
library(doParallel)
rep <- 100
cores <- 50
cl <- makeCluster(cores)
registerDoParallel(cl, cores=cores)
G402C <- foreach(i = 1:rep,.combine='rbind') %dopar% {
  n <- 400
  T <- 5
  k <- 2
  
  #generate Z
  Z <- matrix(rnorm(10*n*k),10*n,k)
  Z <- Z[apply(Z^2,1,sum) < k,][1:n,]
  
  #generate Wt's
  Y <- Z
  for(t in 1:5){
    Wt <- matrix(rnorm(10*n*k),10*n,k)
    Wt <- Wt[apply(Wt^2,1,sum) < k,][1:n,]
    Y <- cbind(Y,Wt)
  }
  Gram0 <- t(Y) %*% Y
  
  #add dependency
  phi <- 0
  rho <- 0
  
  Gram <- matrix(rho,6,6)
  Gram[,1] <- phi
  Gram[1,] <- phi
  diag(Gram) <- 1
  Gram <- kronecker(Gram,diag(rep(1,k)))
  Gram <- (n/(2*sqrt(k)))*Gram
  
  EigG <- eigen(Gram)
  EigG0 <- eigen(Gram0)
  Y <- Y %*% (EigG0$vectors %*% diag(EigG0$values^{-1/2}) %*% t(EigG0$vectors)) %*% (EigG$vectors %*% diag(EigG$values^{1/2}) %*% t(EigG$vectors))
  Z <- Y[,1:k]
  W <- list()
  for(t in 1:T){
    if(t <= 5)
      W[[t]] <- Y[,(cumsum(rep(k,T+1))[t]+1):cumsum(rep(k,T+1))[t+1]]
    else
      W[[t]] <- W[[5]]
  }
  
  #generate A
  P <- array(0,dim = c(n,n,T))
  A <- array(0,dim = c(n,n,T))
  for(t in 1:T){
    P[,,t] <- Z %*% t(Z) + W[[t]] %*% t(W[[t]])
    A[,,t] <- matrix(rnorm(n*n,as.vector(P[,,t])),n,n)
    A[,,t][upper.tri(A[,,t])] <- t(A[,,t])[upper.tri(A[,,t])]
  }
  
  est <- SILR(A,k,rep(k,T),'gaussian')
  O <- (svd(t(est$Z_hat) %*% Z)$v) %*% t(svd(t(est$Z_hat) %*% Z)$u)
  Or <- (svd(t(est$Z_init) %*% Z)$v) %*% t(svd(t(est$Z_init) %*% Z)$u)
  err <- NULL
  err <- c(err,sum((est$Z_init - Z %*% Or)^2))
  err <- c(err,sum((est$Z_hat - Z %*% O)^2))
  
  #compare the errors of G
  library(multiness)
  est_mac <- multiness_fit(A,model="gaussian",tuning="adaptive",refit=FALSE)
  est_macp <- multiness_fit(A,model="gaussian",tuning="adaptive",refit=TRUE)
  err <- c(err,sum((est$Z_init %*% t(est$Z_init) - Z %*% t(Z))^2)/n)
  err <- c(err,sum((est$Z_hat %*% t(est$Z_hat) - Z %*% t(Z))^2)/n)
  err <- c(err,sum((est_mac$F_hat - Z %*% t(Z))^2)/n)
  err <- c(err,sum((est_macp$F_hat - Z %*% t(Z))^2)/n)
  
  #compare the errors of F1
  err <- c(err,sum((est$W_init[[1]] %*% t(est$W_init[[1]]) - W[[1]] %*% t(W[[1]]))^2)/n)
  err <- c(err,sum((est$W_hat[[1]] %*% t(est$W_hat[[1]]) - W[[1]] %*% t(W[[1]]))^2)/n)
  err <- c(err,sum((est_mac$G_hat[[1]] - W[[1]] %*% t(W[[1]]))^2)/n)
  err <- c(err,sum((est_macp$G_hat[[1]] - W[[1]] %*% t(W[[1]]))^2)/n)
  return(err)
}
stopCluster(cl)
colMeans(G402C)
save(G402C,file="G402C.rda")
