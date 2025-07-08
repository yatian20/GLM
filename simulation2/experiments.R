source('GLM/algorithms/latent_vectors_est.R')
source('GLM/algorithms/latent_dimensions_est.R')

##########################################################################
##Simulation for the estimation of k/kt
##########################################################################
#fix n = 400; vary T = 5/10/20/40/80, k = kt = 2/3/4, Case A/B/C
#only shows T = 80 and k = kt = 2 case as the primary example
#results under other T and k values can be easily repeated by changing T <- 80 and k <- 2 to other values 

##### Case A #####
#Poisson model
library(foreach)
library(doParallel)
rep <- 100
cores <- 50
cl <- makeCluster(cores)
registerDoParallel(cl, cores=cores)
EST <- foreach(i = 1:rep,.combine='rbind') %dopar% {
  n <- 400
  T <- 80
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
  
  est <- EST_k(A,'poisson')
  err <- c(est$est_k == k,all(est$est_kw + est$est_k == rep(2*k,T)))
  return(err)
}
stopCluster(cl)
colMeans(EST)

#Bernoulli model
library(foreach)
library(doParallel)
rep <- 100
cores <- 50
cl <- makeCluster(cores)
registerDoParallel(cl, cores=cores)
EST <- foreach(i = 1:rep,.combine='rbind') %dopar% {
  n <- 400
  T <- 80
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
  
  est <- EST_k(A,'bernoulli')
  err <- c(est$est_k == k,all(est$est_kw + est$est_k == rep(2*k,T)))
  return(err)
}
stopCluster(cl)
colMeans(EST)

#Gaussian model
library(foreach)
library(doParallel)
rep <- 100
cores <- 50
cl <- makeCluster(cores)
registerDoParallel(cl, cores=cores)
EST <- foreach(i = 1:rep,.combine='rbind') %dopar% {
  n <- 400
  T <- 80
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
  
  est <- EST_k(A,'gaussian')
  err <- c(est$est_k == k,all(est$est_kw + est$est_k == rep(2*k,T)))
  return(err)
}
stopCluster(cl)
colMeans(EST)

##### Case B #####
#Poisson model
library(foreach)
library(doParallel)
rep <- 100
cores <- 50
cl <- makeCluster(cores)
registerDoParallel(cl, cores=cores)
EST <- foreach(i = 1:rep,.combine='rbind') %dopar% {
  n <- 400
  T <- 80
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
  
  est <- EST_k(A,'poisson')
  err <- c(est$est_k == k,all(est$est_kw + est$est_k == rep(2*k,T)))
  return(err)
}
stopCluster(cl)
colMeans(EST)

#Bernoulli model
library(foreach)
library(doParallel)
rep <- 100
cores <- 50
cl <- makeCluster(cores)
registerDoParallel(cl, cores=cores)
EST <- foreach(i = 1:rep,.combine='rbind') %dopar% {
  n <- 400
  T <- 80
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
  
  est <- EST_k(A,'bernoulli')
  err <- c(est$est_k == k,all(est$est_kw + est$est_k == rep(2*k,T)))
  return(err)
}
stopCluster(cl)
colMeans(EST)

#Gaussian model
library(foreach)
library(doParallel)
rep <- 100
cores <- 50
cl <- makeCluster(cores)
registerDoParallel(cl, cores=cores)
EST <- foreach(i = 1:rep,.combine='rbind') %dopar% {
  n <- 400
  T <- 80
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
  
  est <- EST_k(A,'gaussian')
  err <- c(est$est_k == k,all(est$est_kw + est$est_k == rep(2*k,T)))
  return(err)
}
stopCluster(cl)
colMeans(EST)

##### Case C #####
#Poisson model
library(foreach)
library(doParallel)
rep <- 100
cores <- 50
cl <- makeCluster(cores)
registerDoParallel(cl, cores=cores)
EST <- foreach(i = 1:rep,.combine='rbind') %dopar% {
  n <- 400
  T <- 80
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
  
  est <- EST_k(A,'poisson')
  err <- c(est$est_k == k,all(est$est_kw + est$est_k == rep(2*k,T)))
  return(err)
}
stopCluster(cl)
colMeans(EST)

#Bernoulli model
library(foreach)
library(doParallel)
rep <- 100
cores <- 50
cl <- makeCluster(cores)
registerDoParallel(cl, cores=cores)
EST <- foreach(i = 1:rep,.combine='rbind') %dopar% {
  n <- 400
  T <- 80
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
  
  est <- EST_k(A,'bernoulli')
  err <- c(est$est_k == k,all(est$est_kw + est$est_k == rep(2*k,T)))
  return(err)
}
stopCluster(cl)
colMeans(EST)

#Gaussian model
library(foreach)
library(doParallel)
rep <- 100
cores <- 50
cl <- makeCluster(cores)
registerDoParallel(cl, cores=cores)
EST <- foreach(i = 1:rep,.combine='rbind') %dopar% {
  n <- 400
  T <- 80
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
  
  est <- EST_k(A,'gaussian')
  err <- c(est$est_k == k,all(est$est_kw + est$est_k == rep(2*k,T)))
  return(err)
}
stopCluster(cl)
colMeans(EST)

##########################################################################
##Comparison of Variants of Algorithm 2
##########################################################################
#fix n = 400, k = kt = 2, T=5; vary Case A/B/C, Bernoulli/Gaussian/Poisson distribution

##### Case A / C #####
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
  est2 <- SILR2(A,k,rep(k,T),'poisson')
  est3 <- SILR3(A,k,rep(k,T),'poisson')
  
  err <- NULL
  err <- c(err,sum((est$Z_init %*% t(est$Z_init) - Z %*% t(Z))^2)/n)
  err <- c(err,sum((est$Z_hat %*% t(est$Z_hat) - Z %*% t(Z))^2)/n)
  err <- c(err,sum((est2$Z_hat %*% t(est2$Z_hat) - Z %*% t(Z))^2)/n)
  err <- c(err,sum((est3$Z_hat %*% t(est3$Z_hat) - Z %*% t(Z))^2)/n)
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
  est2 <- SILR2(A,k,rep(k,T),'bernoulli')
  est3 <- SILR3(A,k,rep(k,T),'bernoulli')
  
  err <- NULL
  err <- c(err,sum((est$Z_init %*% t(est$Z_init) - Z %*% t(Z))^2)/n)
  err <- c(err,sum((est$Z_hat %*% t(est$Z_hat) - Z %*% t(Z))^2)/n)
  err <- c(err,sum((est2$Z_hat %*% t(est2$Z_hat) - Z %*% t(Z))^2)/n)
  err <- c(err,sum((est3$Z_hat %*% t(est3$Z_hat) - Z %*% t(Z))^2)/n)
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
  est2 <- SILR2(A,k,rep(k,T),'gaussian')
  est3 <- SILR3(A,k,rep(k,T),'gaussian')
  
  err <- NULL
  err <- c(err,sum((est$Z_init %*% t(est$Z_init) - Z %*% t(Z))^2)/n)
  err <- c(err,sum((est$Z_hat %*% t(est$Z_hat) - Z %*% t(Z))^2)/n)
  err <- c(err,sum((est2$Z_hat %*% t(est2$Z_hat) - Z %*% t(Z))^2)/n)
  err <- c(err,sum((est3$Z_hat %*% t(est3$Z_hat) - Z %*% t(Z))^2)/n)
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
  est2 <- SILR2(A,k,rep(k,T),'poisson')
  est3 <- SILR3(A,k,rep(k,T),'poisson')
  
  err <- NULL
  err <- c(err,sum((est$Z_init %*% t(est$Z_init) - Z %*% t(Z))^2)/n)
  err <- c(err,sum((est$Z_hat %*% t(est$Z_hat) - Z %*% t(Z))^2)/n)
  err <- c(err,sum((est2$Z_hat %*% t(est2$Z_hat) - Z %*% t(Z))^2)/n)
  err <- c(err,sum((est3$Z_hat %*% t(est3$Z_hat) - Z %*% t(Z))^2)/n)
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
  est2 <- SILR2(A,k,rep(k,T),'bernoulli')
  est3 <- SILR3(A,k,rep(k,T),'bernoulli')
  
  err <- NULL
  err <- c(err,sum((est$Z_init %*% t(est$Z_init) - Z %*% t(Z))^2)/n)
  err <- c(err,sum((est$Z_hat %*% t(est$Z_hat) - Z %*% t(Z))^2)/n)
  err <- c(err,sum((est2$Z_hat %*% t(est2$Z_hat) - Z %*% t(Z))^2)/n)
  err <- c(err,sum((est3$Z_hat %*% t(est3$Z_hat) - Z %*% t(Z))^2)/n)
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
  est2 <- SILR2(A,k,rep(k,T),'gaussian')
  est3 <- SILR3(A,k,rep(k,T),'gaussian')
  
  err <- NULL
  err <- c(err,sum((est$Z_init %*% t(est$Z_init) - Z %*% t(Z))^2)/n)
  err <- c(err,sum((est$Z_hat %*% t(est$Z_hat) - Z %*% t(Z))^2)/n)
  err <- c(err,sum((est2$Z_hat %*% t(est2$Z_hat) - Z %*% t(Z))^2)/n)
  err <- c(err,sum((est3$Z_hat %*% t(est3$Z_hat) - Z %*% t(Z))^2)/n)
  return(err)
}
stopCluster(cl)
colMeans(G402B)
save(G402B,file="G402B.rda")

##########################################################################
##Details for MultiNeSS and MultiNeSS+ under Case C
##########################################################################
#fix n = 400, k = kt = 2; vary T = 5/10/20/40/80, Bernoulli/Gaussian distribution
#only shows T = 80 case as the primary example
#results under other T values can be easily repeated by changing T <- 80 to other values 

#Bernoulli model
library(foreach)
library(doParallel)
rep <- 100
cores <- 50
cl <- makeCluster(cores)
registerDoParallel(cl, cores=cores)
B482C <- foreach(i = 1:rep,.combine='rbind') %dopar% {
  n <- 400
  T <- 80
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
  library(multiness)
  est_mac <- multiness_fit(A,model="logistic",tuning="adaptive",refit=FALSE)
  est_macp <- multiness_fit(A,model="logistic",tuning="adaptive",refit=TRUE)
  
  #check convergence
  err <- NULL
  err <- c(err,est_mac$convergence)
  err <- c(err,est_macp$convergence)
  
  #estimated shared latent dimension
  err <- c(err,est_mac$d1)
  err <- c(err,est_macp$d1)
   
  #errors of F_hat
  err <- c(err,sum((est_mac$F_hat - Z %*% t(Z))^2)/n)
  err <- c(err,sum((est_macp$F_hat - Z %*% t(Z))^2)/n)
  err <- c(err,sum((est_mac$F_hat - Z %*% t(Z) - W[[5]] %*% t(W[[5]]))^2)/n)
  err <- c(err,sum((est_macp$F_hat - Z %*% t(Z) - W[[5]] %*% t(W[[5]]))^2)/n)
  
  #errors of Theta_t 
  err_Theta <- NULL
  for(t in 1:T)
    err_Theta <- c(err_Theta,sum((est$Z_init %*% t(est$Z_init) + est$W_init[[t]] %*% t(est$W_init[[t]]) - Z %*% t(Z) - W[[t]] %*% t(W[[t]]))^2)/n)
  err <- c(err,max(err_Theta))
  
  err_Theta <- NULL
  for(t in 1:T)
    err_Theta <- c(err_Theta,sum((est$Z_hat %*% t(est$Z_hat) + est$W_hat[[t]] %*% t(est$W_hat[[t]]) - Z %*% t(Z) - W[[t]] %*% t(W[[t]]))^2)/n)
  err <- c(err,max(err_Theta))
  
  err_Theta <- NULL
  for(t in 1:T)
    err_Theta <- c(err_Theta,sum((est_mac$F_hat + est_mac$G_hat[[t]] - Z %*% t(Z) - W[[t]] %*% t(W[[t]]))^2)/n)
  err <- c(err,max(err_Theta))
  
  err_Theta <- NULL
  for(t in 1:T)
    err_Theta <- c(err_Theta,sum((est_macp$F_hat + est_macp$G_hat[[t]] - Z %*% t(Z) - W[[t]] %*% t(W[[t]]))^2)/n)
  err <- c(err,max(err_Theta))
  
  #estimated heterogeneous latent dimensions
  err <- c(err,est_mac$d2)
  err <- c(err,est_macp$d2)
  return(err)
}
stopCluster(cl)
colMeans(B482C)
save(B482C,file="B482C.rda")

#Gaussian model
library(foreach)
library(doParallel)
rep <- 100
cores <- 50
cl <- makeCluster(cores)
registerDoParallel(cl, cores=cores)
G482C <- foreach(i = 1:rep,.combine='rbind') %dopar% {
  n <- 400
  T <- 80
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
  library(multiness)
  est_mac <- multiness_fit(A,model="gaussian",tuning="adaptive",refit=FALSE)
  est_macp <- multiness_fit(A,model="gaussian",tuning="adaptive",refit=TRUE)
  
  #check convergence
  err <- NULL
  err <- c(err,est_mac$convergence)
  err <- c(err,est_macp$convergence)
  
  #estimated shared latent dimensions 
  err <- c(err,est_mac$d1)
  err <- c(err,est_macp$d1)
  
  #errors of F_hat
  err <- c(err,sum((est_mac$F_hat - Z %*% t(Z))^2)/n)
  err <- c(err,sum((est_macp$F_hat - Z %*% t(Z))^2)/n)
  err <- c(err,sum((est_mac$F_hat - Z %*% t(Z) - W[[5]] %*% t(W[[5]]))^2)/n)
  err <- c(err,sum((est_macp$F_hat - Z %*% t(Z) - W[[5]] %*% t(W[[5]]))^2)/n)
  
  #errors of Theta_t 
  err_Theta <- NULL
  for(t in 1:T)
    err_Theta <- c(err_Theta,sum((est$Z_init %*% t(est$Z_init) + est$W_init[[t]] %*% t(est$W_init[[t]]) - Z %*% t(Z) - W[[t]] %*% t(W[[t]]))^2)/n)
  err <- c(err,max(err_Theta))
  
  err_Theta <- NULL
  for(t in 1:T)
    err_Theta <- c(err_Theta,sum((est$Z_hat %*% t(est$Z_hat) + est$W_hat[[t]] %*% t(est$W_hat[[t]]) - Z %*% t(Z) - W[[t]] %*% t(W[[t]]))^2)/n)
  err <- c(err,max(err_Theta))
 
  err_Theta <- NULL
  for(t in 1:T)
    err_Theta <- c(err_Theta,sum((est_mac$F_hat + est_mac$G_hat[[t]] - Z %*% t(Z) - W[[t]] %*% t(W[[t]]))^2)/n)
  err <- c(err,max(err_Theta))

  err_Theta <- NULL
  for(t in 1:T)
    err_Theta <- c(err_Theta,sum((est_macp$F_hat + est_macp$G_hat[[t]] - Z %*% t(Z) - W[[t]] %*% t(W[[t]]))^2)/n)
  err <- c(err,max(err_Theta))
  
  #estimated heterogeneous latent dimensions
  err <- c(err,est_mac$d2)
  err <- c(err,est_macp$d2)
  return(err)
}
stopCluster(cl)
colMeans(G482C)
save(G482C,file="G482C.rda")
