#This document shows the R code for our simulation studies (under Bernoulli distribution and Case (A))
#You need to download the R file functions.R in your working path.
source("functions.R")

#The following R code is for only one simulation. You can run the file repeatedly to estimate errors more accurately.
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
err <- NULL
err <- c(err,sum((est$Z_init %*% t(est$Z_init) - Z %*% t(Z))^2)/n)
err <- c(err,sum((est$Z_hat %*% t(est$Z_hat) - Z %*% t(Z))^2)/n)
err <- c(err,sum((est$W_init[[1]] %*% t(est$W_init[[1]]) - W[[1]] %*% t(W[[1]]))^2)/n)
err <- c(err,sum((est$W_hat[[1]] %*% t(est$W_hat[[1]]) - W[[1]] %*% t(W[[1]]))^2)/n)
