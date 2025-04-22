# =============================================================================
# SIMULATED DATA GENERATION: BERNOULLI DISTRIBUTION - CASE (A)
# =============================================================================
# Generates simulated binary data from a Bernoulli distribution for Case (A).
# Serves as a template that can be adapted for other simulation scenarios.
# The current implementation performs a single simulation. You can run the file repeatedly to estimate errors more accurately.

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

save(A,file = "simulated_data.rda")
  
