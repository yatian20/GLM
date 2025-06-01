#R function for estimating the latent dimensions
#Input: multiple networks At (n*n*T dimensional array) and data distribution ('poisson', 'bernoulli', or 'gaussian')
#Output: estimates for the shared latent dimension k and the heterogeneous latent dimensions kt
EST_k <- function(A,family){
  n <- dim(A)[1]
  T <- dim(A)[3]
  if(family == 'poisson'){
    g <- function(x) exp(x)
    g_inverse <- function(x) log(x)
  }
  
  if(family == 'bernoulli'){
    g <- function(x) 1/(1+exp(-x))
    g_inverse <- function(x) log(x/(1-x))
  }
  
  if(family == 'gaussian'){
    g <- function(x) x
    g_inverse <- function(x) x
  }
  
  #Stage1: estimate k+kt
  est_ky <- NULL
  Y <- list()
  for(t in 1:T){
    if(family != 'gaussian'){
      svd_At <- svd(A[,,t])
      tau_t <- sqrt(sum(A[,,t])/n)
      P_hat <- svd_At$u %*% diag(svd_At$d * I(svd_At$d > tau_t)) %*% t(svd_At$v)
      P_hat <- g(-5) * I(P_hat < g(-5)) + P_hat * I(g(-5) <= P_hat & P_hat <= g(5)) + g(5) * I(P_hat > g(5))
      Theta_hat <- g_inverse(P_hat)
    }
    else
      Theta_hat <- A[,,t]
    Theta_hat <- eigen(Theta_hat)$vectors %*% diag(eigen(Theta_hat)$values * I(eigen(Theta_hat)$values > 0)) %*% t(eigen(Theta_hat)$vectors)
    EV <- eigen(Theta_hat)$values[1:(n/10 + 1)]
    est_ky <- c(est_ky,which(EV[1:(n/10) + 1]/EV[1:(n/10)] < n^(-1/(4*1:(n/10)+8)))[1])
    Y0 <- eigen(Theta_hat)$vectors %*% rbind(diag(sqrt(eigen(Theta_hat)$values[1:est_ky[t]])),matrix(0,n-est_ky[t],est_ky[t]))
    
    eta <- 0.01 / (svd(Y0)$d[1]^2)
    grad0_Theta <- A[,,t] - g(Y0 %*% t(Y0))
    grad0 <- grad0_Theta %*% Y0
    Y1 <- Y0 + eta * grad0
    for(i in 1:999){
      grad1_Theta <- A[,,t] - g(Y1 %*% t(Y1))
      grad1 <- grad1_Theta %*% Y1
      eta <- -sum(diag(t(Y1-Y0) %*% (grad1-grad0)))/sum(diag(t(grad1-grad0) %*% (grad1-grad0)))
      if(is.na(eta))
        break
      Y0 <- Y1
      grad0 <- grad1
      Y1 <- Y0 + eta * grad0
    }
    Y[[t]] <- Y1
  }
  
  #Stage2: estimate k and kt
  k_est <- NULL
  for(t1 in 1:(T-1))
    for(t2 in (t1+1):T){
      Yc <- cbind(Y[[t1]],Y[[t2]])
      EV <- c(svd(Yc)$d,0)^2
      k_est <- c(k_est,ncol(Yc) - which(EV[2:(ncol(Yc)+1)]/EV[1] < n^(-1/4))[1])
    }
  est_k <- min(k_est)
  est_kw <- est_ky - est_k
  return(list(est_k = est_k,est_kw = est_kw))
}

#R function for estimating the latent dimensions (slightly modified version for data analysis)
#Input: multiple networks At (n*n*T dimensional array) and data distribution ('poisson', 'bernoulli', or 'gaussian')
#Output: estimates for the shared latent dimension k and the heterogeneous latent dimensions kt
EST_k2 <- function(A,family){
  n <- dim(A)[1]
  T <- dim(A)[3]
  if(family == 'poisson'){
    g <- function(x) exp(x)
    g_inverse <- function(x) log(x)
  }
  
  if(family == 'bernoulli'){
    g <- function(x) 1/(1+exp(-x))
    g_inverse <- function(x) log(x/(1-x))
  }
  
  if(family == 'gaussian'){
    g <- function(x) x
    g_inverse <- function(x) x
  }
  
  #Stage1: estimate k+kt
  est_ky <- NULL
  Y <- list()
  for(t in 1:T){
    if(family != 'gaussian'){
      svd_At <- svd(A[,,t])
      tau_t <- sqrt(sum(A[,,t])/n)
      P_hat <- svd_At$u %*% diag(svd_At$d * I(svd_At$d > tau_t)) %*% t(svd_At$v)
      P_hat <- g(-5) * I(P_hat < g(-5)) + P_hat * I(g(-5) <= P_hat & P_hat <= g(5)) + g(5) * I(P_hat > g(5))
      Theta_hat <- g_inverse(P_hat)
    }
    else
      Theta_hat <- A[,,t]
    Theta_hat <- eigen(Theta_hat)$vectors %*% diag(eigen(Theta_hat)$values * I(eigen(Theta_hat)$values > 0)) %*% t(eigen(Theta_hat)$vectors)
    EV <- eigen(Theta_hat)$values[1:(n/5 + 1)]
    est_ky <- c(est_ky,which(EV[1:(n/5) + 1]/EV[1:(n/5)] < n^(-1/(4*1:(n/5)+8)) & cumsum(EV[1:(n/5)])/sum(eigen(Theta_hat)$values) > 0.5)[1])
    Y0 <- eigen(Theta_hat)$vectors %*% rbind(diag(sqrt(eigen(Theta_hat)$values[1:est_ky[t]])),matrix(0,n-est_ky[t],est_ky[t]))
    
    eta <- 0.01 / (svd(Y0)$d[1]^2)
    grad0_Theta <- A[,,t] - g(Y0 %*% t(Y0))
    grad0 <- grad0_Theta %*% Y0
    Y1 <- Y0 + eta * grad0
    for(i in 1:999){
      grad1_Theta <- A[,,t] - g(Y1 %*% t(Y1))
      grad1 <- grad1_Theta %*% Y1
      eta <- -sum(diag(t(Y1-Y0) %*% (grad1-grad0)))/sum(diag(t(grad1-grad0) %*% (grad1-grad0)))
      if(is.na(eta))
        break
      Y0 <- Y1
      grad0 <- grad1
      Y1 <- Y0 + eta * grad0
    }
    Y[[t]] <- Y1
  }
  
  #Stage2: estimate k and kt
  k_est <- NULL
  for(t1 in 1:(T-1))
    for(t2 in (t1+1):T){
      Yc <- cbind(Y[[t1]],Y[[t2]])
      EV <- c(svd(Yc)$d,0)^2
      k_est <- c(k_est,ncol(Yc) - which(EV[2:(ncol(Yc)+1)]/EV[est_ky[t2]] < n^(-1/4))[1])
    }
  est_k <- min(k_est)
  est_kw <- est_ky - est_k
  return(list(est_k = est_k,est_kw = est_kw))
}
