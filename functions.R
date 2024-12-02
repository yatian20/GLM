#R function for implementing SS-Hunting and SS-Refinement (initialized by USVT and PGD)
#Input: multiple networks At (n*n*T dimensional array), shared latent dimension k (numeric), heterogeneous latent dimensions kt (T dimensional vector), and data distribution ('poisson', 'bernoulli', or 'gaussian')
#Output: estimates for Z and Wt by SS-Hunting, estimates for Z and Wt by SS-Refinement
SILR <- function(A,k,kw,family){
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
  
  #Stage1: USVT and PGD for each At
  Y <- list()
  for(t in 1:T){
    if(family != 'gaussian'){
      svd_At <- svd(A[,,t])
      tau_t <- sqrt(sum(A[,,t])/n)
      P_hat <- svd_At$u %*% diag(svd_At$d * (svd_At$d > tau_t)) %*% t(svd_At$v)
      P_hat <- g(-5) * (P_hat < g(-5)) + P_hat * (g(-5) <= P_hat & P_hat <= g(5)) + g(5) * (P_hat > g(5))
      Theta_hat <- g_inverse(P_hat)
    }
    else
      Theta_hat <- A[,,t]
    Theta_hat <- eigen(Theta_hat)$vectors %*% diag(eigen(Theta_hat)$values * (eigen(Theta_hat)$values > 0)) %*% t(eigen(Theta_hat)$vectors)
    Y0 <- eigen(Theta_hat)$vectors %*% rbind(diag(sqrt(eigen(Theta_hat)$values[1:(k+kw[t])])),matrix(0,n-(k+kw[t]),(k+kw[t])))
    
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
  
  #Stage2: distinguish Z and Wt
  num_pair <- 0
  G0 <- matrix(0,n,n)
  for(t1 in 1:(T-1))
    for(t2 in (t1+1):T){
      Yc <- cbind(Y[[t1]],Y[[t2]])
      EV <- c(svd(Yc)$d,0)^2
      if(EV[1]/EV[k+kw[t1]+kw[t2]] < 2*log(n)){
        num_pair <- num_pair + 1
        Vc <- svd(Yc)$v[,(k+kw[t1]+kw[t2]+1):(2*k+kw[t1]+kw[t2])]
        G0 <- G0 + cbind(Y[[t1]],-Y[[t2]]) %*% Vc %*% t(Vc) %*% t(cbind(Y[[t1]],-Y[[t2]]))/2
      }
    }
  G0 <- G0/num_pair
  
  if(k == 1)
    Zr <- eigen(G0)$vectors %*% rbind(diag(as.matrix(sqrt(eigen(G0)$values[1:k]))),matrix(0,n-k,k))
  else
    Zr <- eigen(G0)$vectors %*% rbind(diag(sqrt(eigen(G0)$values[1:k])),matrix(0,n-k,k))
  
  Wr <- list()
  for(t in 1:T){
    Ft <- Y[[t]] %*% t(Y[[t]]) - G0
    if(kw[t] == 1)
      Wr[[t]] <- eigen(Ft)$vectors %*% rbind(diag(as.matrix(sqrt(eigen(Ft)$values[1:kw[t]]))),matrix(0,n-kw[t],kw[t]))
    else 
      Wr[[t]] <- eigen(Ft)$vectors %*% rbind(diag(sqrt(eigen(Ft)$values[1:kw[t]])),matrix(0,n-kw[t],kw[t]))
  }
  
  #Stage3: PGD for Z and Wt
  Z0 <- Zr
  W0 <- Wr
  W1 <- list()
  eta_z <- 0.01 / (T*svd(Z0)$d[1]^2)
  grad0_z <- matrix(0,n,k)
  grad0_w <- list()
  for(t in 1:T){
    eta_wt <- 0.01 / (svd(W0[[t]])$d[1]^2)
    grad0_t <- A[,,t] - g(Z0 %*% t(Z0) + W0[[t]] %*% t(W0[[t]]))
    grad0_w[[t]] <- grad0_t %*% W0[[t]]
    W1[[t]] <- W0[[t]] + eta_wt * grad0_w[[t]]
    grad0_z <- grad0_z + grad0_t %*% Z0
  }
  Z1 <- Z0 + eta_z * grad0_z
  
  for(i in 1:999){
    grad1_z <- matrix(0,n,k)
    grad1_w <- list()
    for(t in 1:T){
      grad1_t <- A[,,t] - g(Z1 %*% t(Z1) + W1[[t]] %*% t(W1[[t]]))
      grad1_w[[t]] <- grad1_t %*% W1[[t]]
      eta_wt <- -sum(diag(t(W1[[t]]-W0[[t]]) %*% (grad1_w[[t]]-grad0_w[[t]])))/sum(diag(t(grad1_w[[t]]-grad0_w[[t]]) %*% (grad1_w[[t]]-grad0_w[[t]])))
      if(is.na(eta_wt))
        break
      W0[[t]] <- W1[[t]]
      grad0_w[[t]] <- grad1_w[[t]]
      W1[[t]] <- W0[[t]] + eta_wt * grad0_w[[t]]
      grad1_z <- grad1_z + grad1_t %*% Z1
    }
    eta_z <- -sum(diag(t(Z1-Z0) %*% (grad1_z-grad0_z)))/sum(diag(t(grad1_z-grad0_z) %*% (grad1_z-grad0_z)))
    if(is.na(eta_wt) | is.na(eta_z))
      break
    Z0 <- Z1
    grad0_z <- grad1_z
    Z1 <- Z0 + eta_z * grad0_z
  }
  return(list(Z_init = Zr,W_init = Wr,Z_hat = Z1,W_hat = W1))
}

#R function for implementing SS-Hunting and SS-Refinement (initialized by USVT and PGD, incorporate the second-order update)
#Input: multiple networks At (n*n*T dimensional array), shared latent dimension k (numeric), heterogeneous latent dimensions kt (T dimensional vector), and data distribution ('poisson', 'bernoulli', or 'gaussian')
#Output: estimates for Z and Wt by SS-Hunting, estimates for Z and Wt by SS-Refinement
SILR2 <- function(A,k,kw,family){
  n <- dim(A)[1]
  T <- dim(A)[3]
  if(family == 'poisson'){
    g <- function(x) exp(x)
    g_inverse <- function(x) log(x)
    nu_pp <- function(x) exp(x)
  }
  
  if(family == 'bernoulli'){
    g <- function(x) 1/(1+exp(-x))
    g_inverse <- function(x) log(x/(1-x))
    nu_pp <- function(x) exp(x)/(1+exp(x))^2
  }
  
  if(family == 'gaussian'){
    g <- function(x) x
    g_inverse <- function(x) x
    nu_pp <- function(x) !is.na(x)
  }
  
  #Stage1: USVT and PGD for each At
  Y <- list()
  for(t in 1:T){
    if(family != 'gaussian'){
      svd_At <- svd(A[,,t])
      tau_t <- sqrt(sum(A[,,t])/n)
      P_hat <- svd_At$u %*% diag(svd_At$d * (svd_At$d > tau_t)) %*% t(svd_At$v)
      P_hat <- g(-5) * (P_hat < g(-5)) + P_hat * (g(-5) <= P_hat & P_hat <= g(5)) + g(5) * (P_hat > g(5))
      Theta_hat <- g_inverse(P_hat)
    }
    else
      Theta_hat <- A[,,t]
    Theta_hat <- eigen(Theta_hat)$vectors %*% diag(eigen(Theta_hat)$values * (eigen(Theta_hat)$values > 0)) %*% t(eigen(Theta_hat)$vectors)
    Y0 <- eigen(Theta_hat)$vectors %*% rbind(diag(sqrt(eigen(Theta_hat)$values[1:(k+kw[t])])),matrix(0,n-(k+kw[t]),(k+kw[t])))
    
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
  
  #Stage2: distinguish Z and Wt
  num_pair <- 0
  G0 <- matrix(0,n,n)
  for(t1 in 1:(T-1))
    for(t2 in (t1+1):T){
      Yc <- cbind(Y[[t1]],Y[[t2]])
      EV <- c(svd(Yc)$d,0)^2
      if(EV[1]/EV[k+kw[t1]+kw[t2]] < 2*log(n)){
        num_pair <- num_pair + 1
        Vc <- svd(Yc)$v[,(k+kw[t1]+kw[t2]+1):(2*k+kw[t1]+kw[t2])]
        G0 <- G0 + cbind(Y[[t1]],-Y[[t2]]) %*% Vc %*% t(Vc) %*% t(cbind(Y[[t1]],-Y[[t2]]))/2
      }
    }
  G0 <- G0/num_pair
  
  if(k == 1)
    Zr <- eigen(G0)$vectors %*% rbind(diag(as.matrix(sqrt(eigen(G0)$values[1:k]))),matrix(0,n-k,k))
  else
    Zr <- eigen(G0)$vectors %*% rbind(diag(sqrt(eigen(G0)$values[1:k])),matrix(0,n-k,k))
  
  Wr <- list()
  for(t in 1:T){
    Ft <- Y[[t]] %*% t(Y[[t]]) - G0
    if(kw[t] == 1)
      Wr[[t]] <- eigen(Ft)$vectors %*% rbind(diag(as.matrix(sqrt(eigen(Ft)$values[1:kw[t]]))),matrix(0,n-kw[t],kw[t]))
    else 
      Wr[[t]] <- eigen(Ft)$vectors %*% rbind(diag(sqrt(eigen(Ft)$values[1:kw[t]])),matrix(0,n-kw[t],kw[t]))
  }
  
  #Stage3: PGD for Z and Wt
  Z0 <- Zr
  W0 <- Wr
  W1 <- list()
  eta_z <- 0.01 / (T*svd(Z0)$d[1]^2)
  grad0_z <- matrix(0,n,k)
  grad0_w <- list()
  for(t in 1:T){
    eta_wt <- 0.01 / (svd(W0[[t]])$d[1]^2)
    grad0_t <- A[,,t] - g(Z0 %*% t(Z0) + W0[[t]] %*% t(W0[[t]]))
    grad0_w[[t]] <- grad0_t %*% W0[[t]]
    W1[[t]] <- W0[[t]] + eta_wt * grad0_w[[t]]
    grad0_z <- grad0_z + grad0_t %*% Z0
  }
  Z1 <- Z0 + eta_z * grad0_z
  
  for(i in 1:999){
    grad1_z <- matrix(0,n,k)
    grad1_w <- list()
    for(t in 1:T){
      grad1_t <- A[,,t] - g(Z1 %*% t(Z1) + W1[[t]] %*% t(W1[[t]]))
      grad1_w[[t]] <- grad1_t %*% W1[[t]]
      eta_wt <- -sum(diag(t(W1[[t]]-W0[[t]]) %*% (grad1_w[[t]]-grad0_w[[t]])))/sum(diag(t(grad1_w[[t]]-grad0_w[[t]]) %*% (grad1_w[[t]]-grad0_w[[t]])))
      if(is.na(eta_wt))
        break
      W0[[t]] <- W1[[t]]
      grad0_w[[t]] <- grad1_w[[t]]
      W1[[t]] <- W0[[t]] + eta_wt * grad0_w[[t]]
      grad1_z <- grad1_z + grad1_t %*% Z1
    }
    eta_z <- -sum(diag(t(Z1-Z0) %*% (grad1_z-grad0_z)))/sum(diag(t(grad1_z-grad0_z) %*% (grad1_z-grad0_z)))
    if(is.na(eta_wt) | is.na(eta_z))
      break
    Z0 <- Z1
    grad0_z <- grad1_z
    Z1 <- Z0 + eta_z * grad0_z
  }

  #Stage4: one-step for Z and Wt
  Z <- Z1
  W <- W1
  Ijoint <- matrix(0,n*(k+sum(kw)),n*(k+sum(kw)))
  for(t in 1:T){
    mu <- nu_pp(Z %*% t(Z) + W[[t]] %*% t(W[[t]]))
    
    Lzz <- matrix(0,n*k,n*k)
    for(i in 1:n)
      for(j in 1:n){
        if(i == j){
          Lzz[((i-1)*k+1):(i*k),((j-1)*k+1):(j*k)] <- 3*mu[i,i] * Z[i,] %*% t(Z[i,])
          for(jj in 1:n)
            Lzz[((i-1)*k+1):(i*k),((j-1)*k+1):(j*k)] <- Lzz[((i-1)*k+1):(i*k),((j-1)*k+1):(j*k)] + mu[i,jj]*Z[jj,] %*% t(Z[jj,])
        }
        else
          Lzz[((i-1)*k+1):(i*k),((j-1)*k+1):(j*k)] <- mu[i,j] * Z[j,] %*% t(Z[i,])
      }
    
    Lzw <- matrix(0,n*k,n*kw[t])
    for(i in 1:n)
      for(j in 1:n){
        if(i == j){
          Lzw[((i-1)*k+1):(i*k),((j-1)*kw[t]+1):(j*kw[t])] <- 3*mu[i,i] * Z[i,] %*% t(W[[t]][i,])
          for(jj in 1:n)
            Lzw[((i-1)*k+1):(i*k),((j-1)*kw[t]+1):(j*kw[t])] <- Lzw[((i-1)*k+1):(i*k),((j-1)*kw[t]+1):(j*kw[t])] + mu[i,jj]*Z[jj,] %*% t(W[[t]][jj,])
        }
        else
          Lzw[((i-1)*k+1):(i*k),((j-1)*kw[t]+1):(j*kw[t])] <- mu[i,j] * Z[j,] %*% t(W[[t]][i,])
      }
    
    Lww <- matrix(0,n*kw[t],n*kw[t])
    for(i in 1:n)
      for(j in 1:n){
        if(i == j){
          Lww[((i-1)*kw[t]+1):(i*kw[t]),((j-1)*kw[t]+1):(j*kw[t])] <- 3*mu[i,i] * W[[t]][i,] %*% t(W[[t]][i,])
          for(jj in 1:n)
            Lww[((i-1)*kw[t]+1):(i*kw[t]),((j-1)*kw[t]+1):(j*kw[t])] <- Lww[((i-1)*kw[t]+1):(i*kw[t]),((j-1)*kw[t]+1):(j*kw[t])] + mu[i,jj]*W[[t]][jj,] %*% t(W[[t]][jj,])
        }
        else
          Lww[((i-1)*kw[t]+1):(i*kw[t]),((j-1)*kw[t]+1):(j*kw[t])] <- mu[i,j] * W[[t]][j,] %*% t(W[[t]][i,])
      }
    
    ind1 <- 1:(n*k)
    ind2 <- (n*(k+cumsum(c(0,kw))[t])+1) : (n*(k+cumsum(c(0,kw))[t+1]))
    Ijoint[ind1,ind1] <- Ijoint[ind1,ind1] + Lzz
    Ijoint[ind2,ind2] <- Ijoint[ind2,ind2] + Lww
    Ijoint[ind1,ind2] <- Ijoint[ind1,ind2] + Lzw 
    Ijoint[ind2,ind1] <- Ijoint[ind2,ind1] + t(Lzw)
  }

  Sjoint <- rep(0,n*(k+sum(kw)))
  for(t in 1:T){
    M <- A[,,t] - g(Z %*% t(Z) + W[[t]] %*% t(W[[t]]))
    
    Lz <- rep(0,n*k)
    for(i in 1:n){
      Lz[((i-1)*k+1):(i*k)] <- Z[i,] * M[i,i]
      for(j in 1:n)
        Lz[((i-1)*k+1):(i*k)] <- Lz[((i-1)*k+1):(i*k)] + Z[j,] * M[i,j]
    }
    
    Lw <- rep(0,n*kw[t])
    for(i in 1:n){
      Lw[((i-1)*kw[t]+1):(i*kw[t])] <- W[[t]][i,] * M[i,i]
      for(j in 1:n)
        Lw[((i-1)*kw[t]+1):(i*kw[t])] <- Lw[((i-1)*kw[t]+1):(i*kw[t])] + W[[t]][j,] * M[i,j]
    }
    
    ind1 <- 1:(n*k)
    ind2 <- (n*(k+cumsum(c(0,kw))[t])+1) : (n*(k+cumsum(c(0,kw))[t+1]))
    Sjoint[ind1] <- Sjoint[ind1] + Lz
    Sjoint[ind2] <- Sjoint[ind2] + Lw
  }
  
  eigI <- eigen(Ijoint)
  rankI <- (n*k - k*(k-1)/2) + sum((n*kw - kw*(kw-1)/2))
  Ipi <- eigI$vectors %*% diag(c(eigI$values[1:rankI]^-1,rep(0,n*(k+sum(kw))-rankI))) %*% t(eigI$vectors)
  
  v <- as.vector(t(Z))
  for(t in 1:T)
    v <- c(v,as.vector(t(W[[t]])))
  v <- v + Ipi %*% Sjoint
  
  Z <- t(matrix(v[1:(n*k)],k,n))
  W <- list()
  for(t in 1:T){
    ind2 <- (n*(k+cumsum(c(0,kw))[t])+1) : (n*(k+cumsum(c(0,kw))[t+1]))
    W[[t]] <- t(matrix(v[ind2],kw[t],n))
  }
  return(list(Z_init = Zr,W_init = Wr,Z_hat = Z,W_hat = W))
}

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
