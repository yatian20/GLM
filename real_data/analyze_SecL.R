library(multiness)
library(doParallel)
library(foreach)

set.seed(1)

# Load Lazega data
Lwork   <- read.table("real_data/raw_data/ELwork.dat")
Lfriend <- read.table("real_data/raw_data/ELfriend.dat")
Ladv    <- read.table("real_data/raw_data/ELadv.dat")
Lattr   <- read.table("real_data/raw_data/ELattr.dat")

n <- 71
m <- 3       # number of layers
T <- m       # number of layers

A_work   <- as.matrix(Lwork)
A_friend <- as.matrix(Lfriend)
A_adv    <- as.matrix(Ladv)

#### multiness only works for symmetric network. 
#### to make the network symmetric, let the edge between i & j = 1 when at least one of A_{ij} and A_{ji} = 1

A_work   <- (pmax(A_work,   t(A_work))   > 0) * 1
A_friend <- (pmax(A_friend, t(A_friend)) > 0) * 1
A_adv    <- (pmax(A_adv,    t(A_adv))    > 0) * 1

Law <- array(0, dim = c(n, n, m))
Law[,,1] <- A_work
Law[,,2] <- A_friend
Law[,,3] <- A_adv

for (k in 1:m) diag(Law[,,k]) <- 0


compute_auc <- function(fpr, tpr) {
  ord <- order(fpr)
  fpr <- fpr[ord]
  tpr <- tpr[ord]
  ok  <- is.finite(fpr) & is.finite(tpr)
  fpr <- fpr[ok]
  tpr <- tpr[ok]
  if (length(fpr) < 2) return(NA_real_)
  df  <- diff(fpr)
  mid <- (tpr[-1] + tpr[-length(tpr)]) / 2
  sum(df * mid)
}

psd_part <- function(M, tol = 0) {
  ee  <- eigen(M, symmetric = TRUE)
  val <- pmax(ee$values, tol)      # truncate below tol
  if (!any(val > 0)) return(matrix(0, nrow(M), ncol(M)))
  Vp  <- ee$vectors[, val > 0, drop = FALSE]
  Dp  <- val[val > 0]
  # Reconstruct PSD part: V diag(D) V'
  Vp %*% (Dp * t(Vp))              # efficient: diag(Dp) %*% t(Vp) = Dp * t(Vp)
}


replications <- 100
pi_obs       <- 0.9               # proportion observed
cutoffs      <- seq(0, 1, by = 0.005)
clen         <- length(cutoffs)

# Parallel backend
cores <- 13
cl <- makeCluster(cores)
registerDoParallel(cl, cores = cores)

# Make sure workers see multiness, SILR_pred4, and other global variables
clusterEvalQ(cl, {
  library(multiness)
  source("algorithms/latent_vectors_est.R")
})
clusterExport(cl, c("compute_auc","psd_part","Law","n","m","T","pi_obs","cutoffs","clen"))

# Layout of each row of Err_all:
#  1:clen              = SS-refinement TPR
#  clen+1:2clen        = SS-refinement FPR
#  2clen+1:3clen       = MultiNeSS+ TPR
#  3clen+1:4clen       = MultiNeSS+ FPR
#  6clen+7             = d1
#  6clen+8:6clen+10   = d2[1:3]
# Total length per row: 4clen+4

Err_all <- foreach(realization = 1:replications,
                   .combine = 'rbind',
                   .inorder = FALSE,
                   .packages = c("multiness")) %dopar% {
                     set.seed(realization*25)
                     
                     # Build mask Omega and partially observed networks Law_ob
                     Omega  <- array(0, dim = c(n, n, m))
                     Law_ob <- array(0, dim = c(n, n, m))
                     for (k in 1:m) {
                       M <- matrix(rbinom(n * n, size = 1, prob = pi_obs), n, n)
                       M[upper.tri(M)] <- t(M)[upper.tri(M)]
                       diag(M) <- 0
                       Omega[,,k]  <- M
                       Law_ob[,,k] <- Law[,,k] * M
                     }
                     
                     missing_pattern <- (Omega == 1)   # TRUE = observed, FALSE = held-out
                     
                     # MultiNeSS+ fitting
                     fit_refit <- multiness_fit(
                       A          = Law_ob,
                       model      = "logistic",
                       self_loops = FALSE,
                       refit      = TRUE, 
                       tuning     = "fixed", 
                       optim_opts = list(
                         missing_pattern = missing_pattern,
                         max_rank        = 20,
                         check_obj       = FALSE,
                         verbose         = FALSE
                       )
                     )
                     
                     test_pos <- (Omega == 0) & (Law == 1)
                     test_neg <- (Omega == 0) & (Law == 0)
                     n_pos    <- sum(test_pos)
                     n_neg    <- sum(test_neg)
                     
                     ## MultiNeSS+ predictions
                     ### only keep the positive kernel part for a fair comparison
                     pred_mn_refit <- array(0, dim = c(n, n, m))
                     F_pos <- psd_part(fit_refit$F_hat, tol = 0)
                     d1_silr <- sum(svd(F_pos)$d>1e-8) #Rank(F_pos)
                     d2_silr <- NULL
                     for (k in 1:m) {
                       G_pos   <- psd_part(fit_refit$G_hat[[k]], tol = 0)
                       d2_silr <- c(d2_silr, sum(svd(G_pos)$d>1e-8)) #Rank(G_pos)
                       Mk_psd  <- F_pos + G_pos
                       Pk_psd  <- plogis(Mk_psd)
                       Pk_psd  <- (Pk_psd + t(Pk_psd))/2
                       diag(Pk_psd) <- 0
                       pred_mn_refit[,,k] <- Pk_psd
                     }
                     
                     mn_tp_refit <- numeric(clen)
                     mn_fp_refit <- numeric(clen)
                     
                     for (i in 1:clen) {
                       thr      <- cutoffs[i]
                       pred_bin <- (pred_mn_refit > thr)
                       
                       mn_tp_refit[i] <- if (n_pos > 0)
                         sum(test_pos & pred_bin) / n_pos else NA_real_
                       mn_fp_refit[i] <- if (n_neg > 0)
                         sum(test_neg & pred_bin) / n_neg else NA_real_
                     }
                     
                     ## --- SS-refinement predictions --- 
                     est <- SILR_pred4(
                       A      = Law_ob,
                       k      = d1_silr,
                       kw     = d2_silr,
                       family = "bernoulli",
                       Omega  = Omega,
                       pi     = pi_obs
                     )

                     pred_1 <- array(0, dim = c(n, n, T)) 
                     for (t in 1:T) {
                       pred_1[,,t] <- 1 / (1 + exp(-est$Z_hat  %*% t(est$Z_hat)  -
                                                     est$W_hat[[t]]  %*% t(est$W_hat[[t]])))
                       diag(pred_1[,,t]) <- 0
                     }

                     silr_tp_ref  <- numeric(clen)
                     silr_fp_ref  <- numeric(clen)
                     
                     for (i in 1:clen) {
                       thr <- cutoffs[i]
                       
                       silr_tp_ref[i]  <- if (n_pos > 0)
                         sum(test_pos & (pred_1 > thr)) / n_pos else NA_real_
                       silr_fp_ref[i]  <- if (n_neg > 0)
                         sum(test_neg & (pred_1 > thr)) / n_neg else NA_real_
                     }
                     
                     ## --- Collect everything into one long vector ---
                     c(
                       silr_tp_ref,
                       silr_fp_ref,
                       mn_tp_refit,
                       mn_fp_refit
                     )
                   }

stopCluster(cl)

LkPdct_SILR_MNplus_Lazega <- Err_all

save(
  LkPdct_SILR_MNplus_Lazega,
  file = "LkPdct_SILR_MNplus_Lazega.rda"
)


