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

##Table S1 (provide only the last three lines as an example)                                                                                                       
load('B482Ak.rda')
load('B483Ak.rda')
load('B484Ak.rda')
load('G482Ak.rda')
load('G483Ak.rda')
load('G484Ak.rda')
load('P482Ak.rda')
load('P483Ak.rda')
load('P484Ak.rda')
load('B482Bk.rda')
load('B483Bk.rda')
load('B484Bk.rda')
load('G482Bk.rda')
load('G483Bk.rda')
load('G484Bk.rda')
load('P482Bk.rda')
load('P483Bk.rda')
load('P484Bk.rda')
load('B482Ck.rda')
load('B483Ck.rda')
load('B484Ck.rda')
load('G482Ck.rda')
load('G483Ck.rda')
load('G484Ck.rda')
load('P482Ck.rda')
load('P483Ck.rda')
load('P484Ck.rda')
                                                                                                                     
c(mean(apply(B482Ak,1,all)),mean(apply(B483Ak,1,all)),mean(apply(B484Ak,1,all)),mean(apply(G482Ak,1,all)),mean(apply(G483Ak,1,all)),mean(apply(G484Ak,1,all)),mean(apply(P482Ak,1,all)),mean(apply(P483Ak,1,all)),mean(apply(P484Ak,1,all)))
c(mean(apply(B482Bk,1,all)),mean(apply(B483Bk,1,all)),mean(apply(B484Bk,1,all)),mean(apply(G482Bk,1,all)),mean(apply(G483Bk,1,all)),mean(apply(G484Bk,1,all)),mean(apply(P482Bk,1,all)),mean(apply(P483Bk,1,all)),mean(apply(P484Bk,1,all)))
c(mean(apply(B482Ck,1,all)),mean(apply(B483Ck,1,all)),mean(apply(B484Ck,1,all)),mean(apply(G482Ck,1,all)),mean(apply(G483Ck,1,all)),mean(apply(G484Ck,1,all)),mean(apply(P482Ck,1,all)),mean(apply(P483Ck,1,all)),mean(apply(P484Ck,1,all))) 







## Figures S13 and S14
## ============================================================
## Six ROC plots (berA, berB, berC) × (T=5, T=80)
## Inputs (files under GLM/simulation2/results/):
##   LkPdct_MN+_<scenario>_T<T>.rda   -> object: MNplus_<scenario>_T<T>
##   LkPdct_SS_<scenario>_T<T>.rda    -> object: SS_<scenario>_T<T>
## Layouts:
##   MNplus_*   : [1:ncut]=MN+ TPR, [(ncut+1):(2*ncut)]=MN+ FPR
##   SS_*       : [1:ncut]=SS-refinement TPR, [ncut+1:2*ncut]=SS-refinement FPR,
##                 [2*ncut+1:3*ncut]=Oracle TPR, [3*ncut+1:4*ncut]=Oracle FPR
## Output:
##   plot_simupred4_<scenario>_mean_T<T>_pi09_MNplus_refit_only.pdf
## ============================================================

get_silr_curves <- function(M) {
  # M: reps x (4*ncut)
  stopifnot(is.matrix(M), ncol(M) %% 4L == 0L)
  ncut  <- ncol(M) %/% 4L
  mmean <- apply(M, 2, mean, na.rm = TRUE)
  list(
    prop = list( # Proposed (SS refinement)
      tpr = mmean[1:ncut],
      fpr = mmean[(ncut + 1):(2L * ncut)]
    ),
    orac = list( # Oracle
      tpr = mmean[(2L * ncut + 1L):(3L * ncut)],
      fpr = mmean[(3L * ncut + 1L):(4L * ncut)]
    )
  )
}

get_mnplus_curves <- function(M) {
  # M: reps x (2*ncut)
  stopifnot(is.matrix(M), ncol(M) %% 2L == 0L)
  ncut  <- ncol(M) %/% 2L
  mmean <- apply(M, 2, mean, na.rm = TRUE)
  list(
    tpr = mmean[1:ncut],
    fpr = mmean[(ncut + 1):(2L * ncut)]
  )
}

col_proposed <- "red";       lty_proposed <- 1; pch_proposed <- 1
col_mnplus   <- "blue";      lty_mnplus   <- 2; pch_mnplus   <- 2
col_oracle   <- "darkgreen"; lty_oracle   <- 3; pch_oracle   <- 4

scenarios <- c("berA", "berB", "berC")
Ts        <- c(5, 80)
base_dir  <- "./"

## ---- Loop: 3 scenarios × 2 T's = 6 plots ----
for (sc in scenarios) {
  for (TT in Ts) {
    
    ## Load SS (Proposed + Oracle)
    ss_file <- file.path(base_dir, sprintf("LkPdct_SS_%s_T%d.rda", sc, TT))
    if (!file.exists(ss_file)) {
      message("Missing SS file: ", ss_file, " (skipping)."); next
    }
    ss_env <- new.env(parent = emptyenv())
    load(ss_file, envir = ss_env)
    ss_obj_name <- sprintf("SS_%s_T%d", sc, TT)
    if (!exists(ss_obj_name, envir = ss_env)) {
      message("Object ", ss_obj_name, " not found in ", ss_file, " (skipping)."); next
    }
    Res_silr <- ss_env[[ss_obj_name]]
    
    ## Load MN+
    mn_file <- file.path(base_dir, sprintf("LkPdct_MN+_%s_T%d.rda", sc, TT))
    if (!file.exists(mn_file)) {
      message("Missing MN+ file: ", mn_file, " (skipping)."); next
    }
    mn_env <- new.env(parent = emptyenv())
    load(mn_file, envir = mn_env)
    mn_obj_name <- sprintf("MNplus_%s_T%d", sc, TT)
    if (!exists(mn_obj_name, envir = mn_env)) {
      message("Object ", mn_obj_name, " not found in ", mn_file, " (skipping)."); next
    }
    Err_mn_refit <- mn_env[[mn_obj_name]]
    
    silr   <- get_silr_curves(Res_silr)       # Proposed + Oracle
    mnplus <- get_mnplus_curves(Err_mn_refit) # MN+ (refit)
    
    out_pdf <- sprintf("plot_simupred4_%s_mean_T%d_pi09_MNplus_refit_only.pdf", sc, TT)
    pdf(out_pdf, width = 4, height = 4)
    par(mar = c(4, 4.5, 1, 1))
    par(pty = "s")
    
    plot(silr$orac$fpr, silr$orac$tpr,
         type = "o", col = col_oracle, pch = pch_oracle, lty = lty_oracle, lwd = 1,
         xlab = "False-positive rate", ylab = "True-positive rate",
         xlim = c(0, 1), ylim = c(0, 1),
         cex.axis = 1.2, cex.lab = 1.5)
    
    lines(mnplus$fpr, mnplus$tpr,
          type = "o", col = col_mnplus, pch = pch_mnplus, lty = lty_mnplus, lwd = 1)
    
    # Finally Proposed (SS refinement)
    lines(silr$prop$fpr, silr$prop$tpr,
          type = "o", col = col_proposed, pch = pch_proposed, lty = lty_proposed, lwd = 1)
    
    abline(0, 1, lty = 3, col = "gray70")
    
    legend("bottomright",
           legend = c("SS-refinement", "MultiNeSS+", "Oracle"),
           col    = c(col_proposed, col_mnplus, col_oracle),
           pch    = c(pch_proposed, pch_mnplus, pch_oracle),
           lty    = c(lty_proposed, lty_mnplus, lty_oracle),
           lwd    = 1,
           bty    = "n",
           cex    = 1.3)
    
    dev.off()
    message("Wrote ", out_pdf)
  }
}

