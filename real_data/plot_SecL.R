## -------------------------------------------------------
## Plot SS-refinement and MultiNeSS+ ROC on Lazega
## -------------------------------------------------------

load("LkPdct_SILR_MNplus_Lazega.rda")
Res <- LkPdct_SILR_MNplus_Lazega

cutoffs <- seq(from = 0, to = 1, by = 0.005)
cl      <- length(cutoffs)

Res_mean <- apply(Res, 2, mean, na.rm = TRUE)

# SS-refine
tpr_silr_ref  <- Res_mean[(1):(cl)]
fpr_silr_ref  <- Res_mean[(cl + 1):(2 * cl)]

# MultiNeSS+
tpr_mn <- Res_mean[(2 * cl + 1):(3 * cl)]
fpr_mn <- Res_mean[(3 * cl + 1):(4 * cl)]

# ---------------- Plot ----------------

pdf("plot_Lazega_SILR_vs_MNplusPSD_pi09.pdf", width = 4, height = 4)

par(mar = c(4.5, 4.5, 1, 1))  
par(pty = "s")                # Force square plotting region

# SS-refine curve and MultiNeSS+ curve 
plot(fpr_mn, tpr_mn,
     type = "o", col = "blue", pch = 2, cex = 1, lty = 3, lwd = 1,
     xlab = "False-positive rate", ylab = "True-positive rate",
     xlim = c(0, 1), ylim = c(0, 1),
     cex.axis = 1.2, cex.lab = 1.5)

lines(fpr_silr_ref, tpr_silr_ref,
      type = "o", col = "red", pch = 1, cex = 1.5, lty = 1, lwd = 1)

# diagonal
abline(0, 1, lty = 3, col = "gray70")

legend("bottomright",
       legend = c("SS-refinement", "MultiNeSS+"),
       col    = c("red", "blue"),
       lty    = c(1, 3),
       pch    = c(1, 2),
       bty    = "n",
       cex    = 1.4)

dev.off()


