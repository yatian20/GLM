library(ggplot2)
load('GLM/real_data/results/lawyers.rda')

##Figure 4
#SS-Refinement
Z_pca <- fit_Law$Z_hat %*% eigen(t(fit_Law$Z_hat) %*% fit_Law$Z_hat)$vectors
Z_Law <- as.data.frame(Z_pca)
names(Z_Law) <- c("Z1","Z2")
Z_Law$office <- factor(office)
Z_Law$practice <- factor(practice)

fig <- ggplot() + geom_point(aes(x=Z1,y=-Z2,col=office,shape=office),data = Z_Law,size = 2)
fig <- fig + labs(x="1st component",y="2nd component") + theme(axis.title = element_text(size=16))
fig <- fig + theme(legend.position = c(0.8,0.89)) + theme(legend.title=element_blank()) + theme(text = element_text(size = 17))
fig <- fig + theme(panel.grid.minor = element_blank())

pdf(file='Law_Z.pdf',width=4,height=3.9)
fig
dev.off()

#MultiNess
Z_Lawp <- as.data.frame(fit_Lawp$V_hat[,2:3])
names(Z_Lawp) <- c("Z1","Z2")
Z_Lawp$office <- factor(office)
Z_Lawp$practice <- factor(practice)

fig <- ggplot() + geom_point(aes(x=Z1,y=-Z2,col=office,shape=office),data = Z_Lawp,size = 2)
fig <- fig + labs(x="1st component",y="2nd component") + theme(axis.title = element_text(size=16))
fig <- fig + theme(legend.position = c(0.8,0.89)) + theme(legend.title=element_blank()) + theme(text = element_text(size = 17))
fig <- fig + theme(panel.grid.minor = element_blank())

pdf(file='Lawp_Z.pdf',width=4,height=3.9)
fig
dev.off()

#MASE
theta <- (1/3) * pi
U_hat <- fit_Lawc$U_hat %*% matrix(c(cos(theta),-sin(theta),sin(theta),cos(theta)),2,2)
U_Lawc <- as.data.frame(U_hat)
names(U_Lawc) <- c("U1","U2")
U_Lawc$office <- factor(office)
U_Lawc$practice <- factor(practice)

fig <- ggplot() + geom_point(aes(x=U1,y=-U2,col=office,shape=office),data = U_Lawc,size = 2)
fig <- fig + labs(x="1st component",y="2nd component") + scale_y_continuous(limits = c(0,0.225)) + theme(axis.title = element_text(size=16))
fig <- fig + theme(legend.position = c(0.8,0.89)) + theme(legend.title=element_blank()) + theme(text = element_text(size = 17))
fig <- fig + theme(panel.grid.minor = element_blank())

pdf(file='Lawc_Z.pdf',width=4,height=3.9)
fig
dev.off()

##Figure 5
#SS-Refinement
W1_pca <- fit_Law$W_hat[[1]] %*% eigen(t(fit_Law$W_hat[[1]]) %*% fit_Law$W_hat[[1]])$vectors
W1_Law <- as.data.frame(W1_pca)
names(W1_Law) <- c("W1","W2","W3","W4","W5","W6")
W1_Law$practice <- factor(practice)

fig <- ggplot() + geom_point(aes(x=W1,y=-W2,col=practice,shape=practice),data = W1_Law,size = 2)
fig <- fig + labs(x="1st component",y="2nd component") + theme(axis.title = element_text(size=16))
fig <- fig + theme(legend.position = c(0.8,0.92)) + theme(legend.title=element_blank()) + theme(text = element_text(size = 19))
fig <- fig + theme(panel.grid.minor = element_blank())

pdf(file='Law_W1.pdf',width=4,height=3.9)
fig
dev.off()

#MultiNess
W1_Lawp <- as.data.frame(fit_Lawp$U_hat[[1]][,c(3,1)])
names(W1_Lawp) <- c("W1","W2")
W1_Lawp$practice <- factor(practice)

fig <- ggplot() + geom_point(aes(x=W1,y=-W2,col=practice,shape=practice),data = W1_Lawp,size = 2)
fig <- fig + labs(x="1st component",y="2nd component") + scale_y_continuous(limits = c(-2.1,0)) + theme(axis.title = element_text(size=16))
fig <- fig + theme(legend.position = c(0.8,0.918)) + theme(legend.title=element_blank()) + theme(text = element_text(size = 19))
fig <- fig + theme(panel.grid.minor = element_blank())

pdf(file='Lawp_W1.pdf',width=4,height=3.9)
fig
dev.off()

##Figure 6
#SS-Refinement
W2_pca <- fit_Law$W_hat[[2]] %*% eigen(t(fit_Law$W_hat[[2]]) %*% fit_Law$W_hat[[2]])$vectors
W2_Law <- as.data.frame(W2_pca)
names(W2_Law) <- c("W1","W2","W3","W4")
W2_Law$status <- factor(status)

fig <- ggplot() + geom_point(aes(x=-W1,y=W2,col=status,shape=status),data = W2_Law,size = 2)
fig <- fig + labs(x="1st component",y="2nd component") + theme(axis.title = element_text(size=16))
fig <- fig + theme(legend.position = c(0.8,0.92)) + theme(legend.title=element_blank()) + theme(text = element_text(size = 19))
fig <- fig + theme(panel.grid.minor = element_blank())

pdf(file='Law_W2.pdf',width=4,height=3.9)
fig
dev.off()

#MultiNess
W2_Lawp <- as.data.frame(fit_Lawp$U_hat[[2]][,2:3])
names(W2_Lawp) <- c("W1","W2")
W2_Lawp$status <- factor(status)

fig <- ggplot() + geom_point(aes(x=-W1,y=W2,col=status,shape=status),data = W2_Lawp,size = 2)
fig <- fig + labs(x="1st component",y="2nd component") + theme(axis.title = element_text(size=16))
fig <- fig + theme(legend.position = c(0.8,0.92)) + theme(legend.title=element_blank()) + theme(text = element_text(size = 19))
fig <- fig + theme(panel.grid.minor = element_blank())

pdf(file='Lawp_W2.pdf',width=4,height=3.9)
fig
dev.off()
