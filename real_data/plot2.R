library(ggplot2)
library(GGally)
load('GLM/real_data/results/lawyers.rda')

##########################################################################
#Table S3: Estimated latent dimensions by MultiNeSS+ for Lawyers Data.
##########################################################################
all_dim = c(fit_Lawp$d1, fit_Lawp$d2)
all_dim

##########################################################################
#Table S4: Correlations (rounded to one digit) between embeddings and nodewise features.
##########################################################################
compute_group_correlations <- function(U_hat_df, groupings) {
  # Ensure U_hat_df is a data frame
  U_hat_df <- as.data.frame(U_hat_df)
  
  # Compute correlation matrix
  cor_matrix <- sapply(groupings, function(group) {
    group_numeric <- as.numeric(factor(group))  # Convert factor to numeric
    sapply(U_hat_df, function(col) cor(col, group_numeric, method = "pearson"))
  })
  
  # Format result as data frame
  cor_df <- as.data.frame(cor_matrix)
  cor_df <- round(abs(cor_df), 4)  # Absolute value of correlation
  
  return(cor_df)
}

# Groupings
groupings <- list(
  office = office,
  practice = practice,
  status = status
)

all_correlation_results <- NULL

# Z
U_hat_df <- as.data.frame(fit_Lawp$V_hat)
correlation_results <- compute_group_correlations(U_hat_df, groupings)
rownames(correlation_results) <- paste0("Z_", seq_len(ncol(U_hat_df)))
all_correlation_results <- rbind(all_correlation_results,correlation_results)

# W1
U_hat_df <- as.data.frame(fit_Lawp$U_hat[[1]])
correlation_results <- compute_group_correlations(U_hat_df, groupings)
rownames(correlation_results) <- paste0("W1_", seq_len(ncol(U_hat_df)))
all_correlation_results <- rbind(all_correlation_results,correlation_results)

# W2
U_hat_df <- as.data.frame(fit_Lawp$U_hat[[2]])
correlation_results <- compute_group_correlations(U_hat_df, groupings)
rownames(correlation_results) <- paste0("W2_", seq_len(ncol(U_hat_df)))
all_correlation_results <- rbind(all_correlation_results,correlation_results)

# W3
U_hat_df <- as.data.frame(fit_Lawp$U_hat[[3]])
correlation_results <- compute_group_correlations(U_hat_df, groupings)
rownames(correlation_results) <- paste0("W3_", seq_len(ncol(U_hat_df)))
all_correlation_results <- rbind(all_correlation_results,correlation_results)

round(all_correlation_results,1)

######################################################################
#Figure S10-S12: Pairwise scatterplots 
######################################################################
plot_component_pairs <- function(U_hat_df, group_var, base_size = 16) {
  # Rename component columns
  comp_cols <- ncol(U_hat_df)
  colnames(U_hat_df) <- paste0("Comp", seq_len(comp_cols))
  
  # Add group variable as factor
  U_hat_df$Group <- factor(group_var)
  
  # Create ggpairs plot
  fig <- ggpairs(
    U_hat_df,
    columns = 1:comp_cols,
    mapping = aes(color = Group, shape = Group),
    diag = list(continuous = "blankDiag"),
    upper = "blank",
    lower = list(continuous = wrap("points", size = 2)),
    axisLabels = "none"
  )
  
  # Remove axis text and ticks from all panels
  for(i in 1:comp_cols){
    fig[i,i] <- fig[i,i] + theme_void()
  }
  
  # Add global theme
  fig <- fig + theme_bw(base_size = base_size) +
    theme(legend.position = "top", legend.title = element_blank())
  
  return(fig)
}

##Figure S10
U_hat_df <- as.data.frame(fit_Lawp$V_hat)
fig_dim = 11
fig <- plot_component_pairs(U_hat_df, office)
  
pdf(file = paste0("Lawp_Z_pair_office.pdf"), width = fig_dim, height = fig_dim)
fig
dev.off()

##Figure S11
U_hat_df <- as.data.frame(fit_Lawp$U_hat[[1]])
fig_dim = 5
fig <- plot_component_pairs(U_hat_df, practice)
 
pdf(file = paste0("Lawp_W1_pair_practice.pdf"), width = fig_dim, height = fig_dim)
fig
dev.off()

##Figure S12
U_hat_df <- as.data.frame(fit_Lawp$U_hat[[2]])
fig_dim = 6
fig <- plot_component_pairs(U_hat_df, status)
 
pdf(file = paste0("Lawp_W2_pair_status.pdf"), width = fig_dim, height = fig_dim)
fig
dev.off()
