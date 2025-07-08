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

# Named list of grouping variables
groupings <- list(
  practice = practice,
  office = office,
  status = status
)

##################################################
#Figure S10: pairwise scatterplots for Z 
# Extract all components from V_hat (Shared)
U_hat_df <- as.data.frame(fit_Lawp$V_hat)

fig_dim = 11

# Loop through each grouping and create + save plot
for (group_name in names(groupings)) {
  fig <- plot_component_pairs(U_hat_df, groupings[[group_name]])
  print(fig)
  
  pdf(file = paste0("Lawp_Z_pair_", group_name, ".pdf"), width = fig_dim, height = fig_dim)
  print(fig)
  dev.off()
}


####################################################
#Figure S11: pairwise scatterplots for W1
# Extract all components from U_hat
U_hat_df <- as.data.frame(fit_Lawp$U_hat[[1]])

fig_dim=5
# Loop through each grouping and create + save plot
for (group_name in names(groupings)) {
  fig <- plot_component_pairs(U_hat_df, groupings[[group_name]])
  print(fig)
  
  pdf(file = paste0("Lawp_W1_pair_", group_name, ".pdf"), width = fig_dim, height = fig_dim)
  print(fig)
  dev.off()
}

####################################################
#Figure S12: pairwise scatterplots for W2  
# Extract all components from U_hat
kvis = all_dim[3]
U_hat_df <- as.data.frame(fit_Lawp$U_hat[[2]][, 1:kvis])

fig_dim = 6

# Loop through each grouping and create + save plot
for (group_name in names(groupings)) {
  fig <- plot_component_pairs(U_hat_df, groupings[[group_name]])
  print(fig)
  
  pdf(file = paste0("Lawp_W2_pair_", group_name, ".pdf"), width = fig_dim, height = fig_dim)
  print(fig)
  dev.off()
}
