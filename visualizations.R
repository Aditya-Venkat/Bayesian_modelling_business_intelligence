# 2. Model Evaluation -----------------------------------------------------
load("mcmc_samples.RData")

# Get true parameters (note true_w is 3x4)
true_C <- simulated_data$z
true_G <- simulated_data$G_true
true_pk <- do.call(rbind, simulated_data$pk_true)  # 4x15 (KxV)
true_pi <- array(unlist(simulated_data$pi_h_true), dim=c(3,15,15))  # HxVxV
true_w <- matrix(unlist(simulated_data$w_true), nrow=3, ncol=4)  # 3x4 (HxK)

# Calculate posterior means
post_mean_C <- apply(samples$C, 2, function(x) as.numeric(names(which.max(table(x)))))
post_mean_G <- apply(samples$G, 2, function(x) as.numeric(names(which.max(table(x)))))
post_mean_pk <- apply(samples$pk, c(2,3), mean)  # KxV
post_mean_pi <- apply(samples$pi, c(2,3,4), mean)  # HxVxV
post_mean_w <- apply(samples$nu, c(2,3), mean)  # Should be 3x4 (HxK)

# Verify dimensions
stopifnot(dim(true_w) == c(3,4))
stopifnot(dim(post_mean_w) == c(3,4))

# Evaluation metrics
metrics <- list(
  ARI_C = adjustedRandIndex(post_mean_C, true_C),
  ARI_G = adjustedRandIndex(post_mean_G, true_G),
  PK_RMSE = sqrt(mean((post_mean_pk - true_pk)^2)),
  PI_RMSE = sapply(1:3, function(h) sqrt(mean((post_mean_pi[h,,] - true_pi[h,,])^2))),
  W_RMSE = sqrt(mean((post_mean_w - true_w)^2))  # Now properly comparing 3x4 matrices
)

par(mfrow=c(2,2))
for(k in 1:4){
  acf(samples$pk[,k,1], main=paste("ACF: pk[",k,",1]"))
}
par(mfrow=c(1,1))
plot.new()



joint_dist <- table(
  True_Cluster=true_C,
  True_Scenario=true_G,
  Estimated_Cluster=post_mean_C,
  Estimated_Scenario=post_mean_G
)

# Flatten for visualization
df_joint <- as.data.frame.table(joint_dist)
ggplot(df_joint, aes(x=interaction(True_Cluster,True_Scenario),
                     y=interaction(Estimated_Cluster,Estimated_Scenario),
                     fill=Freq)) +
  geom_tile() +
  scale_fill_viridis_c() +
  labs(title="Joint Cluster-Scenario Recovery",
       x="True (Cluster.Scenario)", y="Estimated (Cluster.Scenario)") +
  theme(axis.text.x = element_text(angle=90))


