library(MASS)
library(bayesm) # For rdirichlet()
library(pgdraw) # For Polya-Gamma sampling

# Robust log-sum-exp function
logSumExp <- function(x) {
  m <- max(x, na.rm = TRUE)
  if (!is.finite(m)) return(-Inf)
  m + log(sum(exp(x - m), na.rm = TRUE))
}

# Safe probability normalization
safeNormProbs <- function(log_probs) {
  log_probs <- log_probs - logSumExp(log_probs)
  probs <- exp(log_probs)
  probs[is.na(probs) | !is.finite(probs)] <- 0
  if (sum(probs) == 0) probs <- rep(1/length(probs), length(probs))
  probs / sum(probs)
}

# 1. Data Preparation ------------------------------------------------
load("simulated_data.RData")

# Convert y to counts per agency
n_agencies <- length(simulated_data$y)
n <- 200
V <- 15 # products
K <- 4 # Fixed number of clusters (as requested)
H <- 3 # Number of co-subscription scenarios
R <- 2 # Latent space dimension

# Convert to product counts
niv <- t(sapply(simulated_data$y, function(yi) table(factor(yi, levels=1:V))))

# Convert L to edge vectors (upper triangular)
edge_vectors <- t(sapply(simulated_data$L, function(adj) adj[upper.tri(adj)]))
n_edges <- ncol(edge_vectors) # V*(V-1)/2 = 105

# 2. Initialize Parameters -------------------------------------------
set.seed(123)

# Hyperparameters (from paper)
a1 <- 2.5
a2 <- 3.5
sigma2_l <- 10
mu_l <- rep(0, n_edges)
alpha <- rep(1, V) # Dirichlet prior for pk

# Cluster assignments (initialize randomly)
C <- sample(1:K, n_agencies, replace=TRUE)

# Monoproduct probabilities (Dirichlet)
pk <- matrix(NA, K, V)
for(k in 1:K) pk[k,] <- rdirichlet(alpha)

# Co-subscription scenario assignments
G <- sample(1:H, n_agencies, replace=TRUE)

# Mixing weights (cluster-specific)
nu <- matrix(1/H, H, K) # Uniform initial weights
for(k in 1:K) nu[,k] <- rdirichlet(rep(1, H))

# Latent space parameters
X <- array(rnorm(H*V*R, 0, 0.1), dim=c(H, V, R))
lambda <- matrix(1, H, R) # Shrinkage parameters
Z <- rnorm(n_edges, mu_l, sqrt(sigma2_l)) # Global baseline

# Precompute edge pairs for faster indexing
edge_pairs <- which(upper.tri(diag(V)), arr.ind=TRUE)

# 3. MCMC Settings --------------------------------------------------
n_iter <- 1e4
burnin <- 1000
thin <- 5
store_samples <- (n_iter - burnin) / thin
samples <- list(
  C = matrix(NA, store_samples, n_agencies),
  G = matrix(NA, store_samples, n_agencies),
  pk = array(NA, dim=c(store_samples, K, V)),
  pi = array(NA, dim=c(store_samples, H, V, V)),
  nu = array(NA, dim=c(store_samples, H, K))  # Added storage for nu
)
# 4. Gibbs Sampler --------------------------------------------------
for(iter in 1:n_iter) {
  
  # Step 1: Update pk (cluster product probabilities)
  for(k in 1:K) {
    cluster_agencies <- which(C == k)
    if(length(cluster_agencies) > 0) {
      counts <- colSums(niv[cluster_agencies,, drop=FALSE])
      pk[k,] <- rdirichlet(alpha + counts)
    } else {
      pk[k,] <- rdirichlet(alpha)
    }
  }
  
  # Precompute XLambdaXT for all scenarios
  XLambdaXT <- array(NA, dim=c(H, V, V))
  for(h in 1:H) {
    XLambdaXT[h,,] <- X[h,,] %*% diag(lambda[h,]) %*% t(X[h,,])
  }
  
  # Step 2: Update G (co-subscription scenario assignments)
  for(i in 1:n_agencies) {
    log_probs <- numeric(H)
    k <- C[i]
    
    for(h in 1:H) {
      # Get edge probabilities under scenario h
      logit_p <- Z + XLambdaXT[h,,][upper.tri(XLambdaXT[h,,])]
      p <- 1/(1 + exp(-logit_p)) # More stable than plogis for extreme values
      
      # Bernoulli likelihood with safeguards
      log_lik <- sum(dbinom(edge_vectors[i,], 1, p, log=TRUE))
      log_prior <- log(nu[h,k])
      
      # Handle numerical issues
      if(!is.finite(log_lik)) log_lik <- -1e10
      if(!is.finite(log_prior)) log_prior <- -1e10
      
      log_probs[h] <- log_prior + log_lik
    }
    
    # Normalize and sample
    probs <- safeNormProbs(log_probs)
    G[i] <- sample(1:H, 1, prob=probs)
  }
  
  # Step 3: Update nu (mixing weights)
  for(k in 1:K) {
    cluster_agencies <- which(C == k)
    counts <- tabulate(G[cluster_agencies], nbins=H)
    nu[,k] <- rdirichlet(rep(1/H, H) + counts)
  }
  
  # Step 4: Full Latent Space Updates
  
  # 4a. Update Polya-Gamma variables
  omega <- matrix(1, H, n_edges)
  
  for(h in 1:H) {
    agencies_h <- which(G == h)
    n_h <- length(agencies_h)
    
    if(n_h > 0) {
      for(l in 1:n_edges) {
        omega[h,l] <- pgdraw(n_h, Z[l] + XLambdaXT[h,,][upper.tri(XLambdaXT[h,,])][l])
      }
    }
  }
  
  # 4b. Update shared similarity vector Z
  for(l in 1:n_edges) {
    precision_Z <- 1/sigma2_l
    mu_Z <- mu_l[l]/sigma2_l
    
    for(h in 1:H) {
      agencies_h <- which(G == h)
      n_h <- length(agencies_h)
      
      if(n_h > 0) {
        X_term <- XLambdaXT[h,,][upper.tri(XLambdaXT[h,,])][l]
        y_sum <- sum(edge_vectors[agencies_h,l])
        
        precision_Z <- precision_Z + sum(omega[h,l])
        mu_Z <- mu_Z + (y_sum - n_h/2 - omega[h,l]*X_term)
      }
    }
    
    post_var <- 1/precision_Z
    post_mean <- post_var * mu_Z
    Z[l] <- rnorm(1, post_mean, sqrt(post_var))
  }
  
  # 4c. Update latent coordinates X
  for(h in 1:H) {
    agencies_h <- which(G == h)
    n_h <- length(agencies_h)
    
    if(n_h > 0) {
      Lambda_h <- diag(lambda[h,])
      
      for(v in 1:V) {
        # Find all edges involving product v
        edge_mask <- which(edge_pairs[,1] == v | edge_pairs[,2] == v)
        
        # Create design matrix
        other_nodes <- ifelse(edge_pairs[edge_mask,1] == v, 
                              edge_pairs[edge_mask,2], 
                              edge_pairs[edge_mask,1])
        
        X_design <- X[h, other_nodes,] * sqrt(lambda[h,])
        omega_v <- omega[h, edge_mask]
        y_v <- colMeans(edge_vectors[agencies_h, edge_mask, drop=FALSE])
        Z_v <- Z[edge_mask]
        
        # Posterior covariance
        XtOmegaX <- t(X_design) %*% (X_design * omega_v)
        post_cov <- solve(XtOmegaX + diag(1, R))
        
        # Posterior mean
        resid <- y_v - 0.5 - omega_v * Z_v
        post_mean <- post_cov %*% (t(X_design) %*% resid)
        
        # Sample new coordinates
        X[h,v,] <- mvrnorm(1, post_mean, post_cov)
      }
    }
  }
  
  # 4d. Update shrinkage parameters lambda
  for(h in 1:H) {
    for(r in 1:R) {
      shape <- ifelse(r == 1, a1 + V*R/2, a2 + V*(R - r + 1)/2)
      rate <- 1 + 0.5 * sum(X[h,,r]^2 * cumprod(lambda[h,])[r]/lambda[h,r])
      lambda[h,r] <- rgamma(1, shape, rate=rate)
    }
  }
  
  # Step 5: Update edge probabilities π
  pi <- array(0, dim=c(H, V, V))
  for(h in 1:H) {
    logit_p <- Z + XLambdaXT[h,,][upper.tri(XLambdaXT[h,,])]
    p <- 1/(1 + exp(-logit_p))
    
    for(l in 1:n_edges) {
      i <- edge_pairs[l,1]
      j <- edge_pairs[l,2]
      pi[h,i,j] <- p[l]
      pi[h,j,i] <- p[l]
    }
    diag(pi[h,,]) <- 0
  }
  
  # Step 6: Update Cluster Assignments C (fixed K version)
  for(i in 1:n_agencies) {
    log_probs <- numeric(K)
    
    for(k in 1:K) {
      # Calculate product count likelihood
      log_p_y <- sum(niv[i,] * log(pmax(pk[k,], 1e-10)))
      
      # Calculate co-subscription likelihood
      current_G <- G[i]
      log_p_G <- log(pmax(nu[current_G, k], 1e-300))
      
      # Include cluster size probability (Chinese Restaurant Process)
      n_minus_i <- sum(C[-i] == k)
      log_probs[k] <- log(n_minus_i + alpha/K) + log_p_y + log_p_G
    }
    
    # Normalize and sample
    log_probs <- log_probs + rnorm(K, 0, 0.1)
    probs <- safeNormProbs(log_probs)
    C[i] <- sample(1:K, 1, prob=probs)
  }
  
  # Store samples if past burnin and on thinning interval
  if(iter > burnin && iter %% thin == 0) {
    sample_idx <- (iter - burnin) / thin
    samples$C[sample_idx,] <- C
    samples$G[sample_idx,] <- G
    samples$pk[sample_idx,,] <- pk
    samples$pi[sample_idx,,,] <- pi
    samples$nu[sample_idx,,] <- nu
  }
  
  # Progress monitoring
  if(iter %% 50 == 0) {
    cat("Iteration", iter, 
        "Cluster counts:", tabulate(C, nbins=K),
        "Mixing weights:", round(nu[,1], 3), "...\n")  # Show weights for first cluster
  }
}


plot_cluster_trace <- function() {
  df <- data.frame(
    iteration = rep((burnin+1):n_iter, each=n_agencies),
    agency = rep(1:n_agencies, (n_iter-burnin)),
    cluster = factor(c(t(samples$C)))
  )
  
  ggplot(df, aes(x=iteration, y=agency, fill=cluster)) +
    geom_raster() +
    scale_fill_viridis_d() +
    labs(title="Cluster Assignment Trace",
         x="Iteration", y="Agency Index") +
    theme_minimal()
}



plot_pi_trace <- function() {
  plot(samples$pi[,1,1,2], type='l',
       main="Trace: pi[1,1,2] (Scenario 1, Products 1-2)",
       xlab="Iteration", ylab="Probability")
  abline(h=true_pi[1,1,2], col="red", lty=2)
}

plot_pi_trace()


save(samples, file= "mcmc_samples.RData")

