library(MASS)
set.seed(123)

# Parameters
n_agencies <- 200
K <- 4
V <- 15
S <- 500  # monoproduct customers per agency
H <- 3

# Assign clusters
cluster_labels <- rep(1:K, each = n_agencies / K)

# Define monoproduct PMFs (pk0) as described in the paper
pk_list <- vector("list", K)

# Base distribution for clusters 1 and 2 (all products equal except 1 and 9)
base_p <- rep(0.02, V)
base_p[c(1, 9)] <- c(0.3, 0.1)  # product 1 more probable than 9 in cluster 1
base_p <- base_p / sum(base_p)

pk_list[[1]] <- base_p
pk_list[[2]] <- base_p
pk_list[[2]][1] <- base_p[9]  # swap probabilities for products 1 and 9
pk_list[[2]][9] <- base_p[1]

# Base distribution for clusters 3 and 4 (all products equal except 3 and 7)
base_p_34 <- rep(0.02, V)
base_p_34[c(3, 7)] <- c(0.3, 0.1)  # product 3 more probable than 7 in cluster 3
base_p_34 <- base_p_34 / sum(base_p_34)

pk_list[[3]] <- base_p_34
pk_list[[4]] <- base_p_34
pk_list[[4]][3] <- base_p_34[7]  # swap probabilities for products 3 and 7
pk_list[[4]][7] <- base_p_34[3]

# Simulate monoproduct choices
y_list <- lapply(1:n_agencies, function(i) {
  cluster <- cluster_labels[i]
  sample(1:V, S, replace = TRUE, prob = pk_list[[cluster]])
})

# Construct H = 3 co-subscription probability matrices π^h as described
generate_co_prob <- function(type) {
  mat <- matrix(0.01, V, V)  # baseline low probability
  
  if (type == 1) {
    # Scenario 1: one dense community among 10 products
    dense_products <- 1:10
    mat[dense_products, dense_products] <- 0.7
  } else if (type == 2) {
    # Scenario 2: four hub products (2,5,8,12)
    hub_products <- c(2, 5, 8, 12)
    mat[hub_products, ] <- 0.6
    mat[, hub_products] <- 0.6
    diag(mat) <- 0  # no self-loops
  } else if (type == 3) {
    # Scenario 3: similar to 2 but without product 4 as hub
    hub_products <- c(2, 5, 8, 12)  # same as type 2
    mat[hub_products, ] <- 0.6
    mat[, hub_products] <- 0.6
    mat[4, ] <- 0.01  # product 4 held out
    mat[, 4] <- 0.01
    diag(mat) <- 0
  }
  
  # Ensure symmetry and no self-loops
  mat <- (mat + t(mat)) / 2
  diag(mat) <- 0
  return(mat)
}

pi_h_list <- lapply(1:H, generate_co_prob)

# Cluster-specific mixing weights as specified in the paper
w_k <- list(
  c(0.9, 0.05, 0.05),  # cluster 1: mostly scenario 1
  c(0.9, 0.05, 0.05),  # cluster 2: same as cluster 1
  c(0.05, 0.9, 0.05),  # cluster 3: mostly scenario 2
  c(0.05, 0.05, 0.9)   # cluster 4: mostly scenario 3
)

# Simulate co-subscription adjacency matrices (binary, symmetric)
simulate_network <- function(w, pi_h) {
  # First choose which scenario to use
  h <- sample(1:H, 1, prob = w)
  probs <- pi_h[[h]]
  
  # Generate symmetric adjacency matrix
  adj <- matrix(0, V, V)
  for (i in 1:(V-1)) {
    for (j in (i+1):V) {
      adj[i, j] <- rbinom(1, 1, probs[i, j])
      adj[j, i] <- adj[i, j]
    }
  }
  return(adj)
}

L_list <- lapply(1:n_agencies, function(i) {
  cluster <- cluster_labels[i]
  simulate_network(w_k[[cluster]], pi_h_list)
})



# Output is a list of monoproduct vectors and co-subscription networks
simulated_data <- list(
  y = y_list,          # monoproduct choices
  L = L_list,          # co-subscription adjacency matrices
  z = cluster_labels,  # true cluster assignment
  pk_true = pk_list,   # true monoproduct probabilities
  pi_h_true = pi_h_list, # true co-subscription scenarios
  w_true = w_k         # true mixing weights
)
# In your simulation code, right before saving:
# Store the true scenario assignments (G) for each agency
true_G <- sapply(1:n_agencies, function(i) {
  cluster <- simulated_data$z[i]
  sample(1:H, 1, prob=w_k[[cluster]])
})

# Add to the saved data
simulated_data <- c(simulated_data, list(G_true = true_G))

# Then save as before
save(simulated_data, file = "simulated_data.RData")