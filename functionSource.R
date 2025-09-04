####

## data generating mechanism ----

# variables that should actually be variable:
# N , autocorr effects, timepoints, N of x/y related confounders

DGM <- function(N = 20000, 
                autocorr_effects = 0.2, 
                timepoints = 100, 
                n_x_confounders = 5, 
                n_y_confounders = 5,
                beta_focal = 0.1,
                beta_covfocal = 0.1,
                beta_covnofocal = 0.01,
                sigma = 0.5,
                seed = 42) {
  
  # Load required libraries
  if (!require(mvtnorm)) {
    stop("Package 'mvtnorm' is required. Please install it using: install.packages('mvtnorm')")
  }
  
  # Set seed for reproducibility
  set.seed(seed)
  
  # Set parameters based on function inputs
  T <- timepoints
  c1 <- n_x_confounders
  c2 <- n_y_confounders
  alpha <- autocorr_effects
  
  # Helper functions for matrix operations (equivalents to irow/row functions)
  irow <- function(x) matrix(x, nrow = sqrt(length(x)), byrow = TRUE)
  row <- function(x) as.vector(t(x))
  
  # Create objects for the simulation
  ## total number of variables
  v <- c1 + c2 + 2
  
  ## A: Matrix of path coefficients
  A <- matrix(NA, v, v)
  diag(A) <- alpha
  A[lower.tri(A, diag = FALSE)] <- beta_covnofocal 
  A[upper.tri(A, diag = FALSE)] <- beta_covnofocal
  
  c1t <- 3:(2+c1)
  c2t <- (3+c1):v
  A[c1t, 1] <- beta_covfocal
  A[1, c1t] <- beta_covfocal
  A[c2t, 2] <- beta_covfocal
  A[2, c2t] <- beta_covfocal
  A[c1t, c1t][lower.tri(A[c1t, c1t], diag = FALSE)] <- beta_covfocal
  A[c1t, c1t][upper.tri(A[c1t, c1t], diag = FALSE)] <- beta_covfocal
  A[c2t, c2t][lower.tri(A[c2t, c2t], diag = FALSE)] <- beta_covfocal
  A[c2t, c2t][upper.tri(A[c2t, c2t], diag = FALSE)] <- beta_covfocal
  A[1,2] <- beta_focal
  A[2,1] <- beta_focal
  
  ## Check the max eigenvalue (to check convergence)
  max_eigenvalue <- max(eigen(A)$values)
  if (max_eigenvalue >= 1) {
    warning(paste("Maximum eigenvalue is", round(max_eigenvalue, 4), 
                  "- the system may not be stationary"))
  }
  
  ## Sigma: Var-cov matrix of time-specific residuals
  Sigma <- diag(rep(sigma, v))
  
  ## Sigma1: Stationary Var-cov matrix for the first time point
  Sigma1 <- irow(solve(diag(dim(A)[1]^2) - t(A) %x% t(A)) %*% row(Sigma))
  
  ## Mean1: Stationary means for first time point
  Mean1 <- matrix(0, nrow = dim(A)[1], ncol = 1)
  
  # Generate data
  ## Wide-format "bivariate" (X and Y) data df: N x 2*T
  df <- matrix(NA, nrow = N, ncol = 2*T)
  colnames(df) <- c(paste("x", 1:T, sep = ""), paste("y", 1:T, sep = ""))
  
  ## Generate the initial values for X and Y and store them in D
  D <- rmvnorm(N, mean = Mean1, Sigma1)
  
  ## Store the initial data from D to df
  df[, 1] <- D[, 1]
  df[, 1+T] <- D[, 2]
  
  ## Update D and df till time T
  for (i in 2:T) {
    D <- D %*% t(A) + rmvnorm(N, sigma = Sigma)
    df[, i] <- D[, 1]
    df[, i+T] <- D[, 2]
  }
  
  # Return results as a list
  return(list(
    data = df,
    parameters = list(
      N = N,
      timepoints = T,
      n_x_confounders = c1,
      n_y_confounders = c2,
      autocorr_effects = alpha,
      beta_focal = beta_focal,
      beta_covfocal = beta_covfocal,
      beta_covnofocal = beta_covnofocal,
      sigma = sigma,
      max_eigenvalue = max_eigenvalue
    ),
    coefficient_matrix = A
  ))
}

# -------------------------------------------------------------------------------
# Simple function to plot X and Y sequences side by side
plot_xy_sequences <- function(data, timepoints, n_participants = 5) {
  
  # Select random participants
  N <- nrow(data)
  participants <- sample(1:N, n_participants)
  
  # Set up plotting area
  par(mfrow = c(1, 2), mar = c(4, 4, 3, 1))
  
  # Plot X sequences
  plot(1:timepoints, data[participants[1], 1:timepoints], 
       type = "l", col = 1, ylim = range(data[participants, 1:timepoints]),
       main = "X Over Time", xlab = "Time", ylab = "X Value")
  
  for (i in 2:n_participants) {
    lines(1:timepoints, data[participants[i], 1:timepoints], col = i)
  }
  
  # Plot Y sequences  
  y_cols <- (timepoints + 1):(2 * timepoints)
  plot(1:timepoints, data[participants[1], y_cols], 
       type = "l", col = 1, ylim = range(data[participants, y_cols]),
       main = "Y Over Time", xlab = "Time", ylab = "Y Value")
  
  for (i in 2:n_participants) {
    lines(1:timepoints, data[participants[i], y_cols], col = i)
  }
  
  # Reset plotting parameters
  par(mfrow = c(1, 1))
}

# -------------------------------------------------------------------------------
# Create output structure with metadata

create_metadata <- function(){
  list(
    metadata = list(
      timestamp = Sys.time(),
      r_version = R.version.string,
      session_info = sessionInfo()
    ),
    results = data.frame()
    )
}

# -------------------------------------------------------------------------------
# Save results and clear memory function

save_clear <- function(results, filename_prefix = "sim_results"){
  timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
  filename <- paste0(filename_prefix, "_", timestamp, ".rds")
  
  saveRDS(results, file = filename)
  cat("Saved: ", filename, "\n")
  
  gc(verbose = F) # clears memory
  return(filename)
}

# -------------------------------------------------------------------------------
# generate RI-CLPM specification

generate_riclpm <- function(T) {
  
  # --- 1. Input Validation ---
  if (!is.numeric(T) || T <= 1 || T %% 1 != 0) {
    stop("T must be an integer greater than 1.")
  }
  
  # --- 2. Procedurally Generate Model Components ---
  
  # Random Intercepts (RIx and RIy)
  ri_x <- sprintf("RIx =~ %s", paste(sprintf("1*x%d", 1:T), collapse = " + "))
  ri_y <- sprintf("RIy =~ %s", paste(sprintf("1*y%d", 1:T), collapse = " + "))
  
  # Within-person component definitions (wx and wy)
  wx_defs <- paste(sprintf("wx%d =~ 1*x%d", 1:T, 1:T), collapse = "\n")
  wy_defs <- paste(sprintf("wy%d =~ 1*y%d", 1:T, 1:T), collapse = "\n")
  
  # Autoregressive and cross-lagged paths (from T2 to T)
  timepoints_reg <- 2:T
  lagged_timepoints <- 1:(T - 1)
  paths_wx <- paste(sprintf("wx%d ~ ax*wx%d + by*wy%d", timepoints_reg, lagged_timepoints, lagged_timepoints), collapse = "\n")
  paths_wy <- paste(sprintf("wy%d ~ bx*wx%d + ay*wy%d", timepoints_reg, lagged_timepoints, lagged_timepoints), collapse = "\n")
  
  # Constrained correlated residuals (from T2 to T)
  res_covs <- paste(sprintf("wx%d ~~ ur*wy%d", 2:T, 2:T), collapse = "\n")
  
  # Variances of within-person components (all timepoints)
  var_wx <- paste(sprintf("wx%d ~~ wx%d", 1:T, 1:T), collapse = "\n")
  var_wy <- paste(sprintf("wy%d ~~ wy%d", 1:T, 1:T), collapse = "\n")
  
  # Fix observed variances to zero
  zero_var_x <- paste(sprintf("x%d ~~ 0*x%d", 1:T, 1:T), collapse = "\n")
  zero_var_y <- paste(sprintf("y%d ~~ 0*y%d", 1:T, 1:T), collapse = "\n")
  
  
  # --- 3. Assemble Final Model String ---
  
  model_string <- paste(
    "# 1. Random Intercepts", ri_x, ri_y,
    "# 2. Within-person components", wx_defs, wy_defs,
    "# 3. Autoregressive and cross-lagged paths", paths_wx, paths_wy,
    "# 4. Covariances",
    "wx1 ~~ wy1 # Covariance at T1",
    res_covs,
    "# 5. (Co)variances of Random Intercepts",
    "RIx ~~ varRIx*RIx",
    "RIy ~~ varRIy*RIy",
    "RIx ~~ covRI*RIy",
    "# 6. (Residual) variances of within-person components",
    var_wx,
    var_wy,
    "# 7. Fix observed variances to zero",
    zero_var_x,
    zero_var_y,
    "# 8. Fix RI covariances with first state to zero",
    "wx1 ~~ 0*RIx",
    "wx1 ~~ 0*RIy",
    "wy1 ~~ 0*RIx",
    "wy1 ~~ 0*RIy",
    sep = "\n\n" # Separate sections with a blank line for readability
  )
  
  return(model_string)
}


