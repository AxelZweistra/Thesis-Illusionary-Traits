### Functions for illusory trait thesis 

## data generating mechanism - updated to be able to set different CL params

DGM <- function(N = 20000, 
                autocorr_effects = 0.2, 
                timepoints = 100, 
                n_x_confounders = 5, 
                n_y_confounders = 5,
                beta_focal = 0.1,    # Maintained for backward compatibility
                beta_xy = NULL,      # New: X -> Y
                beta_yx = NULL,      # New: Y -> X
                beta_covfocal = 0.1,
                beta_covnofocal = 0.01,
                sigma = 0.5,
                ri_var_x = 0, 
                ri_var_y = 0,    
                ri_cor = 0, 
                seed = 427) 
{
  # Backward Compatibility Logic (update to allow separate beta specifications:
  # If the specific beta_xy/yx are NULL, they default to beta_focal
  if (is.null(beta_xy)) beta_xy <- beta_focal
  if (is.null(beta_yx)) beta_yx <- beta_focal
  
  set.seed(seed)
  
  # rename to match original by Bailey
  T <- timepoints
  c1 <- n_x_confounders
  c2 <- n_y_confounders
  alpha <- autocorr_effects
  
  # should already be in the enviroment, failstop
  irow <- function(x) matrix(x, nrow = sqrt(length(x)), byrow = TRUE)
  row <- function(x) as.vector(t(x))
  
  include_time_varying_confounding <- (c1 > 0 || c2 > 0)
  include_time_invariant_confounding <- (ri_var_x > 0 || ri_var_y > 0)
  
  # time-varying component only needs to get called if necessary
  if (include_time_varying_confounding) {
    v <- c1 + c2 + 2
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
    
    # Asymmetrical paths:
    # In the original matrix A, row 2 col 1 is X -> Y; row 1 col 2 is Y -> X
    # Updated path for focal variables
    A[2, 1] <- beta_xy # X -> Y
    A[1, 2] <- beta_yx # Y -> X
    
    # Use Mod() to handle potential complex values from non-symmetric A
    max_eigenvalue <- max(Mod(eigen(A)$values))
    
    if (max_eigenvalue >= 1) {
      warning(paste("Maximum eigenvalue magnitude is", round(max_eigenvalue, 4),
                    "- the system may not be stationary"))
    }
    
    Sigma <- diag(rep(sigma, v))
    Sigma1 <- irow(solve(diag(dim(A)[1]^2) - t(A) %x% t(A)) %*% row(Sigma))
    Mean1 <- matrix(0, nrow = dim(A)[1], ncol = 1)
    
  } else {
    v <- 2

    A <- matrix(c(alpha, beta_yx, beta_xy, alpha), nrow = 2, byrow = TRUE)
    
    # warning for non-stationarity
    max_eigenvalue <- max(eigen(A)$values)
    if (max_eigenvalue >= 1) {
      warning(paste("Maximum eigenvalue is", round(max_eigenvalue, 4),
                    "- the system may not be stationary"))
    }
    
    Sigma <- diag(rep(sigma, v))
    Sigma1 <- irow(solve(diag(dim(A)[1]^2) - t(A) %x% t(A)) %*% row(Sigma))
    Mean1 <- matrix(0, nrow = dim(A)[1], ncol = 1)
  }
  
  df_within <- matrix(NA, nrow = N, ncol = 2*T)
  D <- rmvnorm(N, mean = Mean1, Sigma1)
  df_within[, 1] <- D[, 1]
  df_within[, 1+T] <- D[, 2]
  
  for (i in 2:T) {
    D <- D %*% t(A) + rmvnorm(N, sigma = Sigma)
    df_within[, i] <- D[, 1]
    df_within[, i+T] <- D[, 2]
  }
  
  # only necessary if time-invariant confounding is present
  if (include_time_invariant_confounding) {
    ri_cov_val <- ri_cor * sqrt(ri_var_x * ri_var_y)
    ri_sigma_mat <- matrix(c(ri_var_x, ri_cov_val, ri_cov_val, ri_var_y), nrow = 2)
    RIs <- rmvnorm(N, mean = c(0, 0), sigma = ri_sigma_mat)
    colnames(RIs) <- c("RI_x", "RI_y")
    df_final <- df_within
    colnames(df_final) <- c(paste("x", 1:T, sep = ""), paste("y", 1:T, sep = ""))
    df_final[, 1:T] <- df_final[, 1:T] + RIs[, 1]
    df_final[, (T+1):(2*T)] <- df_final[, (T+1):(2*T)] + RIs[, 2]
  } else {
    df_final <- df_within
    colnames(df_final) <- c(paste("x", 1:T, sep = ""), paste("y", 1:T, sep = ""))
    RIs <- NULL
  }
  
  result <- list(
    data = df_final,
    parameters = list(
      N = N,
      timepoints = T,
      n_x_confounders = c1,
      n_y_confounders = c2,
      autocorr_effects = alpha,
      beta_focal = beta_focal, 
      beta_xy = beta_xy,       
      beta_yx = beta_yx,
      beta_covfocal = beta_covfocal,
      beta_covnofocal = beta_covnofocal,
      sigma = sigma,
      StatCOVmatrix = Sigma1,
      max_eigenvalue = max_eigenvalue,
      ri_var_x = ri_var_x,
      ri_var_y = ri_var_y,
      ri_cor = ri_cor,
      data_type = paste0(
        ifelse(include_time_varying_confounding, "time-varying confounding", "no time-varying confounding"),
        " + ",
        ifelse(include_time_invariant_confounding, "time-invariant confounding", "no time-invariant confounding")
      )
    ),
    coefficient_matrix = A
  )
  
  if (include_time_invariant_confounding) {
    result$components <- list(within = df_within, between = RIs)
  }
  
  return(result)
}

# -------------------------------------------------------------------------------
# Simple function to plot X and Y sequences side by side
plot_xy_sequences <- function(data, timepoints, n_participants = 5) {
  
  # Sample some random participants
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
# Create output structure with metadata -  helper

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
# Save results and clear memory function - helper

save_clear <- function(results, filename_prefix = "sim_results"){
  timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
  filename <- paste0(filename_prefix, "_", timestamp, ".rds")
  
  # save as RDS file
  saveRDS(results, file = filename)
  cat("Saved: ", filename, "\n")
  
  gc(verbose = F) # clears memory
  return(filename)
}

# -------------------------------------------------------------------------------

# Model specifications

# Generate CLPM model specification
# only parameter is T: number of timepoints

generate_clpm <- function(T){
  
  # Input validation
  if (!is.numeric(T) || T <= 1 || T %% 1 != 0) {
    stop("T must be an integer greater or equal to 2.")
  }
  
  # Generate model components
  
  # Autoregressive and cross-lagged paths (from T2 to T)
  if (T > 1) {
    timepoints_reg <- 2:T
    lagged_timepoints <- 1:(T - 1)
    paths_x <- paste(sprintf("x%d ~ ax*x%d + by*y%d", 
                             timepoints_reg, lagged_timepoints, lagged_timepoints), 
                     collapse = "\n")
    paths_y <- paste(sprintf("y%d ~ bx*x%d + ay*y%d", 
                             timepoints_reg, lagged_timepoints, lagged_timepoints), 
                     collapse = "\n")
  } else {
    paths_x <- ""
    paths_y <- ""
  }
  
  # Variances of observed variables
  # First wave: free variance (exogenous starting values)
  var_x1 <- "x1 ~~ x1"
  var_y1 <- "y1 ~~ y1"
  
  # Subsequent waves: residual/innovation variances (constrained equal)
  if (T > 1) {
    var_x_rest <- paste(sprintf("x%d ~~ ivx*x%d", 2:T, 2:T), collapse = "\n")
    var_y_rest <- paste(sprintf("y%d ~~ ivy*y%d", 2:T, 2:T), collapse = "\n")
  } else {
    var_x_rest <- ""
    var_y_rest <- ""
  }
  
  # Covariance at first wave (exogenous)
  cov_xy1 <- "x1 ~~ y1 # Covariance at T1"
  
  # Correlated residuals at subsequent waves (constrained equal)
  if (T > 1) {
    res_covs <- paste(sprintf("x%d ~~ ur*y%d", 2:T, 2:T), collapse = "\n")
  } else {
    res_covs <- ""
  }
  
  # Assemble final model string
  
  components <- c(
    "# 1. Autoregressive and cross-lagged paths",
    paths_x,
    paths_y,
    "\n# 2. Covariances",
    cov_xy1,
    res_covs,
    "\n# 3. Variances of observed variables",
    "# 3a. First wave variances (exogenous starting values)",
    var_x1,
    var_y1,
    "# 3b. Subsequent wave innovation variances (constrained equal)",
    var_x_rest,
    var_y_rest
  )
  
  # Add together, skip empty strings and add linebreaks (for legibility)
  model_string <- paste(components[components != ""], collapse = "\n")
  
  return(model_string)
}
# ------------------------------------------------------------------------------

# Generate RI-CLPM model specification
# only parameter is T: number of timepoints

generate_riclpm <- function(T) {
  
  # Input validation
  if (!is.numeric(T) || T <= 2 || T %% 1 != 0) {
    stop("T must be an integer greater or equal to 3.")
  }
  
  # Generate model components
  
  # Random intercepts (RIx and RIy)
  ri_x <- sprintf("RIx =~ %s", paste(sprintf("1*x%d", 1:T), collapse = " + "))
  ri_y <- sprintf("RIy =~ %s", paste(sprintf("1*y%d", 1:T), collapse = " + "))
  
  # Within-person component definitions (wx and wy)
  wx_defs <- paste(sprintf("wx%d =~ 1*x%d", 1:T, 1:T), collapse = "\n")
  wy_defs <- paste(sprintf("wy%d =~ 1*y%d", 1:T, 1:T), collapse = "\n")
  
  # Autoregressive and cross-lagged paths (from T2 to T)
  if (T > 1) {
    timepoints_reg <- 2:T
    lagged_timepoints <- 1:(T - 1)
    paths_wx <- paste(sprintf("wx%d ~ ax*wx%d + by*wy%d", 
                              timepoints_reg, lagged_timepoints, lagged_timepoints), 
                      collapse = "\n")
    paths_wy <- paste(sprintf("wy%d ~ bx*wx%d + ay*wy%d", 
                              timepoints_reg, lagged_timepoints, lagged_timepoints), 
                      collapse = "\n")
  } else {
    paths_wx <- ""
    paths_wy <- ""
  }
  
  # Covariances
  # First wave: covariance (exogenous)
  # Subsequent waves: constrained correlated residuals (endogenous)
  if (T > 1) {
    res_covs <- paste(sprintf("wx%d ~~ ur*wy%d", 2:T, 2:T), collapse = "\n")
  } else {
    res_covs <- ""
  }
  
  # Variances of within-person components
  # First wave: variance (exogenous, no predictors)
  var_wx1 <- "wx1 ~~ wx1"
  var_wy1 <- "wy1 ~~ wy1"
  
  # Subsequent waves: residual variances (endogenous, have predictors)
  if (T > 1) {
    var_wx_rest <- paste(sprintf("wx%d ~~ wx%d", 2:T, 2:T), collapse = "\n")
    var_wy_rest <- paste(sprintf("wy%d ~~ wy%d", 2:T, 2:T), collapse = "\n")
  } else {
    var_wx_rest <- ""
    var_wy_rest <- ""
  }
  
  # Fix observed variances to zero
  zero_var_x <- paste(sprintf("x%d ~~ 0*x%d", 1:T, 1:T), collapse = "\n")
  zero_var_y <- paste(sprintf("y%d ~~ 0*y%d", 1:T, 1:T), collapse = "\n")
  
  # Assemble final model string
  
  components <- c(
    "# 1. Random Intercepts", 
    ri_x, 
    ri_y,
    "\n# 2. Within-person components", 
    wx_defs, 
    wy_defs,
    "\n# 3. Autoregressive and cross-lagged paths", 
    paths_wx, 
    paths_wy,
    "\n# 4. Covariances",
    "wx1 ~~ wy1 # Covariance at T1",
    res_covs,
    "\n# 5. (Co)variances of Random Intercepts",
    "RIx ~~ varRIx*RIx",
    "RIy ~~ varRIy*RIy",
    "RIx ~~ covRI*RIy",
    "\n# 6. Variances and residual variances of within-person components",
    "# 6a. First wave variances (exogenous)",
    var_wx1,
    var_wy1,
    "# 6b. Subsequent wave residual variances (endogenous)",
    var_wx_rest,
    var_wy_rest,
    "\n# 7. Fix observed variances to zero",
    zero_var_x,
    zero_var_y,
    "\n# 8. Fix RI covariances with first state to zero",
    "wx1 ~~ 0*RIx", 
    "wx1 ~~ 0*RIy", 
    "wy1 ~~ 0*RIx", 
    "wy1 ~~ 0*RIy"
  )
  
  # Add together, skip empty strings and add linebreaks
  model_string <- paste(components[components != ""], collapse = "\n")
  
  return(model_string)
}

# ------------------------------------------------------------------------------

# Generate SRI-CLPM model specification
# T: total number of timepoints, timespan: timepoints on which RI's are estimated
# across_seg_cov: has to be "zero", "free" or "toeplitz" 

generate_segmented_riclpm <- function(T, timespan = 3, across_seg_cov = "toeplitz") {
  
  # Input Validation
  if (!is.numeric(T) || T <= 2 || T %% 1 != 0) {
    stop("T must be an integer greater or equal to 3.")
  }
  
  if (!is.numeric(timespan) || timespan <= 0 || timespan %% 1 != 0) {
    stop("timespan must be a positive integer.")
  }
  
  if (timespan > T) {
    stop("timespan cannot be greater than T.")
  }
  
  if (!across_seg_cov %in% c("zero", "free", "toeplitz")) {
    stop("across_seg_cov must be either 'zero', 'free', or 'toeplitz'.")
  }
  
  # Determine Number of Segments
  n_segments <- ceiling(T / timespan)
  
  # Create a mapping of timepoints to segments
  timepoint_to_segment <- rep(1:n_segments, each = timespan)[1:T]
  
  # Generate Random Intercepts for Each Segment

  ri_x_definitions <- c()
  ri_y_definitions <- c()
  
  for (seg in 1:n_segments) {
    # Find which timepoints belong to this segment
    tp_in_seg <- which(timepoint_to_segment == seg)
    
    # Define RIx for this segment
    ri_x_seg <- sprintf("RIx%d =~ %s", seg, 
                        paste(sprintf("1*x%d", tp_in_seg), collapse = " + "))
    
    # Define RIy for this segment
    ri_y_seg <- sprintf("RIy%d =~ %s", seg, 
                        paste(sprintf("1*y%d", tp_in_seg), collapse = " + "))
    
    ri_x_definitions <- c(ri_x_definitions, ri_x_seg)
    ri_y_definitions <- c(ri_y_definitions, ri_y_seg)
  }
  
  # Combine: all X's first, then all Y's
  ri_definitions <- c(ri_x_definitions, ri_y_definitions)
  
  # Within-person Component Definitions
  wx_defs <- paste(sprintf("wx%d =~ 1*x%d", 1:T, 1:T), collapse = "\n")
  wy_defs <- paste(sprintf("wy%d =~ 1*y%d", 1:T, 1:T), collapse = "\n")
  
  # Autoregressive and Cross-lagged Paths
  if (T > 1) {
    timepoints_reg <- 2:T
    lagged_timepoints <- 1:(T - 1)
    paths_wx <- paste(sprintf("wx%d ~ ax*wx%d + by*wy%d", 
                              timepoints_reg, lagged_timepoints, lagged_timepoints), 
                      collapse = "\n")
    paths_wy <- paste(sprintf("wy%d ~ bx*wx%d + ay*wy%d", 
                              timepoints_reg, lagged_timepoints, lagged_timepoints), 
                      collapse = "\n")
  } else {
    paths_wx <- ""
    paths_wy <- ""
  }
  
  # Covariances
  # First wave: covariance (exogenous)
  # Subsequent waves: constrained correlated residuals (endogenous)
  if (T > 1) {
    res_covs <- paste(sprintf("wx%d ~~ ur*wy%d", 2:T, 2:T), collapse = "\n")
  } else {
    res_covs <- ""
  }
  
  # Variances of Within-person Components
  # First wave: variance (exogenous, no predictors)
  var_wx1 <- "wx1 ~~ wx1"
  var_wy1 <- "wy1 ~~ wy1"
  
  # Subsequent waves: residual variances (endogenous, have predictors)
  if (T > 1) {
    var_wx_rest <- paste(sprintf("wx%d ~~ wx%d", 2:T, 2:T), collapse = "\n")
    var_wy_rest <- paste(sprintf("wy%d ~~ wy%d", 2:T, 2:T), collapse = "\n")
  } else {
    var_wx_rest <- ""
    var_wy_rest <- ""
  }
  
  # Fix Observed Variances to Zero 
  zero_var_x <- paste(sprintf("x%d ~~ 0*x%d", 1:T, 1:T), collapse = "\n")
  zero_var_y <- paste(sprintf("y%d ~~ 0*y%d", 1:T, 1:T), collapse = "\n")
  
  # Random Intercept Variances and Covariances
  
  # Within-segment variances and within-timespan covariances (RIx with RIy at same segment)
  ri_variances_and_within_seg_covs <- c()
  
  # First: All RIx variances - CONSTRAINED EQUAL
  for (seg in 1:n_segments) {
    ri_variances_and_within_seg_covs <- c(ri_variances_and_within_seg_covs,
                                          sprintf("RIx%d ~~ varRIx*RIx%d", seg, seg))
  }
  
  # Second: All RIy variances - CONSTRAINED EQUAL
  for (seg in 1:n_segments) {
    ri_variances_and_within_seg_covs <- c(ri_variances_and_within_seg_covs,
                                          sprintf("RIy%d ~~ varRIy*RIy%d", seg, seg))
  }
  
  # Third: Within-segment RIx-RIy covariances - CONSTRAINED EQUAL ACROSS SEGMENTS
  for (seg in 1:n_segments) {
    ri_variances_and_within_seg_covs <- c(ri_variances_and_within_seg_covs,
                                          sprintf("RIx%d ~~ covRI*RIy%d", seg, seg))
  }
  
  # Across-segment covariances
  ri_across_seg_covs <- c()
  across_seg_label <- ""
  
  if (n_segments > 1) {
    if (across_seg_cov == "zero") {
      # Fixed to zero: independent confounding across segments
      across_seg_label <- "(Fixed to Zero)"
      for (seg1 in 1:(n_segments - 1)) {
        for (seg2 in (seg1 + 1):n_segments) {
          ri_across_seg_covs <- c(ri_across_seg_covs,
                                  sprintf("RIx%d ~~ 0*RIx%d", seg1, seg2),
                                  sprintf("RIy%d ~~ 0*RIy%d", seg1, seg2),
                                  sprintf("RIx%d ~~ 0*RIy%d", seg1, seg2),
                                  sprintf("RIy%d ~~ 0*RIx%d", seg1, seg2))
        }
      }
      
    } else if (across_seg_cov == "free") {
      # Freely estimate all covariances with unique labels
      across_seg_label <- "(Freely Estimated)"
      for (seg1 in 1:(n_segments - 1)) {
        for (seg2 in (seg1 + 1):n_segments) {
          ri_across_seg_covs <- c(ri_across_seg_covs,
                                  sprintf("RIx%d ~~ covRIx_%d_%d*RIx%d", seg1, seg1, seg2, seg2),
                                  sprintf("RIy%d ~~ covRIy_%d_%d*RIy%d", seg1, seg1, seg2, seg2),
                                  sprintf("RIx%d ~~ covRIxy_%d_%d*RIy%d", seg1, seg1, seg2, seg2),
                                  sprintf("RIy%d ~~ covRIyx_%d_%d*RIx%d", seg1, seg1, seg2, seg2))
        }
      }
      
    } else if (across_seg_cov == "toeplitz") {
      # Block Toeplitz structure: X->X, Y->Y, X->Y, Y->X each have their own lag parameters
      across_seg_label <- "(Block Toeplitz Structure)"
      
      max_lag <- n_segments - 1
      
      # Within-variable Toeplitz (X->X and Y->Y)
      for (lag in 1:max_lag) {
        for (seg1 in 1:(n_segments - lag)) {
          seg2 <- seg1 + lag
          
          ri_across_seg_covs <- c(ri_across_seg_covs,
                                  sprintf("RIx%d ~~ covRIx_lag%d*RIx%d", seg1, lag, seg2))
          ri_across_seg_covs <- c(ri_across_seg_covs,
                                  sprintf("RIy%d ~~ covRIy_lag%d*RIy%d", seg1, lag, seg2))
        }
      }
      
      # Cross-variable Toeplitz (X->Y and Y->X separate)
      for (lag in 1:max_lag) {
        for (seg1 in 1:(n_segments - lag)) {
          seg2 <- seg1 + lag
          
          # X -> Y covariances (all same-lag X->Y constrained equal)
          ri_across_seg_covs <- c(ri_across_seg_covs,
                                  sprintf("RIx%d ~~ covRIxy_lag%d*RIy%d", seg1, lag, seg2))
          
          # Y -> X covariances (all same-lag Y->X constrained equal, but different from X->Y)
          ri_across_seg_covs <- c(ri_across_seg_covs,
                                  sprintf("RIy%d ~~ covRIyx_lag%d*RIx%d", seg1, lag, seg2))
        }
      }
    }
  }
  
  # Fix RI Covariances with First State to Zero
  ri_first_state_zero <- c()
  for (seg in 1:n_segments) {
    ri_first_state_zero <- c(ri_first_state_zero,
                             sprintf("wx1 ~~ 0*RIx%d", seg),
                             sprintf("wx1 ~~ 0*RIy%d", seg),
                             sprintf("wy1 ~~ 0*RIx%d", seg),
                             sprintf("wy1 ~~ 0*RIy%d", seg))
  }
  
  # Assemble Final Model String
  
  components <- c(
    "# 1. Random Intercepts (Segmented) - X then Y",
    paste(ri_definitions, collapse = "\n"),
    "\n# 2. Within-person components",
    wx_defs,
    wy_defs,
    "\n# 3. Autoregressive and cross-lagged paths",
    paths_wx,
    paths_wy,
    "\n# 4. Covariances",
    "wx1 ~~ wy1 # Covariance at T1",
    res_covs,
    "\n# 5. (Co)variances of Random Intercepts",
    "# 5a. Variances (constrained equal) and within-segment covariances (constrained equal)",
    paste(ri_variances_and_within_seg_covs, collapse = "\n"),
    sprintf("\n# 5b. Across-segment covariances %s", 
            ifelse(n_segments > 1, across_seg_label, "")),
    if (length(ri_across_seg_covs) > 0) paste(ri_across_seg_covs, collapse = "\n") else "# (No across-segment covariances for single segment)",
    "\n# 6. Variances and residual variances of within-person components",
    "# 6a. First wave variances (exogenous)",
    var_wx1,
    var_wy1,
    "# 6b. Subsequent wave residual variances (endogenous)",
    var_wx_rest,
    var_wy_rest,
    "\n# 7. Fix observed variances to zero",
    zero_var_x,
    zero_var_y,
    "\n# 8. Fix RI covariances with first state to zero",
    paste(ri_first_state_zero, collapse = "\n")
  )
  
  components <- components[components != ""]
  model_string <- paste(components, collapse = "\n")
  
  return(model_string)
}

# Generate SRI-CLPM with Global RIs model specification
# T: total number of timepoints, timespan: timepoints on which RI's are estimated
# across_seg_cov: has to be "zero", "free" or "toeplitz" 

generate_segmented_riclpm_global <- function(T, timespan = 3, across_seg_cov = "toeplitz") {
  
  #  Input Validation
  if (!is.numeric(T) || T <= 2 || T %% 1 != 0) {
    stop("T must be an integer greater or equal to 3.")
  }
  
  if (!is.numeric(timespan) || timespan <= 0 || timespan %% 1 != 0) {
    stop("timespan must be a positive integer.")
  }
  
  if (timespan > T) {
    stop("timespan cannot be greater than T.")
  }
  
  if (!across_seg_cov %in% c("zero", "free", "toeplitz")) {
    stop("across_seg_cov must be either 'zero', 'free', or 'toeplitz'.")
  }
  
  #  Determine Number of Segments
  n_segments <- ceiling(T / timespan)
  
  # Create a mapping of timepoints to segments
  timepoint_to_segment <- rep(1:n_segments, each = timespan)[1:T]
  
  #  Generate Global Random Intercepts (on ALL timepoints)
  global_ri_x <- sprintf("RIgx =~ %s", paste(sprintf("1*x%d", 1:T), collapse = " + "))
  global_ri_y <- sprintf("RIgy =~ %s", paste(sprintf("1*y%d", 1:T), collapse = " + "))
  
  #  Generate Segmented Random Intercepts

  ri_x_definitions <- c()
  ri_y_definitions <- c()
  
  for (seg in 1:n_segments) {
    # Find which timepoints belong to this segment
    tp_in_seg <- which(timepoint_to_segment == seg)
    
    # Define RIx for this segment
    ri_x_seg <- sprintf("RIx%d =~ %s", seg, 
                        paste(sprintf("1*x%d", tp_in_seg), collapse = " + "))
    
    # Define RIy for this segment
    ri_y_seg <- sprintf("RIy%d =~ %s", seg, 
                        paste(sprintf("1*y%d", tp_in_seg), collapse = " + "))
    
    ri_x_definitions <- c(ri_x_definitions, ri_x_seg)
    ri_y_definitions <- c(ri_y_definitions, ri_y_seg)
  }
  
  # Combine: Global RIs first, then all segmented X's, then all segmented Y's
  ri_definitions <- c(global_ri_x, global_ri_y, ri_x_definitions, ri_y_definitions)
  
  #  Within-person Component Definitions
  wx_defs <- paste(sprintf("wx%d =~ 1*x%d", 1:T, 1:T), collapse = "\n")
  wy_defs <- paste(sprintf("wy%d =~ 1*y%d", 1:T, 1:T), collapse = "\n")
  
  #  Autoregressive and Cross-lagged Paths
  if (T > 1) {
    timepoints_reg <- 2:T
    lagged_timepoints <- 1:(T - 1)
    paths_wx <- paste(sprintf("wx%d ~ ax*wx%d + by*wy%d", 
                              timepoints_reg, lagged_timepoints, lagged_timepoints), 
                      collapse = "\n")
    paths_wy <- paste(sprintf("wy%d ~ bx*wx%d + ay*wy%d", 
                              timepoints_reg, lagged_timepoints, lagged_timepoints), 
                      collapse = "\n")
  } else {
    paths_wx <- ""
    paths_wy <- ""
  }
  
  # Covariances
  # First wave: covariance (exogenous)
  # Subsequent waves: constrained correlated residuals (endogenous)
  if (T > 1) {
    res_covs <- paste(sprintf("wx%d ~~ ur*wy%d", 2:T, 2:T), collapse = "\n")
  } else {
    res_covs <- ""
  }
  
  # Variances of Within-person Components
  # First wave: variance (exogenous, no predictors)
  var_wx1 <- "wx1 ~~ wx1"
  var_wy1 <- "wy1 ~~ wy1"
  
  # Subsequent waves: residual variances (endogenous, have predictors)
  if (T > 1) {
    var_wx_rest <- paste(sprintf("wx%d ~~ wx%d", 2:T, 2:T), collapse = "\n")
    var_wy_rest <- paste(sprintf("wy%d ~~ wy%d", 2:T, 2:T), collapse = "\n")
  } else {
    var_wx_rest <- ""
    var_wy_rest <- ""
  }
  
  # Fix Observed Variances to Zero
  zero_var_x <- paste(sprintf("x%d ~~ 0*x%d", 1:T, 1:T), collapse = "\n")
  zero_var_y <- paste(sprintf("y%d ~~ 0*y%d", 1:T, 1:T), collapse = "\n")
  
  # Random Intercept Variances and Covariances
  
  ri_variances_and_covs <- c()
  
  # Global RI Variances and Covariance
  ri_variances_and_covs <- c(ri_variances_and_covs,
                             "RIgx ~~ varRIgx*RIgx",
                             "RIgy ~~ varRIgy*RIgy",
                             "RIgx ~~ covRIg*RIgy")
  
  # Segmented RI Variances and Within-segment Covariances
  
  # First: All RIx variances - CONSTRAINED EQUAL
  for (seg in 1:n_segments) {
    ri_variances_and_covs <- c(ri_variances_and_covs,
                               sprintf("RIx%d ~~ varRIx*RIx%d", seg, seg))
  }
  
  # Second: All RIy variances - CONSTRAINED EQUAL
  for (seg in 1:n_segments) {
    ri_variances_and_covs <- c(ri_variances_and_covs,
                               sprintf("RIy%d ~~ varRIy*RIy%d", seg, seg))
  }
  
  # Third: Within-segment RIx-RIy covariances - CONSTRAINED EQUAL ACROSS SEGMENTS
  for (seg in 1:n_segments) {
    ri_variances_and_covs <- c(ri_variances_and_covs,
                               sprintf("RIx%d ~~ covRI*RIy%d", seg, seg))
  }
  
  # Fix Global-to-Segmented RI Covariances to Zero
  global_seg_zero_covs <- c()
  for (seg in 1:n_segments) {
    global_seg_zero_covs <- c(global_seg_zero_covs,
                              sprintf("RIgx ~~ 0*RIx%d", seg),
                              sprintf("RIgx ~~ 0*RIy%d", seg),
                              sprintf("RIgy ~~ 0*RIx%d", seg),
                              sprintf("RIgy ~~ 0*RIy%d", seg))
  }
  
  # Across-segment Covariances for Segmented RIs
  ri_across_seg_covs <- c()
  across_seg_label <- ""
  
  if (n_segments > 1) {
    if (across_seg_cov == "zero") {
      # Fixed to zero: independent confounding across segments
      across_seg_label <- "(Fixed to Zero)"
      for (seg1 in 1:(n_segments - 1)) {
        for (seg2 in (seg1 + 1):n_segments) {
          ri_across_seg_covs <- c(ri_across_seg_covs,
                                  sprintf("RIx%d ~~ 0*RIx%d", seg1, seg2),
                                  sprintf("RIy%d ~~ 0*RIy%d", seg1, seg2),
                                  sprintf("RIx%d ~~ 0*RIy%d", seg1, seg2),
                                  sprintf("RIy%d ~~ 0*RIx%d", seg1, seg2))
        }
      }
      
    } else if (across_seg_cov == "free") {
      # Freely estimate all covariances with unique labels
      across_seg_label <- "(Freely Estimated)"
      for (seg1 in 1:(n_segments - 1)) {
        for (seg2 in (seg1 + 1):n_segments) {
          ri_across_seg_covs <- c(ri_across_seg_covs,
                                  sprintf("RIx%d ~~ covRIx_%d_%d*RIx%d", seg1, seg1, seg2, seg2),
                                  sprintf("RIy%d ~~ covRIy_%d_%d*RIy%d", seg1, seg1, seg2, seg2),
                                  sprintf("RIx%d ~~ covRIxy_%d_%d*RIy%d", seg1, seg1, seg2, seg2),
                                  sprintf("RIy%d ~~ covRIyx_%d_%d*RIx%d", seg1, seg1, seg2, seg2))
        }
      }
      
    } else if (across_seg_cov == "toeplitz") {
      # Block Toeplitz structure
      across_seg_label <- "(Block Toeplitz Structure)"
      
      max_lag <- n_segments - 1
      
      # Within-variable Toeplitz (X->X and Y->Y)
      for (lag in 1:max_lag) {
        for (seg1 in 1:(n_segments - lag)) {
          seg2 <- seg1 + lag
          
          ri_across_seg_covs <- c(ri_across_seg_covs,
                                  sprintf("RIx%d ~~ covRIx_lag%d*RIx%d", seg1, lag, seg2))
          ri_across_seg_covs <- c(ri_across_seg_covs,
                                  sprintf("RIy%d ~~ covRIy_lag%d*RIy%d", seg1, lag, seg2))
        }
      }
      
      # Cross-variable Toeplitz (X->Y and Y->X separate)
      for (lag in 1:max_lag) {
        for (seg1 in 1:(n_segments - lag)) {
          seg2 <- seg1 + lag
          
          ri_across_seg_covs <- c(ri_across_seg_covs,
                                  sprintf("RIx%d ~~ covRIxy_lag%d*RIy%d", seg1, lag, seg2))
          ri_across_seg_covs <- c(ri_across_seg_covs,
                                  sprintf("RIy%d ~~ covRIyx_lag%d*RIx%d", seg1, lag, seg2))
        }
      }
    }
  }
  
  # Fix RI Covariances with First State to Zero
  ri_first_state_zero <- c()
  
  # Global RIs with first state
  ri_first_state_zero <- c(ri_first_state_zero,
                           "wx1 ~~ 0*RIgx",
                           "wx1 ~~ 0*RIgy",
                           "wy1 ~~ 0*RIgx",
                           "wy1 ~~ 0*RIgy")
  
  # Segmented RIs with first state
  for (seg in 1:n_segments) {
    ri_first_state_zero <- c(ri_first_state_zero,
                             sprintf("wx1 ~~ 0*RIx%d", seg),
                             sprintf("wx1 ~~ 0*RIy%d", seg),
                             sprintf("wy1 ~~ 0*RIx%d", seg),
                             sprintf("wy1 ~~ 0*RIy%d", seg))
  }
  
  # Assemble Final Model String
  
  components <- c(
    "# 1. Random Intercepts",
    "# 1a. Global RIs (spanning all timepoints)",
    global_ri_x,
    global_ri_y,
    "# 1b. Segmented RIs - X then Y",
    paste(ri_x_definitions, collapse = "\n"),
    paste(ri_y_definitions, collapse = "\n"),
    "\n# 2. Within-person components",
    wx_defs,
    wy_defs,
    "\n# 3. Autoregressive and cross-lagged paths",
    paths_wx,
    paths_wy,
    "\n# 4. Covariances",
    "wx1 ~~ wy1 # Covariance at T1",
    res_covs,
    "\n# 5. (Co)variances of Random Intercepts",
    "# 5a. Global RI variances and covariance",
    paste(ri_variances_and_covs[1:3], collapse = "\n"),
    "# 5b. Segmented RI variances (constrained equal) and within-segment covariances (constrained equal)",
    paste(ri_variances_and_covs[4:length(ri_variances_and_covs)], collapse = "\n"),
    "# 5c. Global-to-Segmented RI covariances (fixed to zero)",
    paste(global_seg_zero_covs, collapse = "\n"),
    sprintf("\n# 5d. Across-segment covariances %s", 
            ifelse(n_segments > 1, across_seg_label, "")),
    if (length(ri_across_seg_covs) > 0) paste(ri_across_seg_covs, collapse = "\n") else "# (No across-segment covariances for single segment)",
    "\n# 6. Variances and residual variances of within-person components",
    "# 6a. First wave variances (exogenous)",
    var_wx1,
    var_wy1,
    "# 6b. Subsequent wave residual variances (endogenous)",
    var_wx_rest,
    var_wy_rest,
    "\n# 7. Fix observed variances to zero",
    zero_var_x,
    zero_var_y,
    "\n# 8. Fix RI covariances with first state to zero",
    paste(ri_first_state_zero, collapse = "\n")
  )
  
  components <- components[components != ""]
  model_string <- paste(components, collapse = "\n")
  
  return(model_string)
}

################################################################################
#                             ---- PLOTTING ----                               #
################################################################################

# Random Intercept variance plotting ---
plot_varRIx <- function(data, x_var, group_var, x_label, group_label, plot_title, 
                        xlim = NULL, ylim = NULL) {
  
  # ggplot call
  varRix_plot <- ggplot(data = data, 
                        aes(x = .data[[x_var]], 
                            y = varRIx, 
                            group = factor(.data[[group_var]]), 
                            color = factor(.data[[group_var]]))) +
    geom_line(linewidth = 1) +  
    geom_point(size = 2) +  
    
    # labels and titles
    labs(title = plot_title,
         x = x_label,
         y = "varRIx Value",
         color = group_label) + 
    
    # theme
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 16), 
      legend.position = "bottom" 
    )
  
  # Add axis limits if specified
  if (!is.null(xlim)) {
    varRix_plot <- varRix_plot + xlim(xlim[1], xlim[2])
  }
  if (!is.null(ylim)) {
    varRix_plot <- varRix_plot + ylim(ylim[1], ylim[2])
  }
  
  return(varRix_plot)
}

# cross-lagged parameter bias plotting ---
plot_crosslag_bias <- function(data, x_var, group_var, true_bx, true_by, 
                               x_label, group_label, title_suffix, 
                               xlim = NULL, ylim = NULL) {
  
  # Extract coefficients and calculate bias
  data$bx <- sapply(data$coefficients, function(coef_list) {
    if(!is.null(coef_list) && "bx" %in% names(coef_list)) {
      return(coef_list["bx"])
    } else {
      return(NA)
    }
  })
  
  data$by <- sapply(data$coefficients, function(coef_list) {
    if(!is.null(coef_list) && "by" %in% names(coef_list)) {
      return(coef_list["by"])
    } else {
      return(NA)
    }
  })
  
  data$bias_bx <- data$bx - true_bx
  data$bias_by <- data$by - true_by
  
  # Plot bias for bx (Y -> X cross-lagged path)
  bias_bx_plot <- ggplot(data = data, 
                         aes(x = .data[[x_var]], 
                             y = bias_bx, 
                             group = factor(.data[[group_var]]), 
                             color = factor(.data[[group_var]]))) +
    geom_line(linewidth = 1) +
    geom_point(size = 2) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "black", alpha = 0.7) +
    labs(title = paste0("Bias in Cross-lagged Path bx (Y→X) - ", title_suffix),
         x = x_label,
         y = "Bias in bx (Estimated - True)",
         color = group_label) +
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 16),
      legend.position = "bottom"
    )
  
  # Add axis limits if specified for bx plot
  if (!is.null(xlim)) {
    bias_bx_plot <- bias_bx_plot + xlim(xlim[1], xlim[2])
  }
  if (!is.null(ylim)) {
    bias_bx_plot <- bias_bx_plot + ylim(ylim[1], ylim[2])
  }
  
  # Plot bias for by (X -> Y cross-lagged path)  
  bias_by_plot <- ggplot(data = data, 
                         aes(x = .data[[x_var]], 
                             y = bias_by, 
                             group = factor(.data[[group_var]]), 
                             color = factor(.data[[group_var]]))) +
    geom_line(linewidth = 1) +
    geom_point(size = 2) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "black", alpha = 0.7) +
    labs(title = paste0("Bias in Cross-lagged Path by (X→Y) - ", title_suffix),
         x = x_label,
         y = "Bias in by (Estimated - True)",
         color = group_label) +
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 16),
      legend.position = "bottom"
    )
  
  # Add axis limits if specified for by plot
  if (!is.null(xlim)) {
    bias_by_plot <- bias_by_plot + xlim(xlim[1], xlim[2])
  }
  if (!is.null(ylim)) {
    bias_by_plot <- bias_by_plot + ylim(ylim[1], ylim[2])
  }
  
  return(list(bx_plot = bias_bx_plot, by_plot = bias_by_plot))
}

# Random intercept variance plotting as a percentage of total variance ---

### SHOULD BE FIXED OR DELETED - RUNNING THE DGM AGAIN IN FUNCTION IS DUMB AND REDUNDANT

plot_varRIx_proportion_auto <- function(data, x_var, group_var, x_label, group_label, plot_title, 
                                        xlim = NULL, ylim = NULL, ...) {
  
  # Source the DGM function
  source("functionSource.R")
  
  # Add a column for total variance in X based on conditions in the data
  data$total_var_x <- sapply(1:nrow(data), function(i) {
    current_row <- data[i, ]
    
    # Create a list of DGM parameters, starting with defaults
    dgm_params <- list(timepoints = 3)  # Always use 3 for population covariance
    
    # Add parameters based on what columns exist in the data
    if("autocorr" %in% colnames(data)) {
      dgm_params$autocorr_effects <- current_row$autocorr
    }
    
    if("N_participants" %in% colnames(data)) {
      dgm_params$N <- current_row$N_participants  # Map N_participants to N
    }
    
    if("N_confounders" %in% colnames(data)) {
      # Split confounders equally between X and Y variables
      dgm_params$n_x_confounders <- current_row$N_confounders
      dgm_params$n_y_confounders <- current_row$N_confounders
    }
    
    # Add any additional parameters passed through ...
    extra_params <- list(...)
    dgm_params <- c(dgm_params, extra_params)
    
    # Generate DGM parameters for this specific condition
    dgm_result <- do.call(DGM, dgm_params)
    
    # Extract total variance in X for this condition
    return(dgm_result$parameters$StatCOVmatrix[1,1])
  })
  
  # Calculate proportion of random intercept variance to total variance
  data$varRIx_proportion <- data$varRIx / data$total_var_x
  
  varRix_plot <- ggplot(data = data, 
                        aes(x = .data[[x_var]], 
                            y = varRIx_proportion, 
                            group = factor(.data[[group_var]]), 
                            color = factor(.data[[group_var]]))) +
    geom_line(linewidth = 1) +  # Add the line layer
    geom_point(size = 2) +     # Add points to mark the actual data points
    
    # --- Customize labels and titles ---
    labs(title = plot_title,
         x = x_label,
         y = "varRIx Proportion of Total Variance in X",
         color = group_label) + # This renames the legend title
    
    # --- Apply a clean theme ---
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 16), # Center the title
      legend.position = "bottom" # Move legend to the bottom
    )
  
  # Add axis limits if specified
  if (!is.null(xlim)) {
    varRix_plot <- varRix_plot + xlim(xlim[1], xlim[2])
  }
  if (!is.null(ylim)) {
    varRix_plot <- varRix_plot + ylim(ylim[1], ylim[2])
  }
  
  return(varRix_plot)
}