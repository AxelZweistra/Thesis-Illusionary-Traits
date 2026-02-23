### Functions for illusory trait thesis 

## data generating mechanism 

DGM <- function(N = 20000, 
                autocorr_effects = 0.2, 
                timepoints = 100, 
                n_x_confounders = 5, 
                n_y_confounders = 5,
                beta_focal = 0.1,
                beta_covfocal = 0.1,
                beta_covnofocal = 0.01,
                sigma = 0.5,
                seed = 427) 
  {

  # Set seed for reproducibility
  set.seed(seed)
  
  # Rename some parameters to match the variables used by the original paper
  T <- timepoints
  c1 <- n_x_confounders
  c2 <- n_y_confounders
  alpha <- autocorr_effects
  
  # Helper functions for matrix operations
  irow <- function(x) matrix(x, nrow = sqrt(length(x)), byrow = TRUE)
  row <- function(x) as.vector(t(x))
  
  # New objects for the simulation
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
  
  ## Check the max eigenvalue (to check stationarity)
  max_eigenvalue <- max(eigen(A)$values)
  if (max_eigenvalue >= 1) {
    warning(paste("Maximum eigenvalue is", round(max_eigenvalue, 4),"- the system may not be stationary"))
  }
  
  ## Sigma: Var-cov matrix of time-specific residuals
  Sigma <- diag(rep(sigma, v))
  
  ## Sigma1: Stationary Var-cov matrix for the first time point
  Sigma1 <- irow(solve(diag(dim(A)[1]^2) - t(A) %x% t(A)) %*% row(Sigma))
  
  ## Mean1: Stationary means for first time point
  Mean1 <- matrix(0, nrow = dim(A)[1], ncol = 1)
  
  # Generate data ---
  
  # create empty dataframe
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
      StatCOVmatrix = Sigma1,
      max_eigenvalue = max_eigenvalue
    ),
    coefficient_matrix = A
  ))
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
  
  # Variances of observed variables at all timepoints
  var_x <- paste(sprintf("x%d ~~ x%d", 1:T, 1:T), collapse = "\n")
  var_y <- paste(sprintf("y%d ~~ y%d", 1:T, 1:T), collapse = "\n")
  
  # Covariances between x and y at all timepoints
  cov_xy <- paste(sprintf("x%d ~~ x%d_y%d_cov*y%d", 1:T, 1:T, 1:T, 1:T), collapse = "\n")
  
  # Assemble final model string
  
  components <- c(
    "# 1. Autoregressive and cross-lagged paths",
    paths_x,
    paths_y,
    "\n# 2. Variances of observed variables",
    var_x,
    var_y,
    "\n# 3. Covariances between x and y at each timepoint",
    cov_xy
  )
  
  # Add together, skip empty strings and add linebreaks
  model_string <- paste(components[components != ""], collapse = "\n")
  
  return(model_string)
}

# ------------------------------------------------------------------------------

# Generate RI-CLPM model specification
# only parameter is T: number of timepoints

generate_riclpm <- function(T) {
  
  # Input validation
  if (!is.numeric(T) || T <= 2 || T %% 1 != 0) {
    stop("T must be an integer greater or equal to 3.") # RI-CLPM is not identified at T < 3
  }
  
  # Generate model components
  
  # Random intercepts (RIx and RIy)
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
  
  
  # Assemble final model string
  
  components <- c(
    "# 1. Random Intercepts", ri_x, ri_y,
    "\n# 2. Within-person components (measurement)", wx_defs, wy_defs,
    "\n# 3. Autoregressive and cross-lagged paths", paths_wx, paths_wy,
    "\n# 4. Covariances",
    "wx1 ~~ T1_cov*wy1 # Covariance at T1",
    res_covs,
    "\n# 5. (Co)variances of Random Intercepts",
    "RIx ~~ varRIx*RIx",
    "RIy ~~ varRIy*RIy",
    "RIx ~~ covRI*RIy",
    "\n# 6. (Residual) variances of within-person components (constrained equal)",
    var_wx,
    var_wy,
    "\n# 7. Fix observed variances to zero",
    zero_var_x,
    zero_var_y,
    "\n# 8. Fix RI covariances with first state to zero",
    "wx1 ~~ 0*RIx", "wx1 ~~ 0*RIy", "wy1 ~~ 0*RIx", "wy1 ~~ 0*RIy"
  )
  
  # add together, skip empty strings and add linebreaks
  model_string <- paste(components[components != ""], collapse = "\n")
  
  return(model_string)
}

# ------------------------------------------------------------------------------

# Generate MRI-CLPM model specification
# T: total number of timepoints, timespan: timepoints on which RI's are estimated
# across_seg_cov: has to be "zero", "free" or "toeplitz" 

generate_segmented_riclpm <- function(T, timespan = 3, across_seg_cov = "toeplitz") {
  
  # --- 1. Input Validation ---
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
  
  # --- 2. Determine Number of Segments ---
  n_segments <- ceiling(T / timespan)
  
  # Create a mapping of timepoints to segments
  timepoint_to_segment <- rep(1:n_segments, each = timespan)[1:T]
  
  # --- 3. Generate Random Intercepts for Each Segment ---
  # Structure: All RIx's first, then all RIy's
  
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
  
  # --- 4. Within-person Component Definitions ---
  wx_defs <- paste(sprintf("wx%d =~ 1*x%d", 1:T, 1:T), collapse = "\n")
  wy_defs <- paste(sprintf("wy%d =~ 1*y%d", 1:T, 1:T), collapse = "\n")
  
  # --- 5. Autoregressive and Cross-lagged Paths ---
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
  
  # --- 6. Constrained Correlated Residuals ---
  if (T > 1) {
    res_covs <- paste(sprintf("wx%d ~~ ur*wy%d", 2:T, 2:T), collapse = "\n")
  } else {
    res_covs <- ""
  }
  
  # --- 7. Variances of Within-person Components ---
  var_wx <- paste(sprintf("wx%d ~~ wx%d", 1:T, 1:T), collapse = "\n")
  var_wy <- paste(sprintf("wy%d ~~ wy%d", 1:T, 1:T), collapse = "\n")
  
  # --- 8. Fix Observed Variances to Zero ---
  zero_var_x <- paste(sprintf("x%d ~~ 0*x%d", 1:T, 1:T), collapse = "\n")
  zero_var_y <- paste(sprintf("y%d ~~ 0*y%d", 1:T, 1:T), collapse = "\n")
  
  # --- 9. Random Intercept Variances and Covariances ---
  
  # Within-segment variances and within-timespan covariances (RIx with RIy at same segment)
  ri_variances_and_within_seg_covs <- c()
  
  # First: All RIx variances
  for (seg in 1:n_segments) {
    ri_variances_and_within_seg_covs <- c(ri_variances_and_within_seg_covs,
                                          sprintf("RIx%d ~~ varRIx%d*RIx%d", seg, seg, seg))
  }
  
  # Second: All RIy variances
  for (seg in 1:n_segments) {
    ri_variances_and_within_seg_covs <- c(ri_variances_and_within_seg_covs,
                                          sprintf("RIy%d ~~ varRIy%d*RIy%d", seg, seg, seg))
  }
  
  # Third: Within-segment RIx-RIy covariances (same timespan only)
  for (seg in 1:n_segments) {
    ri_variances_and_within_seg_covs <- c(ri_variances_and_within_seg_covs,
                                          sprintf("RIx%d ~~ covRI%d*RIy%d", seg, seg, seg))
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
      # Toeplitz structure with cross-lagged covariances (DEFAULT)
      across_seg_label <- "(Toeplitz Structure)"
      
      max_lag <- n_segments - 1
      
      # Within-variable Toeplitz
      for (lag in 1:max_lag) {
        for (seg1 in 1:(n_segments - lag)) {
          seg2 <- seg1 + lag
          
          ri_across_seg_covs <- c(ri_across_seg_covs,
                                  sprintf("RIx%d ~~ covRIx_lag%d*RIx%d", seg1, lag, seg2))
          ri_across_seg_covs <- c(ri_across_seg_covs,
                                  sprintf("RIy%d ~~ covRIy_lag%d*RIy%d", seg1, lag, seg2))
        }
      }
      
      # Cross-variable Toeplitz
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
  
  # --- 10. Fix RI Covariances with First State to Zero ---
  ri_first_state_zero <- c()
  for (seg in 1:n_segments) {
    ri_first_state_zero <- c(ri_first_state_zero,
                             sprintf("wx1 ~~ 0*RIx%d", seg),
                             sprintf("wx1 ~~ 0*RIy%d", seg),
                             sprintf("wy1 ~~ 0*RIx%d", seg),
                             sprintf("wy1 ~~ 0*RIy%d", seg))
  }
  
  # --- 11. Assemble Final Model String ---
  
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
    "# 5a. Variances and within-segment covariances",
    paste(ri_variances_and_within_seg_covs, collapse = "\n"),
    sprintf("\n# 5b. Across-segment covariances %s", 
            ifelse(n_segments > 1, across_seg_label, "")),
    if (length(ri_across_seg_covs) > 0) paste(ri_across_seg_covs, collapse = "\n") else "# (No across-segment covariances for single segment)",
    "\n# 6. (Residual) variances of within-person components",
    var_wx,
    var_wy,
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