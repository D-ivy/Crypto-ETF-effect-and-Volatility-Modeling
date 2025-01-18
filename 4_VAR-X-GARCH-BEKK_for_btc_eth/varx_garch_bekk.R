# Part 1A: Data Preparation and Model Definition

# Load required packages
required_packages <- c("quantmod", "vars", "zoo", "xts", "numDeriv", "sandwich", "gridExtra", "ggplot2", "grid", "gridExtra")
for(pkg in required_packages) {
  if (!require(pkg, character.only = TRUE)) {
    install.packages(pkg)
    library(pkg, character.only = TRUE)
  }
}

# Set seed for reproducibility
set.seed(123)

# Data preparation function
# Data preparation function
prepare_data <- function(start_date, end_date) {
  tryCatch({
    # Get cryptocurrency data and VIX
    getSymbols("BTC-USD", from = start_date, to = end_date, src = "yahoo")
    getSymbols("ETH-USD", from = start_date, to = end_date, src = "yahoo")
    getSymbols("^VIX", from = start_date, to = end_date, src = "yahoo")
    
    # Calculate returns for cryptocurrencies (in percentage)
    btc_returns <- diff(log(`BTC-USD`$`BTC-USD.Adjusted`)) * 1
    eth_returns <- diff(log(`ETH-USD`$`ETH-USD.Adjusted`)) * 1
    
    # Process VIX: First convert to same frequency and fill missing values
    vix_data <- VIX$VIX.Adjusted
    # Create a sequence of dates matching crypto returns
    date_seq <- index(btc_returns)
    
    # Convert VIX to zoo and align with crypto dates
    vix_zoo <- zoo(as.numeric(vix_data), order.by = index(vix_data))
    vix_aligned <- merge(vix_zoo, zoo(, date_seq))
    
    # Fill missing values using previous day's value
    vix_filled <- na.locf(vix_aligned, fromLast = FALSE)
    
    # Ensure all series are zoo objects with matching dates
    btc_zoo <- zoo(btc_returns, order.by = index(btc_returns))
    eth_zoo <- zoo(eth_returns, order.by = index(eth_returns))
    
    # Merge all series to ensure date alignment
    all_series <- merge(btc_zoo, eth_zoo, vix_filled)
    colnames(all_series) <- c("BTC", "ETH", "VIX")
    
    # Remove any remaining NA values
    all_series <- na.omit(all_series)
    
    # Extract final aligned series
    btc_final <- all_series[, "BTC"]
    eth_final <- all_series[, "ETH"]
    vix_final <- all_series[, "VIX"]
    
    # Print diagnostic information
    cat("\nData Summary:\n")
    cat("============\n")
    cat(sprintf("BTC observations: %d\n", length(btc_final)))
    cat(sprintf("ETH observations: %d\n", length(eth_final)))
    cat(sprintf("VIX observations: %d\n", length(vix_final)))
    cat(sprintf("Date range: %s to %s\n", 
                format(start(all_series)), format(end(all_series))))
    
    # Safe summary statistics function
    safe_summary <- function(x, name) {
      if (length(x) == 0) return(NULL)
      stats <- c(
        Mean = mean(x, na.rm = TRUE),
        SD = sd(x, na.rm = TRUE),
        Min = min(x, na.rm = TRUE),
        Max = max(x, na.rm = TRUE),
        Q1 = as.numeric(quantile(x, 0.25, na.rm = TRUE)),
        Median = median(x, na.rm = TRUE),
        Q3 = as.numeric(quantile(x, 0.75, na.rm = TRUE))
      )
      cat(sprintf("\n%s Summary:\n", name))
      print(round(stats, 4))
      return(stats)
    }
    
    # Print summaries
    safe_summary(btc_final, "BTC")
    safe_summary(eth_final, "ETH")
    safe_summary(vix_final, "VIX")
    
    # Print sample data (first 5 rows)
    cat("\nSample Data (First 5 rows):\n")
    cat("=========================\n")
    cat("\nBTC Returns:\n")
    print(head(btc_final, 5))
    cat("\nETH Returns:\n")
    print(head(eth_final, 5))
    cat("\nVIX Values:\n")
    print(head(vix_final, 5))
    
    return(list(
      btc = btc_final,
      eth = eth_final,
      vix = vix_final
    ))
  }, error = function(e) {
    cat("Error in data preparation:", conditionMessage(e), "\n")
    stop(e)
  })
}

# Robust VAR estimation function
robust_var <- function(Y, X, p = 1) {
  tryCatch({
    T <- nrow(Y)
    K <- ncol(Y)
    
    # Input validation
    if (is.null(Y) || is.null(X)) stop("Input data cannot be NULL")
    if (nrow(X) != nrow(Y)) stop("Y and X must have same number of rows")
    
    # Create lagged Y matrices
    Y_lags <- matrix(0, T-p, K*p)
    for(i in 1:p) {
      Y_lags[,(1+K*(i-1)):(K*i)] <- Y[(p-i+1):(T-i),]
    }
    
    # Prepare matrices for estimation
    Y_est <- Y[(p+1):T,]
    X_est <- cbind(1, Y_lags, X[(p+1):T,])
    
    # Check for rank deficiency
    if (qr(X_est)$rank < ncol(X_est)) {
      stop("X matrix is rank deficient")
    }
    
    # QMLE estimation
    beta_qmle <- solve(crossprod(X_est)) %*% crossprod(X_est, Y_est)
    residuals <- Y_est - X_est %*% beta_qmle
    
    # Robust covariance estimation (HC3)
    n <- nrow(X_est)
    k <- ncol(X_est)
    h <- diag(X_est %*% solve(crossprod(X_est)) %*% t(X_est))
    weights <- 1 / (1 - pmin(h, 0.99))^2  # Bounded leverage adjustment
    
    # Calculate robust standard errors
    meat <- matrix(0, k, k)
    for(i in 1:n) {
      xi <- matrix(X_est[i,], ncol=1)
      ri <- matrix(residuals[i,], ncol=1)
      meat <- meat + weights[i] * (xi %*% t(xi)) * (t(ri) %*% ri)
    }
    
    vcov_robust <- solve(crossprod(X_est)) %*% meat %*% solve(crossprod(X_est))
    se_robust <- sqrt(pmax(diag(vcov_robust), 0))
    
    # Calculate statistics
    coef_matrix <- matrix(beta_qmle, nrow=k)
    se_matrix <- matrix(se_robust, nrow=k)
    t_stats <- coef_matrix / se_matrix
    p_values <- 2 * (1 - pnorm(abs(t_stats)))
    
    # R-squared statistics
    r2 <- sapply(1:ncol(Y_est), function(i) {
      1 - sum(residuals[,i]^2) / sum((Y_est[,i] - mean(Y_est[,i]))^2)
    })
    adj_r2 <- 1 - (1 - r2) * (n - 1)/(n - k)
    
    return(list(
      coefficients = coef_matrix,
      std.errors = se_matrix,
      t_statistics = t_stats,
      p_values = p_values,
      residuals = residuals,
      vcov = vcov_robust,
      r.squared = r2,
      adj.r.squared = adj_r2,
      nobs = n,
      ncoef = k
    ))
  }, error = function(e) {
    cat("Error in VAR estimation:", conditionMessage(e), "\n")
    stop(e)
  })
}

# Print function for VAR results
print_var_results <- function(var_results, equation_names, variable_names) {
  tryCatch({
    cat("\nRobust VAR Model Results:\n")
    cat("======================\n")
    
    for(i in 1:length(equation_names)) {
      cat(sprintf("\n%s Equation:\n", equation_names[i]))
      cat(paste(rep("-", nchar(equation_names[i]) + 10), collapse=""), "\n")
      
      for(j in 1:length(variable_names)) {
        coef <- var_results$coefficients[j,i]
        se <- var_results$std.errors[j,i]
        t_stat <- var_results$t_statistics[j,i]
        p_val <- var_results$p_values[j,i]
        
        stars <- ifelse(p_val < 0.01, "***",
                        ifelse(p_val < 0.05, "**",
                               ifelse(p_val < 0.1, "*", "")))
        
        cat(sprintf("%-15s %8.4f (%6.4f) %s\n", 
                    paste0(variable_names[j], ":"), coef, se, stars))
      }
      
      cat(sprintf("\nR-squared:          %8.4f\n", var_results$r.squared[i]))
      cat(sprintf("Adjusted R-squared:  %8.4f\n", var_results$adj.r.squared[i]))
    }
    
    cat("\nSignificance codes: *** p<0.01, ** p<0.05, * p<0.1")
    cat(sprintf("\nObservations: %d", var_results$nobs))
    cat(sprintf("\nCoefficients: %d\n", var_results$ncoef))
    
  }, error = function(e) {
    cat("Error in printing results:", conditionMessage(e), "\n")
  })
}

# Part 1B: BEKK-GARCH Estimation


# Enhanced BEKK function with QMLE estimation
# Enhanced BEKK function with QMLE estimation
custom_bekk <- function(returns, dummy, vix_lag) {
  T <- nrow(returns)
  K <- ncol(returns)
  n_params <- 15
  
  # Scale returns for numerical stability
  scale_factor <- mean(apply(returns, 2, sd))
  scaled_returns <- returns / scale_factor
  
  # Create storage environment with increased stability threshold
  storage <- new.env()
  storage$H <- array(0, c(2, 2, T))
  storage$last_params <- NULL
  
  # Enhanced function to ensure positive definiteness
  ensure_pd <- function(H) {
    tryCatch({
      if(any(is.na(H)) || any(is.infinite(H))) {
        H <- cov(scaled_returns)
      }
      H <- (H + t(H))/2
      # Increased stabilization factor
      H <- H + diag(1e-4, nrow(H))
      
      eig <- eigen(H, symmetric = TRUE)
      if(any(eig$values < 1e-6)) {  # Increased eigenvalue threshold
        values <- pmax(eig$values, 1e-6)
        H <- eig$vectors %*% diag(values) %*% t(eig$vectors)
      }
      return(H)
    }, error = function(e) {
      return(cov(scaled_returns) + diag(1e-4, K))
    })
  }
  
  # Initialize with more stable starting values
  cov_matrix <- cov(scaled_returns)
  init_C <- try(chol(cov_matrix + diag(1e-4, K)), silent = TRUE)
  if(inherits(init_C, "try-error")) {
    init_C <- diag(sqrt(diag(cov_matrix)))
  }
  
  init_params <- c(
    init_C[1,1], init_C[2,1], init_C[2,2],  # C matrix
    0.15, 0.05, 0.05, 0.15,   # A matrix (reduced values)
    0.85, 0.05, 0.05, 0.85,   # G matrix (increased persistence)
    0.03, 0.03,               # A_star matrix (reduced spillover)
    0.03, 0.03                # G_star matrix (reduced spillover)
  )
  
  # Adjusted parameter bounds for better stability
  lower_bounds <- c(
    0.01, -0.5, 0.01,      # C matrix
    0, -0.15, -0.15, 0,    # A matrix
    0.6, -0.15, -0.15, 0.6, # G matrix
    -0.5, -0.5,            # A_star matrix
    -0.5, -0.5             # G_star matrix
  )
  
  upper_bounds <- c(
    2, 0.5, 2,           # C matrix
    0.5, 0.3, 0.3, 0.5,  # A matrix
    0.98, 0.3, 0.3, 0.98, # G matrix
    0.5, 0.5,            # A_star matrix
    0.5, 0.5             # G_star matrix
  )
  
  # Modified QMLE likelihood function with improved numerical stability
  negloglik <- function(params) {
    tryCatch({
      storage$last_params <- params
      
      # Reconstruct matrices from parameters
      C <- matrix(0, 2, 2)
      C[lower.tri(C, diag=TRUE)] <- params[1:3]
      
      A <- matrix(params[4:7], 2, 2)
      G <- matrix(params[8:11], 2, 2)
      
      A_star <- matrix(0, 2, 2)
      A_star[1,2] <- params[12]
      A_star[2,1] <- params[13]
      
      G_star <- matrix(0, 2, 2)
      G_star[1,2] <- params[14]
      G_star[2,1] <- params[15]
      
      # Initialize H array and likelihood with stabilized covariance
      storage$H[,,1] <- ensure_pd(cov(scaled_returns))
      ll <- 0
      
      # Main BEKK recursion with additional checks
      for(t in 2:T) {
        eps_t_1 <- matrix(scaled_returns[t-1,], ncol=1)
        
        # Apply structural break effects
        A_t <- A
        G_t <- G
        if(dummy[t] == 1) {
          A_t <- A + A_star
          G_t <- G + G_star
        }
        
        # BEKK equation with stability checks
        H_next <- crossprod(C) + 
          A_t %*% (eps_t_1 %*% t(eps_t_1)) %*% t(A_t) + 
          G_t %*% storage$H[,,t-1] %*% t(G_t)
        
        storage$H[,,t] <- ensure_pd(H_next)
        
        # Likelihood contribution with additional checks
        eps_t <- matrix(scaled_returns[t,], ncol=1)
        H_inv <- try(solve(storage$H[,,t]), silent=TRUE)
        if(inherits(H_inv, "try-error")) return(1e10)
        
        det_H <- det(storage$H[,,t])
        if(det_H <= 1e-10) return(1e10)
        
        quad_form <- t(eps_t) %*% H_inv %*% eps_t
        if(quad_form < 0 || is.na(quad_form)) return(1e10)
        
        ll <- ll + log(det_H) + quad_form
      }
      
      return(ll/2)
    }, error = function(e) {
      return(1e10)
    })
  }
  
  # Optimize with modified control parameters
  cat("\nStarting BEKK-GARCH optimization...\n")
  opt <- optim(init_params,
               negloglik,
               method = "L-BFGS-B",
               lower = lower_bounds,
               upper = upper_bounds,
               control = list(
                 maxit = 10000,      # Increased maximum iterations
                 factr = 1e8,        # Relaxed convergence criterion
                 pgtol = 1e-6,       # Relaxed gradient tolerance
                 trace = 1
               ))
  
  # Rest of the function remains the same
  final_params <- storage$last_params
  C <- matrix(0, 2, 2)
  C[lower.tri(C, diag=TRUE)] <- final_params[1:3]
  A <- matrix(final_params[4:7], 2, 2)
  G <- matrix(final_params[8:11], 2, 2)
  A_star <- matrix(0, 2, 2)
  A_star[1,2] <- final_params[12]
  A_star[2,1] <- final_params[13]
  G_star <- matrix(0, 2, 2)
  G_star[1,2] <- final_params[14]
  G_star[2,1] <- final_params[15]
  
  # Calculate robust standard errors
  scores <- calc_scores(final_params, scaled_returns, dummy)
  B <- crossprod(scores)/(T-1)
  A_mat <- try(solve(optimHess(final_params, negloglik)), silent=TRUE)
  
  # Initialize SE matrices
  se_C <- matrix(0, 2, 2)
  se_A <- matrix(0, 2, 2)
  se_G <- matrix(0, 2, 2)
  se_A_star <- matrix(0, 2, 2)
  se_G_star <- matrix(0, 2, 2)
  
  if(!inherits(A_mat, "try-error")) {
    robust_vcov <- A_mat %*% B %*% A_mat
    se <- sqrt(pmax(diag(robust_vcov), 0))
    
    se_C[lower.tri(se_C, diag=TRUE)] <- se[1:3]
    se_A[1:2, 1:2] <- matrix(se[4:7], 2, 2)
    se_G[1:2, 1:2] <- matrix(se[8:11], 2, 2)
    se_A_star[1,2] <- se[12]
    se_A_star[2,1] <- se[13]
    se_G_star[1,2] <- se[14]
    se_G_star[2,1] <- se[15]
  }
  
  return(list(
    H = storage$H * scale_factor^2,
    C = C * scale_factor,
    se_C = se_C * scale_factor,
    A = A,
    se_A = se_A,
    G = G,
    se_G = se_G,
    A_star = A_star,
    se_A_star = se_A_star,
    G_star = G_star,
    se_G_star = se_G_star,
    convergence = opt$convergence,
    likelihood = -opt$value,
    scale_factor = scale_factor
  ))
}

# Helper function to calculate scores for robust standard errors
calc_scores <- function(params, returns, dummy) {
  T <- nrow(returns)
  n_params <- length(params)
  scores <- matrix(0, T-1, n_params)
  
  # Create an environment to store H
  storage <- new.env()
  storage$H <- array(0, c(2, 2, T))
  
  # Reconstruct matrices
  C <- matrix(0, 2, 2)
  C[lower.tri(C, diag=TRUE)] <- params[1:3]
  
  A <- matrix(params[4:7], 2, 2)
  G <- matrix(params[8:11], 2, 2)
  
  A_star <- matrix(0, 2, 2)
  A_star[1,2] <- params[12]
  A_star[2,1] <- params[13]
  
  G_star <- matrix(0, 2, 2)
  G_star[1,2] <- params[14]
  G_star[2,1] <- params[15]
  
  # Function to ensure positive definiteness
  ensure_pd <- function(H) {
    tryCatch({
      if(any(is.na(H)) || any(is.infinite(H))) {
        H <- cov(returns)
      }
      H <- (H + t(H))/2
      H <- H + diag(1e-6, nrow(H))
      
      eig <- eigen(H, symmetric = TRUE)
      if(any(eig$values < 1e-8)) {
        values <- pmax(eig$values, 1e-8)
        H <- eig$vectors %*% diag(values) %*% t(eig$vectors)
      }
      return(H)
    }, error = function(e) {
      return(cov(returns) + diag(1e-6, ncol(returns)))
    })
  }
  
  # Calculate H matrices
  storage$H[,,1] <- ensure_pd(cov(returns))
  
  for(t in 2:T) {
    eps_t_1 <- matrix(returns[t-1,], ncol=1)
    A_t <- A
    G_t <- G
    if(dummy[t] == 1) {
      A_t <- A + A_star
      G_t <- G + G_star
    }
    
    H_next <- crossprod(C) + 
      A_t %*% (eps_t_1 %*% t(eps_t_1)) %*% t(A_t) + 
      G_t %*% storage$H[,,t-1] %*% t(G_t)
    
    storage$H[,,t] <- ensure_pd(H_next)
    
    # Calculate score for each parameter
    eps_t <- matrix(returns[t,], ncol=1)
    H_inv <- solve(storage$H[,,t])
    qt <- H_inv %*% eps_t %*% t(eps_t) %*% H_inv - H_inv
    
    score_vec <- numeric(n_params)
    score_vec[1:3] <- c(qt[1,1] * C[1,1],
                        qt[2,1] * C[1,1] + qt[2,2] * C[2,1],
                        qt[2,2] * C[2,2])
    score_vec[4:7] <- c(qt[1,1] * eps_t_1[1]^2,
                        qt[1,2] * eps_t_1[1]*eps_t_1[2],
                        qt[2,1] * eps_t_1[1]*eps_t_1[2],
                        qt[2,2] * eps_t_1[2]^2)
    score_vec[8:11] <- c(qt[1,1] * storage$H[1,1,t-1],
                         qt[1,2] * storage$H[1,2,t-1],
                         qt[2,1] * storage$H[2,1,t-1],
                         qt[2,2] * storage$H[2,2,t-1])
    score_vec[12:13] <- c(qt[1,2] * eps_t_1[1]*eps_t_1[2],
                          qt[2,1] * eps_t_1[1]*eps_t_1[2]) * dummy[t]
    score_vec[14:15] <- c(qt[1,2] * storage$H[1,2,t-1],
                          qt[2,1] * storage$H[2,1,t-1]) * dummy[t]
    
    scores[t-1,] <- 0.5 * score_vec
  }
  
  return(scores)
}

# Function for Wald test of structural break parameters
wald_test_break <- function(bekk_results) {
  # Extract star parameters
  star_params <- c(
    bekk_results$A_star[1,2], bekk_results$A_star[2,1],
    bekk_results$G_star[1,2], bekk_results$G_star[2,1]
  )
  
  # Get standard errors for star parameters
  star_se <- c(
    bekk_results$se_A_star[1,2], bekk_results$se_A_star[2,1],
    bekk_results$se_G_star[1,2], bekk_results$se_G_star[2,1]
  )
  
  # Create variance-covariance matrix (assuming independence for simplicity)
  vcov_star <- diag(star_se^2)
  
  # Calculate Wald statistic
  W <- t(star_params) %*% solve(vcov_star) %*% star_params
  
  # Degrees of freedom = number of restrictions
  df <- length(star_params)
  
  # Calculate p-value
  p_value <- 1 - pchisq(W, df = df)
  
  return(list(
    W_stat = as.numeric(W),
    df = df,
    p_value = p_value,
    star_params = star_params,
    star_se = star_se
  ))
}

# Part 2A: Model Fitting and Initial Analysis

# Grid search function for BEKK parameters
grid_search_bekk <- function(returns, grid_points = 5) {
  # Initialize variables
  n <- ncol(returns)
  best_likelihood <- -Inf
  best_params <- NULL
  
  # Create grid values for A and G matrices
  a_values <- seq(0.1, 0.5, length.out = grid_points)
  g_values <- seq(0.7, 0.99, length.out = grid_points)
  
  # Initialize progress tracking
  total_iterations <- grid_points^4
  current_iteration <- 0
  
  # Calculate initial C matrix
  cov_matrix <- cov(returns)
  init_C <- chol(cov_matrix)
  
  # Grid search over main diagonal elements
  for(a11 in a_values) {
    for(a22 in a_values) {
      for(g11 in g_values) {
        for(g22 in g_values) {
          current_iteration <- current_iteration + 1
          if(current_iteration %% 10 == 0) {
            cat(sprintf("\rProgress: %.1f%%", 100 * current_iteration/total_iterations))
          }
          
          # Create initial parameter vector
          init_params <- c(
            init_C[1,1], init_C[2,1], init_C[2,2],  # C matrix
            a11, 0.1, 0.1, a22,    # A matrix
            g11, 0.05, 0.05, g22,  # G matrix
            0.05, 0.05,            # A_star matrix
            0.05, 0.05             # G_star matrix
          )
          
          # Try these parameters
          try({
            result <- optim(init_params,
                            negloglik,
                            method = "L-BFGS-B",
                            lower = lower_bounds,
                            upper = upper_bounds,
                            control = list(maxit = 1000))
            
            if(-result$value > best_likelihood) {
              best_likelihood <- -result$value
              best_params <- result$par
              cat(sprintf("\nNew best likelihood: %.2f\n", best_likelihood))
            }
          }, silent = TRUE)
        }
      }
    }
  }
  return(list(parameters = best_params, likelihood = best_likelihood))
}

# Format matrix with standard errors function
format_matrix_with_se <- function(mat, se_mat, title, structural = FALSE) {
  cat("\n", title, ":\n")
  for(i in 1:nrow(mat)) {
    for(j in 1:ncol(mat)) {
      if (structural && i == j) {
        cat("---\t")
        next
      }
      
      cat(sprintf("%.3f", mat[i,j]))
      if(!is.na(se_mat[i,j]) && se_mat[i,j] != 0) {
        t_stat <- mat[i,j] / se_mat[i,j]
        p_value <- 2 * (1 - pnorm(abs(t_stat)))
        stars <- ifelse(p_value < 0.01, "***",
                        ifelse(p_value < 0.05, "**",
                               ifelse(p_value < 0.1, "*", "")))
        # Changed here - now showing p-value instead of se
        cat(sprintf(" (%.3f)%s", p_value, stars))
      }
      cat("\t")
    }
    cat("\n")
  }
}

# Set date range
start_date <- as.Date("2022-5-01")
end_date <- as.Date("2024-12-31")
break_date <- as.Date("2024-1-10")

# Prepare data
cat("Preparing data...\n")
# After preparing initial data
data_list <- prepare_data(start_date, end_date)

# Create merged dataset with lagged VIX
all_data <- merge(data_list$btc, data_list$eth)
colnames(all_data) <- c("BTC", "ETH")

# Create lagged VIX series
vix_lagged <- lag(data_list$vix, 1)
all_data <- merge(all_data, vix_lagged)
colnames(all_data)[3] <- "VIX_LAG"

# Create dummy variable for structural break
dummy <- zoo(ifelse(time(all_data) >= break_date, 1, 0), 
             order.by = index(all_data))
all_data <- merge(all_data, dummy)
colnames(all_data)[4] <- "D1"
all_data <- na.omit(all_data)

# Fit VAR model with lagged VIX
cat("Fitting VAR model...\n")
Y <- all_data[, c("BTC", "ETH")]
X <- cbind(all_data$D1, all_data$VIX_LAG)  # Now using lagged VIX
colnames(X) <- c("D1", "VIX_LAG")

var_model <- VAR(Y, p = 1, type = "const", exogen = X)

# Get residuals
residuals <- residuals(var_model)

# Print VAR summary
cat("\nVAR Model Summary:\n")
print(summary(var_model))

# Perform grid search for initial parameters
cat("\nPerforming grid search for optimal starting values...\n")
grid_results <- grid_search_bekk(as.matrix(residuals), grid_points = 5)
cat("\nGrid search complete. Best likelihood:", grid_results$likelihood, "\n")

# Use grid search results for initialization
init_params <- grid_results$parameters

# Fit BEKK model with grid search optimized starting values
cat("\nFitting BEKK model with QMLE estimation using optimized starting values...\n")
# In your main code, update the BEKK call to:
bekk_results <- custom_bekk(as.matrix(residuals), all_data$D1, all_data$VIX_LAG)

# Check convergence and print results
if(bekk_results$convergence == 0) {
  cat("\nBEKK Model converged successfully.\n")
  cat("Scale factor used:", sprintf("%.3f", bekk_results$scale_factor), "\n")
  
  # Print parameters with robust standard errors
  format_matrix_with_se(bekk_results$C, bekk_results$se_C, "C matrix (Constants)")
  format_matrix_with_se(bekk_results$A, bekk_results$se_A, "A matrix (ARCH effects)")
  format_matrix_with_se(bekk_results$G, bekk_results$se_G, "G matrix (GARCH effects)")
  format_matrix_with_se(bekk_results$A_star, bekk_results$se_A_star, 
                        "A_star matrix (Spillover ARCH effects)", structural = TRUE)
  format_matrix_with_se(bekk_results$G_star, bekk_results$se_G_star, 
                        "G_star matrix (Spillover GARCH effects)", structural = TRUE)
  
  cat("\nNote: '---' indicates constrained parameters in structural matrices")
  cat("\nSignificance codes: *** p<0.01, ** p<0.05, * p<0.1")
  cat("\nStandard errors are robust (Bollerslev-Wooldridge)\n")
  
  # Add the Wald test here
  wald_results <- wald_test_break(bekk_results)
  
  # Print Wald test results
  cat("\nWald Test for Structural Break:\n")
  cat("===============================\n")
  cat(sprintf("Wald statistic: %.4f\n", wald_results$W_stat))
  cat(sprintf("Degrees of freedom: %d\n", wald_results$df))
  cat(sprintf("p-value: %.4f\n", wald_results$p_value))
  
  # Print individual parameter contributions
  cat("\nIndividual Parameter Tests:\n")
  for(i in 1:length(wald_results$star_params)) {
    param_name <- c("A12_star", "A21_star", "G12_star", "G21_star")[i]
    t_stat <- wald_results$star_params[i] / wald_results$star_se[i]
    p_val <- 2 * (1 - pnorm(abs(t_stat)))
    cat(sprintf("%s: %.4f (SE: %.4f, p-value: %.4f)\n", 
                param_name, 
                wald_results$star_params[i],
                wald_results$star_se[i],
                p_val))
  }
} else {
  cat("\nWarning: BEKK Model did not converge. Results may not be reliable.\n")
}

# Part 2B: Model Validation and Analysis

# Model Validation Section
cat("\nModel Validation Metrics:\n")
cat("=======================\n")

# Calculate standardized residuals
std_resid <- matrix(0, nrow(residuals), 2)
for(t in 1:nrow(residuals)) {
  H_t <- matrix(c(bekk_results$H[1,1,t], bekk_results$H[1,2,t],
                  bekk_results$H[1,2,t], bekk_results$H[2,2,t]), 2, 2)
  H_t_sqrt <- try(chol(H_t), silent=TRUE)
  if(inherits(H_t_sqrt, "try-error")) {
    H_t_sqrt <- diag(sqrt(diag(H_t)))
  }
  std_resid[t,] <- solve(t(H_t_sqrt)) %*% matrix(residuals[t,], ncol=1)
}

# Ljung-Box Test for Return Series
lb_test_btc_ret <- Box.test(residuals[,1], lag=10, type="Ljung-Box")
lb_test_eth_ret <- Box.test(residuals[,2], lag=10, type="Ljung-Box")

# Ljung-Box Test for Squared Standardized Residuals
lb_test_btc_sq <- Box.test(std_resid[,1]^2, lag=10, type="Ljung-Box")
lb_test_eth_sq <- Box.test(std_resid[,2]^2, lag=10, type="Ljung-Box")

# Calculate Log-Likelihood
log_lik <- -bekk_results$likelihood

# Calculate Information Criteria
n_params <- 15  # Total number of parameters in the BEKK model
n_obs <- nrow(residuals)
aic <- 2 * log_lik + 2 * n_params
bic <- 2 * log_lik + n_params * log(n_obs)

# Print Results
cat("\nLog-Likelihood:", sprintf("%.2f", -log_lik))
cat("\nAIC:", sprintf("%.2f", aic))
cat("\nBIC:", sprintf("%.2f", bic))

cat("\n\nLjung-Box Tests (lag=10):")
cat("\nBTC Returns:")
cat("\n  Q-statistic:", sprintf("%.4f", lb_test_btc_ret$statistic))
cat("\n  p-value:", sprintf("%.4f", lb_test_btc_ret$p.value))

cat("\nETH Returns:")
cat("\n  Q-statistic:", sprintf("%.4f", lb_test_eth_ret$statistic))
cat("\n  p-value:", sprintf("%.4f", lb_test_eth_ret$p.value))

cat("\n\nLjung-Box Tests for Squared Standardized Residuals (lag=10):")
cat("\nBTC:")
cat("\n  Q-statistic:", sprintf("%.4f", lb_test_btc_sq$statistic))
cat("\n  p-value:", sprintf("%.4f", lb_test_btc_sq$p.value))

cat("\nETH:")
cat("\n  Q-statistic:", sprintf("%.4f", lb_test_eth_sq$statistic))
cat("\n  p-value:", sprintf("%.4f", lb_test_eth_sq$p.value))

# Additional Descriptive Statistics for Standardized Residuals
std_resid_stats <- function(x) {
  c(Mean = mean(x),
    SD = sd(x),
    Skewness = sum((x - mean(x))^3)/(length(x) * sd(x)^3),
    Kurtosis = sum((x - mean(x))^4)/(length(x) * sd(x)^4),
    JB_stat = length(x) * (sum((x - mean(x))^3)^2/(6 * length(x) * sd(x)^6) +
                             (sum((x - mean(x))^4)/(length(x) * sd(x)^4) - 3)^2/24))
}

btc_stats <- std_resid_stats(std_resid[,1])
eth_stats <- std_resid_stats(std_resid[,2])

cat("\n\nStandardized Residuals Statistics:")
cat("\nBTC:")
cat("\n  Mean:", sprintf("%.4f", btc_stats["Mean"]))
cat("\n  SD:", sprintf("%.4f", btc_stats["SD"]))
cat("\n  Skewness:", sprintf("%.4f", btc_stats["Skewness"]))
cat("\n  Kurtosis:", sprintf("%.4f", btc_stats["Kurtosis"]))
cat("\n  Jarque-Bera stat:", sprintf("%.4f", btc_stats["JB_stat"]))
cat("\n  JB p-value:", sprintf("%.4f", pchisq(btc_stats["JB_stat"], df=2, lower.tail=FALSE)))

cat("\nETH:")
cat("\n  Mean:", sprintf("%.4f", eth_stats["Mean"]))
cat("\n  SD:", sprintf("%.4f", eth_stats["SD"]))
cat("\n  Skewness:", sprintf("%.4f", eth_stats["Skewness"]))
cat("\n  Kurtosis:", sprintf("%.4f", eth_stats["Kurtosis"]))
cat("\n  Jarque-Bera stat:", sprintf("%.4f", eth_stats["JB_stat"]))
cat("\n  JB p-value:", sprintf("%.4f", pchisq(eth_stats["JB_stat"], df=2, lower.tail=FALSE)))

# Extract and validate conditional measures
cat("\nExtracting conditional measures...\n")
h11 <- sapply(1:nrow(residuals), function(t) bekk_results$H[1,1,t])
h22 <- sapply(1:nrow(residuals), function(t) bekk_results$H[2,2,t])
h12 <- sapply(1:nrow(residuals), function(t) bekk_results$H[1,2,t])

# Ensure positive variances and valid correlations
h11 <- pmax(h11, 1e-10)
h22 <- pmax(h22, 1e-10)

# Create results dataframe with proper dates and measures
results <- data.frame(
  Date = as.POSIXct(index(residuals)),
  BTC_Returns = all_data$BTC[-1],  # Remove first observation to match residuals
  ETH_Returns = all_data$ETH[-1],  # Remove first observation to match residuals
  BTC_Vol = sqrt(h11),
  ETH_Vol = sqrt(h22),
  Correlation = pmin(pmax(h12/sqrt(h11 * h22), -1), 1),
  Spillover = abs(h12) / sqrt(h11 * h22)
)

# Export results to CSV
write.csv(results, "crypto_analysis_results.csv", row.names = FALSE)

# Create the break date as POSIXct
break_date_posix <- as.POSIXct(break_date)

# Theme setting for consistent look
theme_set(theme_minimal() +
            theme(
              plot.title = element_text(size = 12, face = "bold", hjust = 0.5),
              axis.title = element_text(size = 10),
              axis.text = element_text(size = 9),
              legend.position = "bottom",
              legend.title = element_blank(),
              panel.grid.minor = element_line(color = "gray90"),
              panel.grid.major = element_line(color = "gray85")
            ))

# 1. Returns Plot
p1 <- ggplot(results, aes(x = Date)) +
  geom_line(aes(y = BTC_Returns, color = "BTC"), alpha = 0.7) +
  geom_line(aes(y = ETH_Returns, color = "ETH"), alpha = 0.7) +
  geom_vline(xintercept = as.numeric(break_date_posix), 
             linetype = "dashed", color = "red", size = 1) +
  scale_color_manual(values = c("BTC" = "blue", "ETH" = "darkgreen")) +
  labs(title = "Cryptocurrency Returns",
       y = "Returns (%)",
       x = "") +
  scale_x_datetime(date_breaks = "3 months", 
                   date_labels = "%b %Y",
                   expand = c(0.02, 0.02)) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  annotate("text", x = break_date_posix, y = max(results$BTC_Returns, results$ETH_Returns),
           label = "Break Point", color = "red", vjust = -1)

# 2. Conditional Volatility Plot
p2 <- ggplot(results, aes(x = Date)) +
  geom_line(aes(y = BTC_Vol, color = "BTC"), alpha = 0.7) +
  geom_line(aes(y = ETH_Vol, color = "ETH"), alpha = 0.7) +
  geom_vline(xintercept = as.numeric(break_date_posix), 
             linetype = "dashed", color = "red", size = 1) +
  scale_color_manual(values = c("BTC" = "blue", "ETH" = "darkgreen")) +
  labs(title = "Conditional Volatility",
       y = "Volatility",
       x = "") +
  scale_x_datetime(date_breaks = "3 months", 
                   date_labels = "%b %Y",
                   expand = c(0.02, 0.02)) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# 3. Dynamic Correlation Plot
p3 <- ggplot(results, aes(x = Date)) +
  geom_line(aes(y = Correlation), color = "purple", alpha = 0.7) +
  geom_vline(xintercept = as.numeric(break_date_posix), 
             linetype = "dashed", color = "red", size = 1) +
  geom_hline(yintercept = mean(results$Correlation, na.rm = TRUE), 
             linetype = "dotted", color = "darkgrey") +
  labs(title = "Dynamic Correlation",
       y = "Correlation",
       x = "Date") +
  scale_y_continuous(limits = c(-1, 1)) +
  scale_x_datetime(date_breaks = "3 months", 
                   date_labels = "%b %Y",
                   expand = c(0.02, 0.02)) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# 4. Spillover Intensity Plot
p4 <- ggplot(results, aes(x = Date)) +
  geom_line(aes(y = Spillover), color = "darkred", alpha = 0.7) +
  geom_vline(xintercept = as.numeric(break_date_posix), 
             linetype = "dashed", color = "red", size = 1) +
  labs(title = "Spillover Intensity",
       y = "Intensity",
       x = "Date") +
  scale_x_datetime(date_breaks = "3 months", 
                   date_labels = "%b %Y",
                   expand = c(0.02, 0.02)) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Save plots with higher resolution
pdf("crypto_analysis_plots.pdf", width = 12, height = 10)
grid.arrange(p1, p2, p3, p4, ncol = 2, 
             top = textGrob("Cryptocurrency Market Dynamics with Structural Break",
                            gp = gpar(fontsize = 14, font = 2)))
dev.off()

# Also save as PNG for easier viewing
png("crypto_analysis_plots.png", width = 1200, height = 1000, res = 120)
grid.arrange(p1, p2, p3, p4, ncol = 2, 
             top = textGrob("Cryptocurrency Market Dynamics with Structural Break",
                            gp = gpar(fontsize = 14, font = 2)))
dev.off()

# Display plots in R
grid.arrange(p1, p2, p3, p4, ncol = 2, 
             top = textGrob("Cryptocurrency Market Dynamics with Structural Break",
                            gp = gpar(fontsize = 14, font = 2)))

# Structural break analysis
cat("\nStructural Break Analysis:\n")
pre_break <- results$Date < break_date_posix
post_break <- results$Date >= break_date_posix

# Period statistics calculation
calc_period_stats <- function(data, period) {
  data_subset <- data[period]
  data_subset <- data_subset[is.finite(data_subset)]
  
  if(length(data_subset) == 0) {
    return(c(Mean = NA, SD = NA, Min = NA, Max = NA))
  }
  
  c(Mean = mean(data_subset, na.rm = TRUE),
    SD = sd(data_subset, na.rm = TRUE),
    Min = min(data_subset, na.rm = TRUE),
    Max = max(data_subset, na.rm = TRUE))
}

# Calculate statistics for each measure
measures <- list(
  Correlation = results$Correlation,
  BTC_Vol = results$BTC_Vol,
  ETH_Vol = results$ETH_Vol,
  Spillover = results$Spillover
)

# Create summary statistics dataframe
summary_stats <- data.frame(
  Measure = character(),
  Period = character(),
  Mean = numeric(),
  SD = numeric(),
  Min = numeric(),
  Max = numeric(),
  stringsAsFactors = FALSE
)

# Calculate and store statistics
for(measure_name in names(measures)) {
  pre_stats <- calc_period_stats(measures[[measure_name]], pre_break)
  post_stats <- calc_period_stats(measures[[measure_name]], post_break)
  
  summary_stats <- rbind(summary_stats,
                         data.frame(
                           Measure = measure_name,
                           Period = "Pre-Break",
                           t(pre_stats)
                         ))
  summary_stats <- rbind(summary_stats,
                         data.frame(
                           Measure = measure_name,
                           Period = "Post-Break",
                           t(post_stats)
                         ))
}

# Export summary statistics
write.csv(summary_stats, "summary_statistics.csv", row.names = FALSE)

# Print statistics for each measure
for(measure_name in names(measures)) {
  pre_stats <- calc_period_stats(measures[[measure_name]], pre_break)
  post_stats <- calc_period_stats(measures[[measure_name]], post_break)
  
  if(!any(is.na(pre_stats)) && !any(is.na(post_stats))) {
    cat(sprintf("\n%s:\n", measure_name))
    cat("Pre-break  - ")
    cat(sprintf("Mean: %.3f, SD: %.3f, Range: [%.3f, %.3f]\n",
                pre_stats["Mean"], pre_stats["SD"], 
                pre_stats["Min"], pre_stats["Max"]))
    
    cat("Post-break - ")
    cat(sprintf("Mean: %.3f, SD: %.3f, Range: [%.3f, %.3f]\n",
                post_stats["Mean"], post_stats["SD"], 
                post_stats["Min"], post_stats["Max"]))
    
    # Calculate changes with handling for zero denominators
    if(abs(pre_stats["Mean"]) > 1e-10) {
      mean_change <- (post_stats["Mean"] / pre_stats["Mean"] - 1) * 100
      vol_change <- (post_stats["SD"] / pre_stats["SD"] - 1) * 100
      cat(sprintf("Changes - Mean: %.1f%%, Volatility: %.1f%%\n",
                  mean_change, vol_change))
    }
  }
}

# Calculate persistence measures
calc_persistence <- function(A, G, A_star = NULL, G_star = NULL) {
  if(is.null(A_star) && is.null(G_star)) {
    M <- kronecker(A, A) + kronecker(G, G)
  } else {
    M <- kronecker(A + A_star, A + A_star) + kronecker(G + G_star, G + G_star)
  }
  max(abs(eigen(M)$values))
}

pre_pers <- calc_persistence(bekk_results$A, bekk_results$G)
post_pers <- calc_persistence(bekk_results$A, bekk_results$G, 
                              bekk_results$A_star, bekk_results$G_star)

cat("\nPersistence Analysis:\n")
cat(sprintf("Pre-break:  %.3f %s\n", pre_pers,
            ifelse(pre_pers < 1, "(Stationary)", "(Non-stationary)")))
cat(sprintf("Post-break: %.3f %s\n", post_pers,
            ifelse(post_pers < 1, "(Stationary)", "(Non-stationary)")))

cat("\nAnalysis complete.\n")