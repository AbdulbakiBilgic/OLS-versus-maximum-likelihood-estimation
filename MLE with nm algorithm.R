# *********************************************************************************************************
# Course: Econometrics
# Topic : ML estimation with the hand-written Nelder-Mead (NM) algorithm 
# Instructor: Abdulbaki Bilgic
# *********************************************************************************************************
rm(list=ls()); cat("\f"); graphics.off()
# *********************************************************************************************************
nelder_mead_regression <- function(data, yvar, xvars, 
                                   tol = 1e-6, max_iter = 100, 
                                   start_fraction = 1, verbose = TRUE) {
  # Prepare data
  y <- as.matrix(data[[yvar]])
  X <- as.matrix(cbind(1, data[, xvars]))
  colnames(X) <- c("Intercept", xvars)
  n <- nrow(X)
  k <- ncol(X)
  
  # OLS starting values
  beta_ols     <- solve(t(X) %*% X, t(X) %*% y)
  sigma_ols    <- sqrt(sum((y - X %*% beta_ols)^2) / n)
  theta        <- c(as.vector(start_fraction * beta_ols), log(sigma_ols))
  param_names  <- c(colnames(X), "logSigma")
  names(theta) <- param_names
  
  # Log-likelihood function
  ll_fun <- function(theta) {
    b     <- theta[1:(length(theta)-1)]
    l     <- theta[length(theta)]
    sigma <- exp(l)
    e     <- y - X %*% b
    logl  <- -0.5 * n * log(2 * pi * sigma^2) - 0.5 * sum(e^2) / (sigma^2)
    return(as.numeric(logl))
  }
  
  # Hessian matrix (for final variance calculation)
  hess_fun <- function(theta) {
    b     <- theta[1:(length(theta)-1)]
    l     <- theta[length(theta)]
    sigma <- exp(l)
    e     <- as.vector(y - X %*% b)
    H_bb  <- - (1 / (sigma^2)) * (t(X) %*% X)
    H_bl  <- - 2 * (1 / (sigma^2)) * (t(X) %*% e)
    H_lb  <- t(H_bl)
    H_ll  <- -2 * sum(e^2) / (sigma^2)
    
    H     <- matrix(0, nrow = length(theta), ncol = length(theta))
    H[1:(length(theta)-1), 1:(length(theta)-1)] <- H_bb
    H[1:(length(theta)-1), length(theta)]       <- H_bl
    H[length(theta), 1:(length(theta)-1)]       <- H_lb
    H[length(theta), length(theta)]             <- H_ll
    return(H)
  }
  
  # Nelder-Mead optimization
  cat("Nelder-Mead Optimization\n")
  cat("========================\n")
  
  p <- length(theta)
  
  # Create initial simplex
  simplex <- matrix(rep(theta, each = p + 1), nrow = p + 1, byrow = TRUE)
  for(i in 2:(p + 1)) {
    simplex[i, i - 1] <- simplex[i, i - 1] + 0.1  # Perturbation
  }
  
  final_iter <- max_iter
  
  for (iter in 1:max_iter) {
    # Calculate and sort values
    values    <- apply(simplex, 1, ll_fun)
    order_idx <- order(values, decreasing = TRUE)
    simplex   <- simplex[order_idx, ]
    values    <- values[order_idx]
    
    best         <- simplex[1, ]
    worst        <- simplex[p + 1, ]
    second_worst <- simplex[p, ]
    
    # Centroid (average of all points except worst)
    centroid <- colMeans(simplex[1:p, ])
    
    # Reflection
    reflected <- centroid + 1.0 * (centroid - worst)
    val_reflected <- ll_fun(reflected)
    
    if (val_reflected > values[p]) {
      # Expansion
      if (val_reflected > values[1]) {
        expanded     <- centroid + 2.0 * (centroid - worst)
        val_expanded <- ll_fun(expanded)
        if (val_expanded > val_reflected) {
          simplex[p + 1, ] <- expanded
        } else {
          simplex[p + 1, ] <- reflected
        }
      } else {
        simplex[p + 1, ] <- reflected
      }
    } else {
      # Contraction
      if (val_reflected >= values[p + 1]) {
        # Outside contraction
        contracted <- centroid + 0.5 * (centroid - worst)
      } else {
        # Inside contraction
        contracted <- centroid - 0.5 * (centroid - worst)
      }
      val_contracted <- ll_fun(contracted)
      
      if (val_contracted > values[p + 1]) {
        simplex[p + 1, ] <- contracted
      } else {
        # Shrink
        for(i in 2:(p + 1)) {
          simplex[i, ] <- best + 0.5 * (simplex[i, ] - best)
        }
      }
    }
    
    # Convergence check
    simplex_range <- max(apply(simplex, 2, function(x) diff(range(x))))
    best_val      <- values[1]
    
    if (verbose && iter %% 10 == 0) {
      cat(sprintf("Iter %d: logLik = %.6f, simplex range = %.6e\n", 
                  iter, best_val, simplex_range))
    }
    
    if (simplex_range < tol) {
      theta      <- best
      final_iter <- iter
      cat(sprintf("Converged after %d iterations\n", iter))
      break
    }
    
    theta <- best
    final_iter <- iter
  }
  
  if (iter == max_iter) {
    cat(sprintf("Maximum iterations (%d) reached\n", max_iter))
  }
  
  # Final estimates and results
  b_hat         <- theta[1:(length(theta)-1)]
  log_sigma_hat <- theta[length(theta)]
  sigma_hat     <- exp(log_sigma_hat)
  
  # Reassign parameter names (important for Nelder-Mead):
  if(is.null(names(b_hat))) {
    names(b_hat) <- param_names[1:(length(param_names)-1)]
  }
  
  H_final  <- hess_fun(theta)
  vcov_mat <- -solve(H_final)
  
  se_b         <- sqrt(diag(vcov_mat))[1:length(b_hat)]
  se_log_sigma <- sqrt(diag(vcov_mat)[length(theta)])
  se_sigma     <- sigma_hat * se_log_sigma
  
  estimates  <- c(b_hat, sigma_hat)
  std_errors <- c(se_b, se_sigma)
  t_values   <- estimates / std_errors
  
  df <- n - k
  p_values <- c(
    2 * pt(-abs(t_values[1:(length(t_values)-1)]), df = df),
    2 * pnorm(-abs(t_values[length(t_values)]))
  )
  
  breakpoints <- c(0, 0.001, 0.01, 0.05, 0.1, 1)
  labels      <- c("***", "**", "*", ".", " ")
  stars       <- cut(p_values, breaks = breakpoints, labels = labels, include.lowest = TRUE)
  stars       <- as.character(stars)
  
  results <- data.frame(
    Estimate    = estimates,
    `Std. Error`= std_errors,
    `t value`   = t_values,
    `Pr(>|t|)`  = p_values,
    Signif      = stars,
    check.names = FALSE
  )
  
  row_names         <- c(names(b_hat), "sigma")
  rownames(results) <- row_names
  
  cat(sprintf("\nFinal Results:\n"))
  cat(sprintf("Iterations: %d\n", final_iter))
  cat(sprintf("Log-Likelihood: %.4f\n", ll_fun(theta)))
  cat("\nCoefficients:\n")
  printCoefmat(results[,1:4], 
               signif.stars = TRUE, 
               digits       = 4, 
               P.values     = TRUE, 
               has.Pvalue   = TRUE)
  
  return(list(
    coefficients = results,
    sigma        = sigma_hat,
    logLik       = ll_fun(theta),
    iterations   = final_iter,
    vcov         = vcov_mat,
    method       = "Nelder-Mead"
  ))
}
# *********************************************************************************************************
# Test function
test_optimization_methods <- function() {
  set.seed(123)
  n <- 500
  data <- data.frame(
    X1 = rnorm(n),
    X2 = rnorm(n)
  )
  data$y <- 1 + 0.5 * data$X1 - 0.7 * data$X2 + rnorm(n, sd = 1)
  result <- nelder_mead_regression(data, "y", c("X1","X2"), max_iter = 1000, verbose = FALSE)
}

# Run tests
test_optimization_methods()
# *********************************************************************************************************
# Compare with lm()
set.seed(123)
n <- 500
data <- data.frame(
  X1 = rnorm(n),
  X2 = rnorm(n)
)
data$y  <- 1 + 0.5 * data$X1 - 0.7 * data$X2 + rnorm(n, sd = 1)
summary(lm(y ~ X1 + X2, data = data))
# *********************************************************************************************************