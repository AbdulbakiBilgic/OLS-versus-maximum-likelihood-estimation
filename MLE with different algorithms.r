# *********************************************************************************************************
# Course: Econometrics
# Topic : ML estimation with different hand written algorithms 
# Instructor: Abdulbaki Bilgic
# *********************************************************************************************************
rm(list=ls()); cat("\f"); graphics.off()
# *********************************************************************************************************
mle_with_different_algorithms <- function(data, 
                                          yvar, 
                                          xvars, 
                                          method         = "NR",                # NR, BFGS, BFGSR, BHHH, NM
                                          tol            = 1e-6, 
                                          max_iter       = 100, 
                                          start_fraction = 1, 
                                          verbose        = TRUE) {
  
  # Prepare data:
  y <- as.matrix(data[[yvar]])
  X <- as.matrix(cbind(1, data[, xvars]))
  colnames(X) <- c("Intercept", xvars)
  n <- nrow(X)
  k <- ncol(X)
  
  # OLS starting values:
  beta_ols     <- solve(t(X) %*% X, t(X) %*% y)
  sigma_ols    <- sqrt(sum((y - X %*% beta_ols)^2) / n)
  theta        <- c(as.vector(start_fraction * beta_ols), log(sigma_ols))
  param_names  <- c(colnames(X), "logSigma")
  names(theta) <- param_names
  
  # Log-likelihood function:
  ll_fun <- function(theta) {
    b     <- theta[1:(length(theta)-1)]
    l     <- theta[length(theta)]
    sigma <- exp(l)
    e     <- y - X %*% b
    logl  <- -0.5 * n * log(2 * pi * sigma^2) - 0.5 * sum(e^2) / (sigma^2)
    return(as.numeric(logl))
  }
  
  # Gradient function:
  grad_fun <- function(theta) {
    b     <- theta[1:(length(theta)-1)]
    l     <- theta[length(theta)]
    sigma <- exp(l)
    e     <- as.vector(y - X %*% b)
    g_b   <- as.vector(t(X) %*% e) / (sigma^2)
    g_l   <- - n + sum(e^2) / (sigma^2)
    return(c(g_b, g_l))
  }
  # Hessian matrix:
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
  
  # Individual gradients for BHHH (observation-based): works differently
  ind_grad_fun <- function(theta) {
    b     <- theta[1:(length(theta)-1)]
    l     <- theta[length(theta)]
    sigma <- exp(l)
    e     <- as.vector(y - X %*% b)
    
    # Gradient for each observation:
    grad_matrix <- matrix(0, n, length(theta))
    for(i in 1:n) {
      g_b             <- as.vector(X[i,] * e[i]) / (sigma^2)
      g_l             <- -1 + (e[i]^2) / (sigma^2)
      grad_matrix[i,] <- c(g_b, g_l)
    }
    return(grad_matrix)
  }
  
  # Helper functions for BFGS
  bfgs_update <- function(H, delta, gamma) {
    # BFGS formula: H_{k+1} = (I - ρδγ')H(I - ργδ') + ρδδ'
    # Where: δ = θ_{k+1} - θ_k, γ = g_{k+1} - g_k
    rho   <- 1 / as.numeric(t(gamma) %*% delta)
    I     <- diag(length(delta))
    
    term1 <- I - rho * delta %*% t(gamma)
    term2 <- I - rho * gamma %*% t(delta)
    
    H_new <- term1 %*% H %*% term2 + rho * delta %*% t(delta)
    return(H_new)
  }
  
  # Line search function
  line_search <- function(theta, direction, grad, ll_current, alpha = 1, beta = 0.5, sigma = 1e-4) {
    # Armijo condition: f(θ + αd) ≥ f(θ) + σαg'd
    max_ls_iter <- 20
    for(ls_iter in 1:max_ls_iter) {
      theta_new <- theta + alpha * direction
      ll_new    <- ll_fun(theta_new)
      
      if(ll_new >= ll_current + sigma * alpha * sum(grad * direction)) {
        return(alpha)
      }
      alpha <- beta * alpha
    }
    return(alpha)
  }
  
  # Optimization algorithms:
  final_iter <- 0                                                               # Initialize final_iter
  
  if(method == "NR") {
    # Newton-Raphson loop:
    for (iter in 1:max_iter) {
      g         <- grad_fun(theta)
      H         <- hess_fun(theta)
      step      <- solve(H, g)
      theta_new <- theta - step
      
      if (verbose) {
        cat(sprintf("Iter %d: logLik=%.6f, step norm=%.6e\n", 
                    iter, ll_fun(theta_new), sqrt(sum(step^2))))
      }
      
      if (sqrt(sum((theta_new - theta)^2)) < tol) {
        theta      <- theta_new
        final_iter <- iter
        break
      }
      theta      <- theta_new
      final_iter <- iter
    }
    
  } else if(method == "BFGS") {
    # BFGS algorithm (manual implementation):
    H             <- diag(length(theta))                                        # Initial Hessian approximation (identity matrix)
    theta_current <- theta
    g_current     <- grad_fun(theta_current)
    ll_current    <- ll_fun(theta_current)
    
    for (iter in 1:max_iter) {
      # Direction vector:
      direction <- solve(H, g_current)
      
      # Line search:
      alpha     <- line_search(theta_current, direction, g_current, ll_current)
      
      # Parameter update:
      theta_new <- theta_current + alpha * direction
      g_new     <- grad_fun(theta_new)
      ll_new    <- ll_fun(theta_new)
      
      # BFGS update:
      delta <- theta_new - theta_current
      gamma <- g_new - g_current
      
      if (sum(delta * gamma) > 1e-10) {                                         # Curvature condition
        H <- bfgs_update(H, delta, gamma)
      }
      
      if (verbose) {
        cat(sprintf("Iter %d: logLik=%.6f, step norm=%.6e, alpha=%.4f\n", 
                    iter, ll_new, sqrt(sum((theta_new - theta_current)^2)), alpha))
      }
      
      if (sqrt(sum((theta_new - theta_current)^2)) < tol) {
        theta_current <- theta_new
        final_iter    <- iter
        break
      }
      
      theta_current <- theta_new
      g_current     <- g_new
      ll_current    <- ll_new
      final_iter    <- iter
    }
    theta <- theta_current
    
  } else if(method == "BFGSR") {
    # BFGS with restarts
    best_theta <- theta
    best_ll    <- ll_fun(theta)
    total_iter <- 0
    
    for(restart in 1:3) {
      H             <- diag(length(theta))
      theta_current <- best_theta + rnorm(length(theta), 0, 0.1)                # Random perturbation
      g_current     <- grad_fun(theta_current)
      ll_current    <- ll_fun(theta_current)
      
      for (iter in 1:(max_iter %/% 3)) {
        direction <- solve(H, g_current)
        alpha     <- line_search(theta_current, direction, g_current, ll_current)
        theta_new <- theta_current + alpha * direction
        g_new     <- grad_fun(theta_new)
        ll_new    <- ll_fun(theta_new)
        
        delta     <- theta_new - theta_current
        gamma     <- g_new - g_current
        
        if (sum(delta * gamma) > 1e-10) {
          H <- bfgs_update(H, delta, gamma)
        }
        
        if (verbose && restart == 1) {
          cat(sprintf("Restart %d, Iter %d: logLik=%.6f\n", restart, iter, ll_new))
        }
        
        if (sqrt(sum((theta_new - theta_current)^2)) < tol) {
          total_iter <- total_iter + iter
          break
        }
        
        theta_current <- theta_new
        g_current     <- g_new
        ll_current    <- ll_new
        total_iter    <- total_iter + 1
      }
      
      if (ll_current > best_ll) {
        best_ll    <- ll_current
        best_theta <- theta_current
      }
    }
    theta      <- best_theta
    final_iter <- total_iter
    
  } else if(method == "BHHH") {
    # BHHH algorithm
    learning_rate <- 0.01
    for (iter in 1:max_iter) {
      g <- grad_fun(theta)
      G <- ind_grad_fun(theta)
      
      # BHHH Hessian approximation: G'G
      H_bhhh <- t(G) %*% G
      H_bhhh <- H_bhhh + diag(1e-8, ncol(H_bhhh))                               # Regularization
      
      step      <- solve(H_bhhh, g)
      theta_new <- theta + learning_rate * step
      
      if (verbose) {
        cat(sprintf("Iter %d: logLik=%.6f, step norm=%.6e\n", 
                    iter, ll_fun(theta_new), sqrt(sum(step^2))))
      }
      
      if (sqrt(sum((theta_new - theta)^2)) < tol) {
        theta      <- theta_new
        final_iter <- iter
        break
      }
      theta      <- theta_new
      final_iter <- iter
    }
    
  } else if(method == "NM") {
    # Nelder-Mead algorithm (manual implementation)
    p <- length(theta)
    
    # Create simplex:
    simplex <- matrix(rep(theta, each = p + 1), nrow = p + 1, byrow = TRUE)
    for(i in 2:(p + 1)) {
      simplex[i, i - 1] <- simplex[i, i - 1] + 0.1                              # Perturbation
    }
    
    final_iter <- max_iter                                                      # Default value
    
    for (iter in 1:max_iter) {
      # Calculate and sort values:
      values    <- apply(simplex, 1, ll_fun)
      order_idx <- order(values, decreasing = TRUE)
      simplex   <- simplex[order_idx, ]
      values    <- values[order_idx]
      
      best         <- simplex[1, ]
      worst        <- simplex[p + 1, ]
      second_worst <- simplex[p, ]
      
      # Centroid (average of all points except worst):
      centroid <- colMeans(simplex[1:p, ])
      
      # Reflection:
      reflected     <- centroid + 1.0 * (centroid - worst)
      val_reflected <- ll_fun(reflected)
      
      if (val_reflected > values[p]) {
        # Expansion:
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
        # Contraction:
        if (val_reflected >= values[p + 1]) {
          # Outside contraction:
          contracted <- centroid + 0.5 * (centroid - worst)
        } else {
          # Inside contraction:
          contracted <- centroid - 0.5 * (centroid - worst)
        }
        val_contracted <- ll_fun(contracted)
        
        if (val_contracted > values[p + 1]) {
          simplex[p + 1, ] <- contracted
        } else {
          # Shrink:
          for(i in 2:(p + 1)) {
            simplex[i, ] <- best + 0.5 * (simplex[i, ] - best)
          }
        }
      }
      
      # Convergence check:
      simplex_range <- max(apply(simplex, 2, function(x) diff(range(x))))
      best_val      <- values[1]
      
      if (verbose && iter %% 10 == 0) {
        cat(sprintf("Iter %d: logLik=%.6f, simplex range=%.6e\n", 
                    iter, best_val, simplex_range))
      }
      
      if (simplex_range < tol) {
        theta      <- best
        final_iter <- iter
        break
      }
      
      theta      <- best
      final_iter <- iter
    }
    
  } else {
    stop("Invalid method. Available methods: NR, BFGS, BFGSR, BHHH, NM")
  }
  
  # Final estimates (common for all methods):
  b_hat         <- theta[1:(length(theta)-1)]
  log_sigma_hat <- theta[length(theta)]
  sigma_hat     <- exp(log_sigma_hat)
  
  # Reassign parameter names (especially needed for NM):
  if(is.null(names(b_hat))) {
    names(b_hat) <- param_names[1:(length(param_names)-1)]
  }
  
  # Variance-covariance matrix (Hessian inverse):
  H_final  <- hess_fun(theta)
  vcov_mat <- -solve(H_final)
  
  # Standard errors:
  se_b         <- sqrt(diag(vcov_mat))[1:length(b_hat)]
  se_log_sigma <- sqrt(diag(vcov_mat)[length(theta)])
  se_sigma     <- sigma_hat * se_log_sigma
  
  # t/z values and p-values:
  estimates  <- c(b_hat, sigma_hat)
  std_errors <- c(se_b, se_sigma)
  t_values   <- estimates / std_errors
  
  df <- n - k
  p_values <- c(
    2 * pt(-abs(t_values[1:(length(t_values)-1)]), df = df),
    2 * pnorm(-abs(t_values[length(t_values)]))
  )
  
  # Significance stars:
  breakpoints <- c(0, 0.001, 0.01, 0.05, 0.1, 1)
  labels      <- c("***", "**", "*", ".", " ")
  stars       <- cut(p_values, breaks = breakpoints, labels = labels, include.lowest = TRUE)
  stars       <- as.character(stars)
  
  # Results table - safe row name assignment:
  results <- data.frame(
    Estimate    = estimates,
    `Std. Error`= std_errors,
    `t value`   = t_values,
    `Pr(>|t|)`  = p_values,
    Signif      = stars,
    check.names = FALSE
  )
  
  # Safe row name assignment:
  row_names <- c(names(b_hat), "sigma")
  if(length(row_names) == nrow(results)) {
    rownames(results) <- row_names
  } else {
    # If row names length doesn't match, use generic names:
    rownames(results) <- c(paste0("beta", 0:(length(b_hat)-1)), "sigma")
  }
  cat(sprintf("\nOptimization Method: %s\n", method))
  cat(sprintf("Iterations: %d\n", final_iter))
  cat(sprintf("Log-Likelihood: %.4f\n", ll_fun(theta)))
  cat("\nCoefficients:\n")
  printCoefmat(results[,1:4], 
               signif.stars = TRUE, 
               digits       = 4, 
               P.values     = TRUE, 
               has.Pvalue   = TRUE, 
               na.print     = "")
  
  return(list(
    coefficients = results,
    sigma        = sigma_hat,
    logLik       = ll_fun(theta),
    iterations   = final_iter,
    vcov         = vcov_mat,
    method       = method
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
  data$y  <- 1 + 0.5 * data$X1 - 0.7 * data$X2 + rnorm(n, sd = 1)
  
  methods <- c("NR", "BFGS", "BFGSR", "BHHH", "NM")
  
  for(method in methods) {
    cat(sprintf("\n%s Testing %s %s\n", paste(rep("=", 20), collapse=""), 
                method, paste(rep("=", 20), 
                              collapse="")))
    result <- mle_with_different_algorithms(data, "y", c("X1","X2"), 
                                            method = method, verbose = FALSE)
  }
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
