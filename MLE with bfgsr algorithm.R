# *********************************************************************************************************
# Course: Econometrics
# Topic : ML estimation with the hand-written Broyden-Fletcher-Goldfarb-Shanno wiht restarts (BFGSR) 
#         algorithm 
# Instructor: Abdulbaki Bilgic
# *********************************************************************************************************
rm(list=ls()); cat("\f"); graphics.off()
# *********************************************************************************************************
bfgsr_regression <- function(data, yvar, xvars, 
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
  
  # Gradient function
  grad_fun <- function(theta) {
    b     <- theta[1:(length(theta)-1)]
    l     <- theta[length(theta)]
    sigma <- exp(l)
    e     <- as.vector(y - X %*% b)
    g_b   <- as.vector(t(X) %*% e) / (sigma^2)
    g_l   <- - n + sum(e^2) / (sigma^2)
    return(c(g_b, g_l))
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
  
  # BFGS update function
  bfgs_update <- function(H, delta, gamma) {
    rho   <- 1 / as.numeric(t(gamma) %*% delta)
    I     <- diag(length(delta))
    
    term1 <- I - rho * delta %*% t(gamma)
    term2 <- I - rho * gamma %*% t(delta)
    
    H_new <- term1 %*% H %*% term2 + rho * delta %*% t(delta)
    return(H_new)
  }
  
  # Line search function
  line_search <- function(theta, direction, grad, ll_current, alpha = 1, beta = 0.5, sigma = 1e-4) {
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
  
  # Single BFGS run
  bfgs_run <- function(theta_start, max_iter_run, tol_run, verbose_run) {
    H             <- diag(length(theta_start))
    theta_current <- theta_start
    g_current     <- grad_fun(theta_current)
    ll_current    <- ll_fun(theta_current)
    
    for (iter in 1:max_iter_run) {
      direction <- solve(H, g_current)
      alpha     <- line_search(theta_current, direction, g_current, ll_current)
      theta_new <- theta_current + alpha * direction
      g_new     <- grad_fun(theta_new)
      ll_new    <- ll_fun(theta_new)
      
      delta <- theta_new - theta_current
      gamma <- g_new - g_current
      
      if (sum(delta * gamma) > 1e-10) {
        H <- bfgs_update(H, delta, gamma)
      }
      
      if (verbose_run) {
        cat(sprintf("  Iter %d: logLik = %.6f\n", iter, ll_new))
      }
      
      if (sqrt(sum((theta_new - theta_current)^2)) < tol_run) {
        return(list(theta = theta_new, ll = ll_new, converged = TRUE, iterations = iter))
      }
      
      theta_current <- theta_new
      g_current     <- g_new
      ll_current    <- ll_new
    }
    
    return(list(theta = theta_current, ll = ll_current, converged = FALSE, iterations = max_iter_run))
  }
  
  # BFGS with restarts optimization
  cat("BFGS with Restarts Optimization\n")
  cat("================================\n")
  
  best_theta <- theta
  best_ll    <- ll_fun(theta)
  total_iter <- 0
  n_restarts <- 3
  
  for(restart in 1:n_restarts) {
    cat(sprintf("Restart %d:\n", restart))
    
    if (restart > 1) {
      # Random perturbation for subsequent restarts
      theta_start <- best_theta + rnorm(length(theta), 0, 0.1)
    } else {
      theta_start <- best_theta
    }
    
    result     <- bfgs_run(theta_start, max_iter %/% n_restarts, tol, verbose)
    total_iter <- total_iter + result$iterations
    
    if (result$ll > best_ll) {
      best_ll    <- result$ll
      best_theta <- result$theta
      cat(sprintf("  New best solution found: logLik = %.6f\n", best_ll))
    }
    
    if (result$converged) {
      cat(sprintf("  Converged after %d iterations\n", result$iterations))
      break
    }
  }
  
  theta      <- best_theta
  final_iter <- total_iter
  
  cat(sprintf("\nBest solution after %d restarts\n", n_restarts))
  
  # Final estimates and results
  b_hat         <- theta[1:(length(theta)-1)]
  log_sigma_hat <- theta[length(theta)]
  sigma_hat     <- exp(log_sigma_hat)
  
  H_final  <- hess_fun(theta)
  vcov_mat <- -solve(H_final)
  
  se_b <- sqrt(diag(vcov_mat))[1:length(b_hat)]
  se_log_sigma <- sqrt(diag(vcov_mat)[length(theta)])
  se_sigma     <- sigma_hat * se_log_sigma
  
  estimates  <- c(b_hat, sigma_hat)
  std_errors <- c(se_b, se_sigma)
  t_values   <- estimates / std_errors
  
  df       <- n - k
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
  
  row_names <- c(names(b_hat), "sigma")
  rownames(results) <- row_names
  
  cat(sprintf("\nFinal Results:\n"))
  cat(sprintf("Total Iterations: %d\n", final_iter))
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
    method       = "BFGS with Restarts"
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
  result <- bfgsr_regression(data, "y", c("X1","X2"), verbose = FALSE)
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