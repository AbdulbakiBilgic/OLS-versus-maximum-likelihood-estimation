# *********************************************************************************************************
# Course    : Econometrics
# Topic     : Optimization techniques
# Instructor: Abdulbaki Bilgic
# *********************************************************************************************************
# mle optimization methods.R
# Full R script: MLE optimizer suite with NR, BFGS, BFGSR, BHHH, NM
# Includes numerical derivative wrappers, logit, poisson, gaussian examples, and comparison with glm/lm
# *********************************************************************************************************
# rm(list=ls()); cat("\f"); graphics.off()
# *********************************************************************************************************
if(!requireNamespace("numDeriv", quietly = TRUE)) {
  stop("Please install 'numDeriv' package: install.packages('numDeriv')")
}
# *********************************************************************************************************
# Armijo line search for maximization
armijo_line_search <- function(ll_fun, grad_fun, theta, direction,
                               alpha_init = 1, rho = 0.5, c = 1e-4, max_iter = 50) {
  alpha <- alpha_init
  ll_current <- ll_fun(theta)
  g <- grad_fun(theta)
  for (i in seq_len(max_iter)) {
    theta_new <- theta + alpha * direction
    ll_new <- ll_fun(theta_new)
    if (ll_new >= ll_current + c * alpha * sum(g * direction)) return(alpha)
    alpha <- alpha * rho
  }
  return(alpha)
}
# *********************************************************************************************************
# BFGS inverse Hessian update:
bfgs_inverse_update <- function(B, s, y) {
  sty   <- as.numeric(crossprod(s, y))
  if (sty <= 1e-12) return(B)
  rho   <- 1 / sty
  I     <- diag(length(s))
  V     <- I - rho * s %*% t(y)
  B_new <- V %*% B %*% t(V) + rho * s %*% t(s)
  B_new
}
# *********************************************************************************************************
# Numerical gradient and Hessian wrapper
make_grad_hess_wrappers <- function(ll_fun, theta) {
  grad_fun <- function(theta) numDeriv::grad(ll_fun, theta)
  hess_fun <- function(theta) numDeriv::hessian(ll_fun, theta)
  # Individual gradients for BHHH
  ind_grad_fun <- function(theta) {
    n <- length(ll_fun(theta, individual = TRUE))
    p <- length(theta)
    G <- matrix(0, n, p)
    for(i in 1:n) G[i,] <- numDeriv::grad(function(t) ll_fun(t, individual = TRUE)[i], theta)
    G
  }
  attr(grad_fun, "ind_grad") <- ind_grad_fun
  list(grad_fun = grad_fun, hess_fun = hess_fun)
}
# *********************************************************************************************************
# Summary table for MLE
mle_summary_table <- function(theta, hess_fun, n_obs, param_names = NULL) {
  vcov_mat  <- -solve(hess_fun(theta))
  std_errors<- sqrt(diag(vcov_mat))
  t_values  <- theta / std_errors
  p_values  <- 2 * pnorm(-abs(t_values))
  
  p_display <- ifelse(p_values < 2e-16, "<2e-16", format(p_values, digits = 4, scientific = TRUE))
  # Significance stars:
  breakpoints <- c(0, 0.001, 0.01, 0.05, 0.1, 1)
  labels      <- c("***", "**", "*", ".", " ")
  stars       <- cut(p_values, breaks = breakpoints, labels = labels, include.lowest = TRUE)
  stars       <- as.character(stars)
  
  results <- data.frame(
    Estimate    = theta,
    `Std. Error`= std_errors,
    `t value`   = t_values,
    `Pr(>|t|)`  = p_display,
    ` `         = stars,
    check.names = FALSE
  )
  if(!is.null(param_names)) rownames(results) <- param_names
  results
}
# *********************************************************************************************************
# Optimizer methods
optimizer_methods <- function(theta_init, ll_fun, grad_fun, hess_fun,
                              method  = c("NR","BFGS","BFGSR","BHHH","NM"),
                              tol      = 1e-6, 
                              max_iter = 200, 
                              verbose  = TRUE) {
  method     <- match.arg(method)
  theta      <- as.numeric(theta_init)
  p          <- length(theta)
  final_iter <- 0
  
  if (method == "NR") {
    for (iter in 1:max_iter) {
      g         <- grad_fun(theta)
      H         <- hess_fun(theta)
      step      <- tryCatch(solve(H, g), error = function(e) solve(H + diag(1e-8, p), g))
      theta_new <- theta - step
      if (verbose) cat(sprintf("[NR] Iter %d | logLik=%.6f | stepNorm=%.3e\n", iter, 
                               ll_fun(theta_new), sqrt(sum(step^2))))
      if (sqrt(sum((theta_new - theta)^2)) < tol) { theta <- theta_new; final_iter <- iter; break }
      theta <- theta_new; final_iter <- iter
    }
  } else if (method == "BFGS") {
    B             <- diag(p)
    theta_current <- theta
    g_current     <- grad_fun(theta_current)
    for (iter in 1:max_iter) {
      direction <- as.numeric(B %*% g_current)
      alpha     <- armijo_line_search(ll_fun, grad_fun, theta_current, direction)
      theta_new <- theta_current + alpha * direction
      g_new     <- grad_fun(theta_new)
      s         <- theta_new - theta_current
      y         <- g_new - g_current
      if (sum(s * y) > 1e-12) B <- bfgs_inverse_update(B, s, y)
      if (verbose) cat(sprintf("[BFGS] Iter %d | logLik=%.6f | stepNorm=%.3e | alpha=%.4f\n", iter, 
                               ll_fun(theta_new), sqrt(sum((theta_new - theta_current)^2)), alpha))
      if (sqrt(sum((theta_new - theta_current)^2)) < tol) { theta_current <- theta_new; final_iter <- iter; break }
      theta_current <- theta_new; 
      g_current     <- g_new; 
      final_iter    <- iter
    }
    theta <- theta_current
  } else if (method == "BFGSR") {
    best_theta <- theta
    best_ll    <- ll_fun(theta)
    total_iter <- 0
    n_restarts <- 3
    for (r in 1:n_restarts) {
      theta_current <- best_theta + rnorm(p, 0, 0.1)
      B             <- diag(p)
      g_current     <- grad_fun(theta_current)
      for (iter in 1:ceiling(max_iter / n_restarts)) {
        direction <- as.numeric(B %*% g_current)
        alpha     <- armijo_line_search(ll_fun, grad_fun, theta_current, direction)
        theta_new <- theta_current + alpha * direction
        g_new     <- grad_fun(theta_new)
        s         <- theta_new - theta_current
        y         <- g_new - g_current
        if (sum(s * y) > 1e-12) B <- bfgs_inverse_update(B, s, y)
        theta_current <- theta_new; g_current <- g_new
        total_iter    <- total_iter + 1
      }
      ll_current <- ll_fun(theta_current)
      if (ll_current > best_ll) { best_ll <- ll_current; best_theta <- theta_current }
    }
    theta <- best_theta; final_iter <- total_iter
  } else if (method == "BHHH") {
    if (is.null(attr(grad_fun, "ind_grad"))) stop("BHHH requires an 'ind_grad' attribute on grad_fun (function returning n x p matrix).")
    learning_rate <- 0.01
    for (iter in 1:max_iter) {
      g         <- grad_fun(theta)
      G         <- attr(grad_fun, "ind_grad")(theta)
      H_bhhh    <- t(G) %*% G + diag(1e-8, p)
      step      <- tryCatch(solve(H_bhhh, g), error = function(e) solve(H_bhhh + diag(1e-8, p), g))
      theta_new <- theta + learning_rate * step
      if (verbose) cat(sprintf("[BHHH] Iter %d | logLik=%.6f | stepNorm=%.3e\n", iter, 
                               ll_fun(theta_new), sqrt(sum(step^2))))
      if (sqrt(sum((theta_new - theta)^2)) < tol) { theta <- theta_new; final_iter <- iter; break }
      theta <- theta_new; final_iter <- iter
    }
  } else if (method == "NM") {
    p       <- length(theta)
    simplex <- matrix(rep(theta, each = p+1), nrow = p+1, byrow = TRUE)
    for (i in 2:(p+1)) simplex[i, i-1] <- simplex[i, i-1] + 0.1
    for (iter in 1:max_iter) {
      vals          <- apply(simplex, 1, ll_fun)
      ord           <- order(vals, decreasing = TRUE)
      simplex       <- simplex[ord, , drop = FALSE]
      vals          <- vals[ord]
      best          <- simplex[1,]; worst <- simplex[p+1,]
      centroid      <- colMeans(simplex[1:p, , drop = FALSE])
      reflected     <- centroid + 1.0 * (centroid - worst)
      val_reflected <- ll_fun(reflected)
      if (val_reflected > vals[p]) {
        if (val_reflected > vals[1]) {
          expanded     <- centroid + 2.0 * (centroid - worst)
          val_expanded <- ll_fun(expanded)
          if (val_expanded > val_reflected) simplex[p+1, ] <- expanded else simplex[p+1, ] <- reflected
        } else simplex[p+1, ] <- reflected
      } else {
        if (val_reflected >= vals[p+1]) {
          contracted <- centroid + 0.5 * (centroid - worst)
        } else {
          contracted <- centroid - 0.5 * (centroid - worst)
        }
        val_contracted <- ll_fun(contracted)
        if (val_contracted > vals[p+1]) simplex[p+1, ] <- contracted
        else { for (i in 2:(p+1)) simplex[i, ] <- best + 0.5 * (simplex[i, ] - best) }
      }
      simplex_range <- max(apply(simplex, 2, function(x) diff(range(x))))
      if (iter %% 10 == 0 && verbose) cat(sprintf("[NM] Iter %d | best logLik=%.6f | simplexRange=%.3e\n", 
                                                  iter, vals[1], simplex_range))
      if (simplex_range < tol) { theta <- simplex[1, ]; final_iter <- iter; break }
      theta <- simplex[1, ]; final_iter <- iter
    }
  }
  list(theta = theta, logLik = ll_fun(theta), iterations = final_iter, method = method)
}
# *********************************************************************************************************