#===============================================================================
# Course: Compare OLS with maxlik estimates
# Author: Abdulbaki Bilgic
#===============================================================================
rm(list = ls(all = TRUE))
library(maxLik)
#===============================================================================
# To construct data:
set.seed(123)
n <- 500
x <- rnorm(n, mean = 0, sd = 1)
y <- rnorm(n, mean = 500, sd = 10)   # sample data
X <- cbind(1, x)
colnames(X) <- c("Constant", "x1")
#===============================================================================
# 1) Closed-form OLS starting values:
beta_ols_closed  <- solve(t(X) %*% X) %*% t(X) %*% y                            # Obtain OLS estimates with matrix form
resid_ols        <- y - X %*% beta_ols_closed                                   # Obtaine OLS residuals
sigma_mle_closed <- sqrt(sum(resid_ols^2) / n)                                  # Construct sigma estimate

# OLS coefficients and sigma in a table
ols_results <- data.frame(
  Variable = c(colnames(X), "sigma"),
  Estimate = c(as.vector(beta_ols_closed), sigma_mle_closed)
)
print(ols_results)
#===============================================================================
# 2) Log-likelihood function:
ll_fun <- function(theta) {
  b     <- theta[1:2]
  l     <- theta[3]
  sigma <- exp(l)
  e     <- y - X %*% b
  logl  <- -0.5 * length(y) * log(2 * pi * sigma^2) - 0.5 * sum(e^2) / (sigma^2)
  return(as.numeric(logl))
}

# 3) Gradient (scores):
grad_fun <- function(theta) {
  b     <- theta[1:2]
  l     <- theta[3]
  sigma <- exp(l)
  e     <- as.vector(y - X %*% b)
  g_b   <- as.vector(t(X) %*% e) / (sigma^2)   # length 2
  g_l   <- - length(y) + sum(e^2) / (sigma^2)  # scalar
  return(c(g_b, g_l))
}

# 4) Hessian matrix:
hess_fun <- function(theta) {
  b     <- theta[1:2]
  l     <- theta[3]
  sigma <- exp(l)
  e     <- as.vector(y - X %*% b)
  # constructing each block:
  H_bb  <- - (1 / (sigma^2)) * (t(X) %*% X)                                     # 2x2
  H_bl  <- - 2 * (1 / (sigma^2)) * (t(X) %*% e)                                 # 2x1
  H_lb  <- t(H_bl)                                                              # 1x2
  H_ll  <- -2 * sum(e^2) / (sigma^2)                                            # scalar
  # build full 3x3
  H           <- matrix(0, nrow = 3, ncol = 3)
  H[1:2, 1:2] <- H_bb
  H[1:2, 3]   <- H_bl
  H[3, 1:2]   <- H_lb
  H[3, 3]     <- H_ll
  return(H)
}

# 5) Starting values (OLS coefficients + log sigma MLE):
startv        <- c(as.vector(0.5*beta_ols_closed), log(sigma_mle_closed))
names(startv) <- c("Constant", "x1", "logSigma")

# 6) maxLik estimation:
fit <- maxLik(logLik = ll_fun, 
              grad   = grad_fun, 
              hess   = hess_fun, 
              start  = startv, 
              method = "NR",                                                    # change your optimization choice here 
              print.level = 1)

# 7) Display maxLik results:
summary(fit)

# 8) Compare with lm() results:
lm.model <- lm(y ~ x)
summary(lm.model)
#===============================================================================