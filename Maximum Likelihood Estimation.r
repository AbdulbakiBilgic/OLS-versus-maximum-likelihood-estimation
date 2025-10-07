#===============================================================================
# Course: Econometrics
# Topic : ML estimation
# Instructor: Abdulbaki Bilgic
#===============================================================================
rm(list=ls()); cat("\f"); graphics.off()
library(maxLik); library(corrr)
#===============================================================================
# Do not change here:
ml_regression <- function(data, 
                          yvar, 
                          xvars, 
                          method         = "NR", 
                          start_fraction = 1,
                          print.level    = 0) {
  #===================================
  # Prepare y and X from data:
  y <- as.matrix(data[[yvar]])
  X <- as.matrix(cbind(1, data[, xvars]))
  colnames(X) <- c("Intercept", xvars)
  #===================================
  # OLS estimates for starting values:
  beta_ols  <- solve(t(X) %*% X, t(X) %*% y)
  sigma_ols <- sqrt(sum((y - X %*% beta_ols)^2) / length(y))
  #===================================
  # 1) Log-likelihood:
  ll_fun <- function(theta) {
    b     <- theta[1:(length(theta)-1)]
    l     <- theta[length(theta)]
    sigma <- exp(l)
    e     <- y - X %*% b
    logl  <- -0.5 * length(y) * log(2 * pi * sigma^2) - 0.5 * sum(e^2) / (sigma^2)
    return(as.numeric(logl))
  }
  # 2) Gradient (or score):
  grad_fun <- function(theta) {
    b     <- theta[1:(length(theta)-1)]
    l     <- theta[length(theta)]
    sigma <- exp(l)
    e     <- as.vector(y - X %*% b)
    g_b   <- as.vector(t(X) %*% e) / (sigma^2)
    g_l   <- - length(y) + sum(e^2) / (sigma^2)
    return(c(g_b, g_l))
  }
  # 3) Hessian matrix:
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
  #===================================
  # Starting values (OLS + log sigma):
  startv        <- c(as.vector(start_fraction * beta_ols), log(sigma_ols))
  names(startv) <- c(colnames(X), "logSigma")
  #===================================
  # Maximum Likelihood Estimation:
  mle_fit <- maxLik(logLik = ll_fun,
                    grad   = grad_fun,
                    hess   = hess_fun,
                    start  = startv,
                    method = method,
                    print.level = 3)
  #===================================
  # Return results:
  return(summary(mle_fit))
}
#===============================================================================
# Example usage:
set.seed(123)
n  <- 1000
data <- data.frame(
  X1 = rnorm(n),
  X2 = rnorm(n)
)
data$y <- 1 + 0.5 * data$X1 - 0.7 * data$X2 + rnorm(n, sd = 1)
# Get some diagnostics about your data:
print(psych::describe(data))                                                    # get your descriptive statistics (pay attention to the mean and sd of y)
cor_matrix <- cor(data)                                                         # get information about correlations among all variables
print(round(cor_matrix, 3))
#===============================================================================
# Obtain maximum likelihood (ML) results:
yvar  <- c("y")
xvars <- c("X1", "X2")
result <- ml_regression(data, 
                        yvar  = yvar, 
                        xvars = xvars, 
                        method= "BFGS",
                        start_fraction = 0.05,                                  # Change start_fraction freely, or set it to 1 if you want OLS estimates as starting values
                        print.level    = 3)
print(result)
#===============================================================================
# Compare with lm() results:
ols_fit  <- lm(y ~ ., data = data)
print("lm() summary:")
summary(ols_fit)
#===============================================================================

#===============================================================================

