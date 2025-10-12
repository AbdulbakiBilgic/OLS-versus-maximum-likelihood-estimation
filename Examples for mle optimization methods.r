# *********************************************************************************************************
# Different examples for mle optimization methods:
# *********************************************************************************************************
rm(list=ls()); cat("\f"); graphics.off()
source('F:/R-Programs/mle optimization methods.r')                              # save "mle optimization methods.r" and then call it
# *********************************************************************************************************
# Examples: Gaussian, Poisson, Logistic
set.seed(123)
n    <- 500
data <- data.frame(X1 = rnorm(n), X2 = rnorm(n))
# *********************************************************************************************************
# Gaussian Example
data$y_gauss <- 1 + 0.5 * data$X1 - 0.7 * data$X2 + rnorm(n, 0, 1)
ll_gauss <- function(theta) {
  beta   <- theta[1:(length(theta)-1)]
  sigma  <- exp(theta[length(theta)])
  mu     <- beta[1] + beta[2]*data$X1 + beta[3]*data$X2
  sum(dnorm(data$y_gauss, mean=mu, sd=sigma, log=TRUE))                         # you can set up your own log-likelihood function by hand if you wish
}

Chosen_method <- "BFGS"                                                         # you can call a different algorithm method if you wish
wrappers      <- make_grad_hess_wrappers(ll_gauss, rep(0,4))
res_gauss     <- optimizer_methods(rep(0,4), 
                                   ll_gauss, 
                                   wrappers$grad_fun, 
                                   wrappers$hess_fun, 
                                   method  = Chosen_method, 
                                   verbose = FALSE)

cat(sprintf("Gaussian regression MLE (%s):\n", Chosen_method))
print(mle_summary_table(res_gauss$theta, wrappers$hess_fun, n, c("Intercept","X1","X2","logSigma")))
cat("\nCheck with lm():\n")

print(summary(lm(y_gauss ~ X1 + X2, data=data)))
# *********************************************************************************************************
# Poisson Example
data$y_pois <- rpois(n, lambda=exp(0.2 + 0.3*data$X1 - 0.5*data$X2))
ll_pois  <- function(theta) {
  beta   <- theta
  lambda <- exp(beta[1] + beta[2]*data$X1 + beta[3]*data$X2)
  sum(dpois(data$y_pois, lambda=lambda, log=TRUE))                              # you can set up your own log-likelihood function by hand if you wish
}
wrappers <- make_grad_hess_wrappers(ll_pois, rep(0,3))
res_pois <- optimizer_methods(rep(0,3), ll_pois, wrappers$grad_fun, wrappers$hess_fun, 
                              method = Chosen_method,
                              verbose= FALSE)

cat(sprintf("Gaussian regression MLE (%s):\n", Chosen_method))
print(mle_summary_table(res_pois$theta, wrappers$hess_fun, n, c("Intercept","X1","X2")))
cat("\nCheck with glm(family=poisson):\n")

print(summary(glm(y_pois ~ X1 + X2, data=data, family=poisson)))
# *********************************************************************************************************
# Logistic Example
data$y_logit <- rbinom(n, size=1, prob=1/(1+exp(-( -0.5 + 0.7*data$X1 - 1*data$X2))))
ll_logit <- function(theta) {
  beta <- theta
  eta <- beta[1] + beta[2]*data$X1 + beta[3]*data$X2
  sum(data$y_logit * eta - log(1 + exp(eta)))                                   # you can set up your own log-likelihood function by hand if you wish
}
wrappers  <- make_grad_hess_wrappers(ll_logit, rep(0,3))
res_logit <- optimizer_methods(rep(0,3), ll_logit, wrappers$grad_fun, wrappers$hess_fun, 
                               method = Chosen_method, verbose = FALSE)

cat(sprintf("Gaussian regression MLE (%s):\n", Chosen_method))
print(mle_summary_table(res_logit$theta, wrappers$hess_fun, n, c("Intercept","X1","X2")))
cat("\nCheck with glm(family=binomial):\n")
print(summary(glm(y_logit ~ X1 + X2, data=data, family=binomial)))
# *********************************************************************************************************