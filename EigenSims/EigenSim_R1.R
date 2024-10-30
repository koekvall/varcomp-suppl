# -------------------------------------------------
# set n
# -------------------------------------------------
uu <- as.numeric(Sys.getenv('SLURM_ARRAY_TASK_ID'))
params <-  expand.grid(n=c(seq(200, 2000, length=9)), h2 = c(.0001, .01, 0.1, 0.5), rho = c(0.1, 0.5, 0.95, 0.999))

n <- params[uu,1]
h2 <- params[uu,2]
rho <- params[uu,3]
savename <- paste("~/blue/lmmvar/Simulations/Eigen/Results_R1/Coverage_", uu, "_R1.RDS", sep="")

library(MASS)
library(profvis)
library(matrixcalc)
library(lmmvar)

# -------------------------------------------------
# simulate y with given Z
# -------------------------------------------------
simData_K1 <- function(sigma2_g, sigma2_e = 1, X, Ksqrt) {
  n <- dim(X)[1]
  p <- dim(X)[2]
  beta <- rep(0, p)
  e <- rnorm(n, mean = 0, sd = sqrt(sigma2_e))
  u <- Ksqrt%*%rnorm(n)
  y <- X %*% beta + sqrt(sigma2_g)*u + e
  return(y)
}

# -------------------------------------------------
# negative log-likelihood function for h2, s2p
# -------------------------------------------------
neg_loglik2 <- function(par, X, y, lambda) {
  n <- nrow(X)
  Sigma <- (par[1]*lambda + rep(1-par[1],n)) * par[2]
  Sigma_inv <- 1 / Sigma
  XSX <- crossprod(X, Sigma_inv * X)
  betahat <- solve(XSX, crossprod(X, Sigma_inv * y))
  XB <- crossprod(t(X), betahat)
  l <- sum(log(Sigma)) / 2 +
    determinant(XSX)$modulus[1] / 2 +
    sum((y-XB) * Sigma_inv * (y-XB)) / 2
  return(l)
}

# -----------------------------------
# true parameter values; beta = 0 
# ----------------------------------
s2p <- 1

# ----------------------------------
# Model setup
# ---------------------------------
set.seed(1)
p <- 5
X <- matrix(rnorm(n*p), nrow=n)
K <- matrix(0, nrow=n, ncol=n)
for(j in 1:n){
  for(k in 1:n){
    K[j,k] <- rho^abs(j-k)
  }
}
eoK <- eigen(K)
V <- eoK$vec
lambda <- eoK$val
Ksqrt <- V%*%diag(lambda^0.5)%*%t(V)
Xnew <- crossprod(V, X)

reps <- 1e4
lrt <- matrix(0, nrow = reps, ncol = 3)
score <- matrix(0, nrow = reps, ncol = 3)
CIwidth <- numeric(reps)
CIcovered <- numeric(reps)
ALBIwidth <- numeric(1e3)
ALBIcovered <- numeric(1e3)
Sigma0 <- h2*K + (1-h2)*diag(1,n)
XSigmaInvX <- crossprod(X, solve(Sigma0, X))
h2s <- rep(0, 1e4)

for (jj in 1:reps) {
  
  # ---------------------------------------------------------------
  # Generate response
  # ---------------------------------------------------------------
  set.seed(jj + reps)
  y <- simData_K1(sigma2_g = h2*s2p, sigma2_e = (1-h2)*s2p, X, Ksqrt)
  ynew <- crossprod(V, y)
  partrue <- c(h2, s2p) # initialize at the truth
  opt <- optim(par=partrue, neg_loglik2, X=Xnew, y=ynew, lambda=lambda, method = "L-BFGS-B", lower = c(1e-4, 1e-4), upper = c(1-1e-4, Inf))#, hessian = FALSE) # , control = list(trace = 10)
  parhat <- opt$par
  h2s[jj] <- opt$par[1]

  # --------------------------------------------------------------------
  # LRT; need marginal maximizer tilde-beta, and sp2(h_0^2) at h_0^2
  # --------------------------------------------------------------------
  tildebeta <- solve(XSigmaInvX,crossprod(X, solve(Sigma0, y)))
  XBtmp <- X%*%tildebeta
  hat.sigmap2 <- (1/(n-p))*t(y - XBtmp)%*%solve(Sigma0, y - XBtmp)

  tmp <- 2*(neg_loglik2(par = c(h2, hat.sigmap2), X=Xnew, y=ynew, lambda=lambda) - opt$value)
  lrt[jj,1] <- 1*(tmp < qchisq(0.95, df = 1))
  lrt[jj,2] <- 1*(tmp < qchisq(0.9, df = 1))
  lrt[jj,3] <- 1*(tmp < qchisq(0.80, df = 1))

  tmp <- varRatioTest1d(h2, ynew, Xnew, lambda)
  score[jj,1] <- 1*(tmp < qchisq(0.95, df = 1))
  score[jj,2] <- 1*(tmp < qchisq(0.9, df = 1))
  score[jj,3] <- 1*(tmp < qchisq(0.80, df = 1))

  tmpCI <- confInv(ynew, Xnew, lambda, c(0, 1), tolerance = 1e-12, confLevel = 0.95) 
  CIwidth[jj] <- tmpCI[2] - tmpCI[1]
  CIcovered[jj] <- 1*(h2 >= tmpCI[1] & h2 <= tmpCI[2])
  cat(jj, "\n")
  
  result <- list("lrtCoverage" = lrt,
               "scoreCoverage" = score,
               "CIwidth" = CIwidth,
               "CIcovered" = CIcovered)
  
  saveRDS(result, file=savename)

}

