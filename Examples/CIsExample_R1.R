# -------------------------------------------------
# set n
# -------------------------------------------------
uu <- as.numeric(Sys.getenv('SLURM_ARRAY_TASK_ID'))
n <-  c(300, 1000)
n <- n[uu]
set.seed(1)
savename <- paste("~/blue/lmmvar/Examples/Coverage_n", n, "_R1.RDS", sep="")
library(MASS)
library(profvis)
library(matrixcalc)
library(lmmvar)

# -------------------------------------------------
# negative log-likelihood function for h2, s2p
# -------------------------------------------------
neg_loglik2 <- function(par, X, y, lambda) {
  n <- nrow(X)
  Sigma <- (par[1]*lambda + rep(1-par[1],n)) * par[2]
  Sigma_inv <- 1 / Sigma
  XSX <- crossprod(X, Sigma_inv * X)
  betahat <- solve(XSX, crossprod(X, Sigma_inv * y))
  l <- sum(log(Sigma)) / 2 +
    determinant(XSX)$modulus[1] / 2 +
    sum((y-X%*%betahat) * Sigma_inv * (y-X%*%betahat)) / 2
  return(l)
}



# -----------------------------------
# true parameter values; beta = 0 
# ----------------------------------
h2vec <- c(10^seq(-4, -1, length=19), 0.25, 0.5, 0.75, rev(1-10^seq(-4, -1, length=19)))
s2p <- 1

waldCIs <- matrix(0, nrow=2, ncol=length(h2vec))
lrtCIs <- matrix(0, nrow=2, ncol=length(h2vec))
scoreCIs <- matrix(0, nrow=2, ncol=length(h2vec))

# ---------------------------------------
# Make K
# ---------------------------------------
K <- matrix(0, nrow=n, ncol=n)
for(jj in 1:n){
  for(kk in 1:n){
    K[jj,kk] <- 0.95^abs(jj - kk)
  }
}
eo <- eigen(K)
V <- eo$vec
lambda <- eo$values
Ksqrt <- eo$vec%*%diag(eo$val^0.5)%*%t(eo$vec)


for(h2 in h2vec){
  
  # ----------------------------------
  # Model setup
  # ---------------------------------
  p <- 5
  set.seed(1)
  X <- matrix(rnorm(n*p), nrow=n, ncol=p)
  Xnew <- crossprod(V, X)
  # ---------------------------------------------------------------
  # Generate response
  # ---------------------------------------------------------------
  reps <- 5e3
  wald <- numeric(reps)
  lrt <- numeric(reps)
  score <- numeric(reps)
  
  Sigma0 <- h2*K + (1-h2)*diag(1,n)
  XSigmaInvX <- crossprod(X, solve(Sigma0, X))
  
  for (jj in 1:reps) {
    
    # --------------------------
    set.seed(jj)
    y <- crossprod(Ksqrt, rnorm(n, sd = sqrt(h2*s2p))) + rnorm(n, sd = sqrt((1-h2)*s2p))
    ynew <- crossprod(V, y)
    partrue <- c(h2, s2p) # initialize at the truth
    opt <- optim(par=partrue, neg_loglik2, X=Xnew, y=ynew, lambda=lambda, method = "L-BFGS-B", lower = c(1e-4, 1e-4), upper = c(1-1e-4, Inf))#, hessian = FALSE) # , control = list(trace = 10)
    parhat <- opt$par
    # Ihat <- opt$hessian # fisherhat(parhat, X, y, V, lambda, K)
    
    # Observed information at MLE
    S <- parhat[2]*(parhat[1]*K + (1 - parhat[1])*diag(1,n))
    eo <- eigen(S)
    SinvSqrt <- tcrossprod(eo$vec*matrix(eo$val^(-0.5), nrow= dim(eo$vec)[1], ncol=dim(eo$vec)[2], byrow=T), eo$vec)
    # ------------------------------
    # Using I_{11} from Prop 1
    # -------------------------------
    t0 <- crossprod(X, SinvSqrt)
    P <- crossprod(t0, solve(crossprod(X, solve(S, X)), t0))
    Q <- diag(1,n) - P
    D <- diag(1/(parhat[1]*lambda + (1 - parhat[1])))*diag(lambda-1,n)
    H <- tcrossprod(V*matrix(diag(D), nrow=dim(V)[1], ncol=dim(V)[2], byrow=T),V)
    tmp <- crossprod(t(Q), H)
    I11 <- 0.5*sum(diag(crossprod(t(tmp), tmp)))
    I12 <- 1/(2*parhat[2])*sum(diag(tmp))
    I22 <- (n-p)/(2*parhat[2]^2)
    I <- I22/(I11*I22 - I12^2)
    
    # --------------------------------------------------
    # wald test (hat-h^2 - h^_0) ~ N(0, I_{11}(...))
    # --------------------------------------------------
    wald[jj] <- 1*(2*pnorm(-abs(parhat[1]-h2)/sqrt(I)) < 0.05)
    
    # --------------------------------------------------
    # LRT; need marginal maximizer tilde-beta, and sp2(h_0^2)
    # at h_0^2
    # ---------------------------------------------------
    tildebeta <- solve(XSigmaInvX,crossprod(X, solve(Sigma0, y)))
    tmp <- X%*%tildebeta
    hat.sigmap2 <- (1/(n-p))*crossprod(y - tmp, solve(Sigma0, y - tmp))
    lrt[jj] <- 1*(2*(neg_loglik2(par = c(h2, hat.sigmap2), X=Xnew, y=ynew, lambda=lambda) - opt$value) > qchisq(0.95, df = 1))
    if(jj%%100==0){
      cat(jj, "\n")
    }
    score[jj] <- 1*(varRatioTest1d(h2, ynew, Xnew, lambda) > qchisq(0.95, 1))
  }



  cat(paste("Through ", which(h2vec == h2), " of ", length(h2vec)), "\n")
  waldCIs[,which(h2vec == h2)] <- 1 - prop.test(x = sum(wald), n = length(wald), correct=FALSE)$conf.int[1:2]
  lrtCIs[,which(h2vec == h2)] <- 1 - prop.test(x = sum(lrt), n = length(lrt), correct=FALSE)$conf.int[1:2]
  scoreCIs[,which(h2vec==h2)] <- 1 - prop.test(x = sum(score, na.rm=TRUE), n = length(score), correct =FALSE)$conf.int[1:2]
  
}

result <- list("waldCIs" = waldCIs, 
               "lrtCIs" = lrtCIs, 
               "scoreCIs" = scoreCIs)

saveRDS(result, file=savename)