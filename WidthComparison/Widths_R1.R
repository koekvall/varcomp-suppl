# -------------------------------------------------
# set n
# -------------------------------------------------
uu <- as.numeric(Sys.getenv('SLURM_ARRAY_TASK_ID'))
params <-  expand.grid(n = c(200, 500, 1000, 2000), 
  h2 = c(0.0001, .001, .01, .1, 0.25, 0.5, 0.75, 0.9, 0.99, 0.999))

n <- params[uu,1]
h2 <- params[uu,2]
savename <- paste("~/blue/lmmvar/Width/Results_R1/Width_", uu, ".RDS", sep="")
# library(devtools)
# devtools::install_github("yqzhang5972/lmmvar")
library(MASS)
library(profvis)
library(matrixcalc)
library(lmmvar)


# negative log-likelihood function for h2, s2p
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
s2p <- 1

# ----------------------------------
# Model setup
# ---------------------------------
set.seed(0)
p <- 5
X <- matrix(rnorm(n*p), nrow=n)
K <- matrix(0, nrow=n, ncol=n)
for(j in 1:n){
  for(k in 1:n){
    K[j,k] <- 0.95 ^ abs(j-k)
  }
}
eoK <- eigen(K)
V <- eoK$vec
lambda <- eoK$val
Ksqrt <- V%*%diag(lambda^0.5)%*%t(V)
Xnew <- crossprod(V, X)
write.table(lambda, paste0("AuxillaryFiles1_R1/evalue_",uu,".txt"), sep="\n",col.names = FALSE, row.names = FALSE)  # requirement: one per line
write.table(V, paste0("AuxillaryFiles1_R1/evector_",uu,".txt"), col.names = FALSE, row.names = FALSE)
write.table(X, paste0("AuxillaryFiles1_R1/covs_",uu,".txt"), col.names = FALSE, row.names = FALSE)

reps <- 1e3
score <- matrix(0, nrow = reps, ncol = 3)
CIwidth <- numeric(reps)
CIcovered <- numeric(reps)
ALBIwidth <- numeric(reps)
ALBIcovered <- numeric(reps)
h2s <- rep(0, reps)


for (jj in 1:reps) {
  
  # ---------------------------------------------------------------
  # Generate response
  # ---------------------------------------------------------------
  set.seed(jj)
  y <- crossprod(Ksqrt, rnorm(n, sd = sqrt(h2*s2p))) + rnorm(n, sd = sqrt((1-h2)*s2p))
  ynew <- crossprod(V, y)
  partrue <- c(h2, s2p) # initialize at the truth
  opt <- optim(par=partrue, neg_loglik2, X=Xnew, y=ynew, lambda=lambda, method = "L-BFGS-B", lower = c(0, 1e-8), upper = c(1, Inf))#, hessian = FALSE) # , control = list(trace = 10)
  parhat <- opt$par
  h2s[jj] <- opt$par[1]

  tmp <- varRatioTest1d(h2, ynew, Xnew, lambda)
  score[jj,1] <- 1*(tmp < qchisq(0.95, df = 1))
  score[jj,2] <- 1*(tmp < qchisq(0.9, df = 1))
  score[jj,3] <- 1*(tmp < qchisq(0.80, df = 1))

  tmpCI <- confInv(ynew, Xnew, lambda) 
  CIwidth[jj] <- tmpCI[2] - tmpCI[1]
  CIcovered[jj] <- 1*(h2 >= tmpCI[1] & h2 <= tmpCI[2])

  result <- list("scoreCoverage" = score,
               "CIwidth" = CIwidth,
               "CIcovered" = CIcovered)

  saveRDS(result, file=savename)

}

  write.table(h2s, paste0("AuxillaryFiles1_R1/ests_", uu, ".txt"), sep="\n", col.names = FALSE, row.names = FALSE)  # requirement: one per line
  # command line run FIESTA
  sysoutput <- try(system(paste0("python3 albi_lib/fiesta.py -k AuxillaryFiles1_R1/evalue_", uu, ".txt -v AuxillaryFiles1_R1/evector_", uu,
                                  ".txt -x AuxillaryFiles1_R1/covs_", uu, ".txt -f AuxillaryFiles1_R1/ests_", uu, ".txt"), intern = T))
  #sysoutput <- try(system(paste0("python3 albi_lib/fiesta.py -k AuxillaryFiles1/evalue_", uu, ".txt -v AuxillaryFiles1/evector_", uu,
  #                                ".txt -f AuxillaryFiles1/ests_", uu, ".txt"), intern = T))
  t0 <- lapply(sysoutput, function(x){strsplit(x[length(x)], "\t")})
  for (kk in 1:reps) {
    ALBIwidth[kk] <- as.numeric(t0[[kk+2]][[1]][3]) -  as.numeric(t0[[kk+2]][[1]][2])
    ALBIcovered[kk] <- 1*(h2 <= as.numeric(t0[[kk+2]][[1]][3]) & h2 >= as.numeric(t0[[kk+2]][[1]][2]))
  }

  result <- list("lrtCoverage" = lrt, 
             "scoreCoverage" = score,
             "CIwidth" = CIwidth,
             "CIcovered" = CIcovered,
             "ALBIwidth" = ALBIwidth,
             "ALBIcovered" = ALBIcovered)

  saveRDS(result, file=savename)

