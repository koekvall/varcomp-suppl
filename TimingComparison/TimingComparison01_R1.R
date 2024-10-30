# ----------------------------------------------
# Run timing simulations
# ----------------------------------------------
# setwd("TimingComparison")
uu <- as.numeric(Sys.getenv('SLURM_ARRAY_TASK_ID'))
library(Rcpp)
library(lmmvar)
# run simulation
h2 <- 0.01   # c("0.0001", "0.001", "0.01", "0.1", "0.25", "0.5", "0.75", "0.9", "0.99", "0.999")
s2p <- 1
partrue <- c(h2, s2p) # initialize at the truth
nseq <- c(seq(200, 2000, length=9), c(2250, 2500, 2750, 3000))      # size of K
time <- list()
time$albi <- rep(0, length(nseq))
time$proposed <- rep(0, length(nseq))
time$fiesta <- rep(0, length(nseq))
time$rlrt <- rep(0, length(nseq))

ciseq <- seq(from=0, to=5, length.out=200)  # sequence of tested values for sigma^2g / sigma^e

neg_loglik2 <- function(par, X, y, eigens) {
  n = nrow(X)
  Sigma = (par[1]*eigens + rep(1-par[1],n)) * par[2]
  Sigma_inv = 1 / Sigma
  XSX = crossprod(X, Sigma_inv * X)
  XSX_inv = chol2inv(chol(XSX))
  betahat = XSX_inv %*% crossprod(X, Sigma_inv * y)
  l <- sum(log(Sigma)) / 2 +
    determinant(XSX)$modulus[1] / 2 +
    sum((y-X%*%betahat) * Sigma_inv * (y-X%*%betahat)) / 2
  return(l)
}




for (i in 1:length(nseq)) {

  dat <- readRDS(paste("data_R1/dat_", nseq[i], ".RDS", sep=""))
  X <- dat$X; dat$X <- NULL
  K <- dat$K; dat$K <- NULL
  V <- dat$V; dat$V <- NULL
  eigens <- dat$lambda; dat <- NULL
  Ksqrt <- V%*%diag(eigens^0.5)%*%t(V)

  set.seed(uu)
  y <- crossprod(Ksqrt, rnorm(nseq[i], sd = sqrt(h2*s2p))) + rnorm(nseq[i], sd = (1-h2)*s2p)

  # proposed score test method
  starttime <- Sys.time()
  Xnew <- crossprod(V, X)
  ynew <- crossprod(V, y)  # this step is O(n^2)
  ci <- confInv(ynew, Xnew, eigens)
  time$proposed[i] <- difftime(Sys.time(), starttime, units = "secs") # time of y transformation and confInv function

  # albi and fiesta methods
  # find mle first using optim
  opt <- optim(par=partrue, neg_loglik2, X=Xnew, y=ynew, eigens=eigens, method = "L-BFGS-B", lower = c(0, 1e-4), upper = c(1, Inf))#, hessian = FALSE) # , control = list(trace = 10)
  write.table(c(opt$par[1]), paste0("data_R1/ests_",nseq[i],"_", uu, ".txt"), sep="\n",col.names = FALSE, row.names = FALSE)  # requirement: one per line


  # command line run albi
  starttime <- Sys.time()
  sysoutput <- try(system(paste0("python3 ../albi_lib/albi.py -k data_R1/evalue_", nseq[i], ".txt -v data_R1/evector_", nseq[i],
                                 ".txt -x data_R1/covs_", nseq[i], ".txt -f data_R1/ests_", nseq[i],"_", uu, ".txt"), intern = T))      # executable file cannot upload, so use albi_lib/albi.py directly
  time$albi[i] <- difftime(Sys.time(), starttime, units = "secs")

  # command line run fiesta
  starttime <- Sys.time()
  sysoutput <- try(system(paste0("python3 ../albi_lib/fiesta.py -k data_R1/evalue_", nseq[i], ".txt -v data_R1/evector_", nseq[i],
                                 ".txt -x data_R1/covs_", nseq[i], ".txt -f data_R1/ests_", nseq[i],"_", uu, ".txt"), intern = T))
  time$fiesta[i] <- difftime(Sys.time(), starttime, units = "secs")


  # RLRT method
  starttime <- Sys.time()
  SimDists <- list()
  for (j in 1:length(ciseq)) {
    # SimDists[[j]] <- RLRTSim(Xnew, diag(eigens^0.5), lambda0 = ciseq[j]) # got sim-distribution for each tested value, grid=200
    SimDists[[j]] <- RLRTSim(X, Ksqrt, lambda0 = ciseq[j])
  }
  # time$RLRsim[i] = time$RLRsim[i] + difftime(Sys.time(), starttime, units = "secs")
  # starttime <- Sys.time()
  ci <- RLRTCI(SimDists=SimDists, ciseq=ciseq, Xnew=Xnew, ynew=ynew, eigens=eigens)
  time$rlrt[i] = difftime(Sys.time(), starttime, units = "secs")

  results <- list("time" = time)
  saveRDS(results, paste("~/Results_R1/h2_01_Replicate", uu, ".RDS", sep=""))
}



