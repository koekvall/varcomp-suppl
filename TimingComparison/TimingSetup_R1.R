# ---------------------------------------
# Data setup for ALBI comparisons
# ---------------------------------------
p <- 5
nseq <- c(seq(200, 2000, length=9), c(2250, 2500, 2750, 3000))
for (i in 1:length(nseq)) {
	set.seed(0)
	X <- matrix(rnorm(nseq[i]*p), nrow=nseq[i])
	K <- matrix(0, nrow = nseq[i], ncol = nseq[i])
	for (kk in 1:nseq[i]) {
		for (jj in 1:nseq[i]) {
	  		K[jj,kk] <- 0.95 ^ abs(jj - kk)
		}
	}
	eo <- eigen(K, symmetric = T)
	V <- eo$vectors
	lambda <- eo$values
	write.table(lambda, paste0("data_R1/evalue_",nseq[i],".txt"), sep="\n",col.names = FALSE, row.names = FALSE)  # requirement: one per line
	write.table(V, paste0("data_R1/evector_",nseq[i],".txt"), col.names = FALSE, row.names = FALSE)
	write.table(X, paste0("data_R1/covs_",nseq[i],".txt"), col.names = FALSE, row.names = FALSE)
	dat <- list("X" = X, "K" = K, "V" = V, "lambda" = lambda)
	saveRDS(dat, paste("data_R1/dat_", nseq[i], ".RDS", sep = ""))
}
