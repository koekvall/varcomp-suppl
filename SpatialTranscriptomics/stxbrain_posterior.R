source("helper_functions.R")
library(Rcpp)
library(lmmvar)
library(Seurat)
library(SeuratData)
library(ggplot2)
library(patchwork)
library(dplyr)
library(MASS)
library(latex2exp)
library(matrixcalc)
library(SpatialExperiment)
library(STexampleData)
library(scran)
library(nnSVG)
library(cowplot)

# -----------------------------------
# InstallData("stxBrain")
# -----------------------------------
brainp.obj <- LoadData("stxBrain", type = "posterior1")    # 31053 features across 3353 samples

# -----------------------------------
# find boundary points
# -----------------------------------
indexp_bd <- findboundary(brainp.obj@images$posterior1@coordinates$row, brainp.obj@images$posterior1@coordinates$col)
length(indexp_bd)  # 342 spots

# -----------------------------------------------------
# filter genes with low expression and log-normalize
# -----------------------------------------------------
brainp <- GetAssayData(brainp.obj, slot = 'counts') # 31053  3353
# brain_norm <- as.matrix(LogNormalize(brain))       # normalized on each cell across features, than log1p
bp <- brainp[which(rowSums(brainp) >= 3 & rowMeans(brainp>0)*100>1), -indexp_bd]  # filter out low expression genes that have count less than 3 or appear in less than 1% of cells
bp <- as.matrix(LogNormalize(bp))                   # 14738  3011
coords <- matrix(c(brainp.obj@images$posterior1@coordinates$imagerow[-indexp_bd], brainp.obj@images$posterior1@coordinates$imagecol[-indexp_bd]), ncol = 2) # change to pixel coordinates

# -----------------------------------
# compute K and X as intercept
# -----------------------------------
l <- 0.02  
n <- nrow(coords)  # number of slots
p <- 1             # only intercepts as coveriate
X <- rep(1,n)
K <- spatialCor(coords = coords, l, scale = TRUE)
#save(K, file = "K_stxbrainp.RData")
#save(bp, file = "bp_stxbrainp.RData")
timep <- list("Eigendecomp" = 0, "Xtrans" = 0, "Ytrans" = rep(0, nrow(bp)), "CI" = rep(0, nrow(bp)), "total" = 0)

# -----------------------------------
# start
# -----------------------------------
starttime0 <- Sys.time()
eigen_decomp = eigen(K, symmetric = T)
timep$Eigendecomp = difftime(Sys.time(), starttime0, units = "secs") # 30.46s

V = eigen_decomp$vectors
lambda = eigen_decomp$values

# -----------------------------------
# calculating ci
# -----------------------------------
starttime = Sys.time()
Xnew = crossprod(V, X)
timep$Xtrans = difftime(Sys.time(), starttime, units = "secs")  # 0.015s

scoreCIp = matrix(0, nrow = nrow(bp), ncol = 6)
for (i in 1:nrow(bp)) {
  starttime = Sys.time()
  ynew = crossprod(V, bp[i,])
  timep$Ytrans[i] = difftime(Sys.time(), starttime, units = "secs")

  starttime = Sys.time()
  scoreCIp[i,5:6] = confInv(ynew, Xnew, lambda, confLevel = 0.95, tolerance = 1e-8)
  timep$CI[i] = difftime(Sys.time(), starttime, units = "secs")   # store running time for 1 CI

  if (i%%100 == 0) print(i)
}

# -----------------------------------
# total time used:
# -----------------------------------
timep$total = difftime(Sys.time(), starttime0, units = "secs")   #  1CI: 270.4814 s
sum(timep$Ytrans) # 209.04s
sum(timep$CI) # 28.69184s
#save(timep, file = "time_stxbrainp.Rdata")

scoreCIp = round(scoreCIp, digits = 8)
highp <- order(scoreCIp[,5], decreasing = T, na.last = NA) # rank genes according to 95% CI lower bound

# -----------------------------------
# confidence interval plot
# -----------------------------------
plotCIp <- data.frame(rank = 1:length(highp), l95 = scoreCIp[highp,5], u95 = scoreCIp[highp,6]) %>%
  ggplot(aes(x=rank, ymin = l95, ymax = u95)) +
  geom_errorbar(aes(ymin=l95, ymax=u95, width=.8)) +
  ggtitle("")  + theme_bw() +
  labs(y=TeX("$h^2$"))+
  theme(plot.title = element_text(hjust = 0.5))+
  theme(axis.title=element_text(size=14), plot.title = element_text(size = 16), axis.text = element_text(size=10))

# -----------------------------------
# histogram of 95% lower bound
# -----------------------------------
scorelowerp = rep(0,nrow(bp))
for (i in 1:nrow(bp)) {
  ynew = crossprod(V, bp[i,])
  scorelowerp[i] = confInv(ynew, Xnew, lambda, confLevel = 0.95, tolerance = 1e-8, type = "lower_bd")[1]
  if (i%%100 == 0) print(i)
}
#stxbrainpostresult <- list(scoreCIp, scorelowerp)
#save(stxbrainpostresult, file = "stxbrainpostresult.RData")
#scoreCIp <- stxbrainpostresult[[1]]

highhistp <- order(scorelowerp, decreasing = T, na.last = NA) # rank genes according to 95% CI lower bound
histp <- data.frame(lowerbound = scorelowerp[highhistp]) %>%
  ggplot(aes(x=lowerbound)) +
  geom_histogram(binwidth=0.02, alpha=1) +
  ggtitle("")  + theme_bw()+ xlab("lower bound") +
  theme(plot.title = element_text(hjust = 0.5))+
  theme(axis.title=element_text(size=14), plot.title = element_text(size = 16), axis.text = element_text(size=10))
rownames(bp)[highhistp[1:5]]                     # name of the 5 genes with highest 95% CI lower bound

# -----------------------------------
# spatial plots for top gene.
# -----------------------------------

## remove outside-of-tissue dots
r <- brainp.obj@images$posterior1@coordinates$row
c <- brainp.obj@images$posterior1@coordinates$col
plot(c[c>100 & c <120 & r<20],r[c>100 & c <120& r<20])
which(c>100 & c <120 & r<13) #indices: 1906, 2821

fea1p <- SpatialFeaturePlot(brainp.obj[, -c(1906, 2821)], features = rownames(bp)[highp[4]])

pdf("brainp_last.pdf", width=12, height=3.4) # 3.4
plot_grid(plotCIp , histp, fea1p, nrow=1, align="h", axis = "tblr")
dev.off()




