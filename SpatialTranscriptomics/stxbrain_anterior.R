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
# Install dataset from package
# InstallData("stxBrain")
# -----------------------------------
brain.obj <- LoadData("stxBrain", type = "anterior1")    # 31053 features across 2696 samples

# -----------------------------------
# find boundary points
# -----------------------------------
index_bd <- findboundary(brain.obj@images$anterior1@coordinates$row, brain.obj@images$anterior1@coordinates$col)
length(index_bd)  # 316 spots

# -----------------------------------
# filter genes with low expression and log-normalize
# -----------------------------------
brain <- GetAssayData(brain.obj, slot = 'counts') # 31053  2696
b <- brain[which(rowSums(brain) >= 3 & rowMeans(brain>0)*100>1), -index_bd]  # filter out low expression genes that have count less than 3 or appear in less than 1% of cells
b <- as.matrix(LogNormalize(b))                   # 15117 2380
coords <- matrix(c(brain.obj@images$anterior1@coordinates$imagerow[-index_bd], brain.obj@images$anterior1@coordinates$imagecol[-index_bd]), ncol = 2) # change to pixel coordinates

# -----------------------------------
# compute K and X as intercept
# -----------------------------------
l = 0.02    # based on result from NNSVG paper
n <- nrow(coords)  # number of slots
p <- 1             # only intercepts as coveriate
X <- rep(1,n)
K <- spatialCor(coords = coords, l, scale = TRUE)
#save(K, file = "K_stxbrain.RData")
#save(b, file = "b_stxbrain.RData")
time <- list("Eigendecomp" = 0, "Xtrans" = 0, "Ytrans" = rep(0, nrow(b)), "CI" = rep(0, nrow(b)), "total" = 0)

# -----------------------------------
# start
# -----------------------------------
starttime0 <- Sys.time()
eigen_decomp = eigen(K, symmetric = T)
time$Eigendecomp = difftime(Sys.time(), starttime0, units = "secs")

V = eigen_decomp$vectors
lambda = eigen_decomp$values

# -----------------------------------
# calculating ci
# -----------------------------------
starttime = Sys.time()
Xnew = crossprod(V, X)
time$Xtrans = difftime(Sys.time(), starttime, units = "secs")

scoreCI = matrix(0, nrow = nrow(b), ncol = 6)
for (i in 1:nrow(b)) {
  starttime = Sys.time()
  ynew = crossprod(V, b[i,])
  time$Ytrans[i] = difftime(Sys.time(), starttime, units = "secs")

  starttime = Sys.time()
  scoreCI[i,5:6] = confInv(ynew, Xnew, lambda, confLevel = 0.95, tolerance = 1e-8)
  time$CI[i] = difftime(Sys.time(), starttime, units = "secs")   # store running time for 1 CI
  if (i%%100 == 0) print(i)
}

# -----------------------------------
# total time used:
# -----------------------------------
time$total = difftime(Sys.time(), starttime0, units = "secs") # 3CI: 185.71s, 1CI: 160.62s, 3CI one-sided: 169.36s
# Time difference of 173.015 secs
sum(time$Ytrans) # 134.005s
sum(time$CI)     # 23.12813s
# save(time, file = "time_stxbrain.Rdata")

scoreCI = round(scoreCI, digits = 8)
high <- order(scoreCI[,5], decreasing = T, na.last = NA) # rank genes according to 95% CI lower bound

# -----------------------------------
# confidence interval plot
# -----------------------------------
plotCI <- data.frame(rank = 1:length(high), l95 = scoreCI[high,5], u95 = scoreCI[high,6]) %>%
  ggplot(aes(x=rank, ymin = l95, ymax = u95)) +
  geom_errorbar(aes(ymin=l95, ymax=u95, width=.8)) +
  ggtitle("two-sided confidence intervals")  + theme_bw() +
  labs(y=TeX("$h^2$")) +
  theme(plot.title = element_text(hjust = 0.5))+
  theme(axis.title=element_text(size=14), plot.title = element_text(size = 16), axis.text = element_text(size=10))

# -----------------------------------
# histogram of 95% lower bound
# -----------------------------------
scorelower = rep(0,nrow(b))
for (i in 1:nrow(b)) {
  ynew = crossprod(V, b[i,])
  scorelower[i] = confInv(ynew, Xnew, lambda, confLevel = 0.95, tolerance = 1e-8, type = "lower_bd")[1]
  if (i%%100 == 0) print(i)
}
#stxbrainresult <- list(scoreCI, scorelower)
#save(stxbrainresult, file = "stxbrainresult.RData")

highhist <- order(scorelower, decreasing = T, na.last = NA) # rank genes according to 95% CI lower bound
hist <- data.frame(lowerbound = scorelower[highhist]) %>%
  ggplot(aes(x=lowerbound)) +
  geom_histogram(binwidth=0.02, alpha=1) +
  ggtitle("one-sided lower bounds")  + theme_bw() + xlab("lower bound") +
  theme(plot.title = element_text(hjust = 0.5))+
  theme(axis.title=element_text(size=14), plot.title = element_text(size = 16), axis.text = element_text(size=10))
rownames(b)[highhist[1:5]]                     # name of the 5 genes with highest 95% CI lower bound

# -----------------------------------
# spatial plots for top three gene
# -----------------------------------
fea1 <- SpatialFeaturePlot(brain.obj, features = rownames(b)[highhist[3]])
#ggsave("1stfeature_brain.png", plot = fea1, dpi = "print", width = 10, height = 10, units = "cm")

pdf("brain_last.pdf", width=12, height=3.4) # 3.4
plot_grid(plotCI, hist, fea1, nrow=1, align="h", axis = "tblr")
dev.off()


