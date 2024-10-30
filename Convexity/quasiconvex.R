library(lmmvar)
library(ggplot2)
library(MASS)
library(latex2exp)
library(cowplot)

# -------------------
# simulate y
# -------------------
simData <- function(sigma2_g, sigma2_e = 1, X, Ksqrt) {
  n <- dim(X)[1]
  p <- dim(X)[2]
  beta <- rep(0, p)
  e <- rnorm(n, mean = 0, sd = sqrt(sigma2_e))
  u <- rnorm(n, mean = 0, sd = sqrt(sigma2_g))
  y <- X %*% beta + Ksqrt %*% u + e
  return(y)
}


set.seed(1)
n <- 1000
p <- 5
lout <- 200
t1 <- rep(0, lout)
t2 <- rep(0, lout)
t3 <- rep(0, lout)
h2seq <- seq(0,1,length.out = lout)

K <- matrix(0, nrow = n, ncol = n)
for (kk in 1:n) {
  for (jj in 1:n) {
    K[jj,kk] <- 0.95 ^ abs(jj - kk)
  }
}
eo <- eigen(K)
Ksqrt <- eo$vec %*% diag(eo$val^0.5) %*% t(eo$vec)

### repeat simulation from here
eigen_decomp <- eigen(K, symmetric = T)
V <- eigen_decomp$vectors
lambda <- eigen_decomp$values
X <- matrix(rnorm(n*p, 0, 1), nrow = n)
Xnew <- crossprod(V, X)

h2true <- 0.0001
y <- simData(sigma2_g = h2true/(1-h2true), sigma2_e = 1, X, Ksqrt) # h2=0
ynew <- crossprod(V, y)
h2seq1 <- seq(0,0.1,length.out = lout)
for (k in 1:length(h2seq1)) {
  t1[k] <- varRatioTest1d(h2 = h2seq1[k], y=ynew, X=Xnew, lambda=lambda)
}

h2true <- 0.5
y <- simData(sigma2_g = h2true/(1-h2true), sigma2_e = 1, X, Ksqrt) # h2=0
ynew <- crossprod(V, y)
h2seq2 <- seq(0.35, 0.7,length.out = lout)
for (k in 1:length(h2seq2)) {
  t2[k] <- varRatioTest1d(h2 = h2seq2[k], y=ynew, X=Xnew, lambda=lambda)
}

h2true <- 0.9999 # cannot be 1
y <- simData(sigma2_g = h2true/(1-h2true), sigma2_e = 1, X, Ksqrt) # h2=0
ynew <- crossprod(V, y)
h2seq3 <- seq(0.98, 1,length.out = lout)

for (k in 1:length(h2seq3)) {
  t3[k] <- varRatioTest1d(h2 = h2seq3[k], y=ynew, X=Xnew, lambda=lambda)
}

plot1 <- data.frame(h2 = h2seq1, teststatistic = t1) %>%
  ggplot(aes(x = h2, y = teststatistic))  + theme_bw() + #geom_point(size=0.8) +
  geom_line() + xlab(expression(h^2)) + ylab(TeX("$T^R_n(h^2)$")) +
  ggtitle(TeX("$h^2_*=0.0001$")) + theme(plot.title = element_text(hjust = 0.5)) +
  geom_hline(yintercept=qchisq(0.95, 1), color = "red")+
  theme(axis.title=element_text(size=16), plot.title = element_text(size = 18), axis.text = element_text(size=12))

plot2 <- data.frame(h2 = h2seq2, teststatistic = t2) %>%
  ggplot(aes(x = h2, y = teststatistic))  + theme_bw() + #geom_point(size=0.8) +
  geom_line() + xlab(expression(h^2)) + ylab("") +
  ggtitle(TeX("$h^2_*=0.5$")) + theme(plot.title = element_text(hjust = 0.5)) +
  geom_hline(yintercept=qchisq(0.95, 1), color = "red")+
  theme(axis.title=element_text(size=16), plot.title = element_text(size = 18), axis.text = element_text(size=12))

plot3 <- data.frame(h2 = h2seq3, teststatistic = t3) %>%
  ggplot(aes(x = h2, y = teststatistic)) + theme_bw() + #geom_point(size=0.8) +
  geom_line() + xlab(expression(h^2)) + ylab("") +
  ggtitle(TeX("$h^2_*=0.9999$")) + theme(plot.title = element_text(hjust = 0.5)) +
  geom_hline(yintercept=qchisq(0.95, 1), color = "red")+
  theme(axis.title=element_text(size=16), plot.title = element_text(size = 18), axis.text = element_text(size=12))

pdf("ConvexityPlot.pdf", width=12, height=3.4)
plot_grid(plot1 + theme(legend.position=""), plot2, plot3, nrow=1, align="h")
dev.off()


