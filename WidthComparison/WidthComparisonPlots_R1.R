params <-  expand.grid(n = c(200, 500, 1000, 2000), 
                       h2 = c(0.0001, .001, .01, .1, 0.25, 0.5, 0.75, 0.9, 0.99, 0.999))
library(ggplot2)
library(reshape2)

uu <- 1
results <- matrix(0, 1000, 10)
for(kk in which(params[,1] == 200)){
  tmp <- readRDS(paste("/WidthResults/Width_", kk, ".RDS", sep=""))
  cat(uu, "\n")
  results[,uu] <- tmp$ALBIwidth/tmp$CIwidth
  uu <- uu + 1
}
colnames(results) <- c(0.0001, .001, .01, .1, 0.25, 0.5, 0.75, 0.9, 0.99, 0.999)

results.vec <- c(results)
h2.vec <- rep(c("0.0001", "0.001", "0.01", "0.1", "0.25", "0.5", "0.75", "0.9", "0.99", "0.999"), each=1000)
dat <- data.frame(results.vec, h2.vec)
# Create side-by-side violin plots
p1 <- ggplot(dat, aes(x=h2.vec, y=results.vec)) + 
  geom_violin(trim = FALSE) +  # Create violin plot
  stat_summary(fun = mean, geom = "point", shape = 20, size = 3, color = "black", position = position_dodge(width = 0.9)) +  # Mean
  stat_summary(fun = median, geom = "point", shape = 18, size = 3, color = "grey", position = position_dodge(width = 0.9)) +  # Median
  coord_cartesian(ylim = c(0, 3)) + ggtitle("n = 200") + xlab(expression(h^2)) + ylab("ALBI to Score width ratio") + 
  theme_bw() +   theme(plot.title = element_text(hjust = 0.5)) +   geom_hline(yintercept = 1, color = "red", linetype = "solid", size = 0.3) 



uu <- 1
results <- matrix(0, 1000, 10)
for(kk in which(params[,1] == 500)){
  tmp <- readRDS(paste("/WidthResults/Width_", kk, ".RDS", sep=""))
  cat(uu, "\n")
  results[,uu] <- tmp$ALBIwidth/tmp$CIwidth
  uu <- uu + 1
}
colnames(results) <- c(0.0001, .001, .01, .1, 0.25, 0.5, 0.75, 0.9, 0.99, 0.999)

results.vec <- c(results)
h2.vec <- rep(c("0.0001", "0.001", "0.01", "0.1", "0.25", "0.5", "0.75", "0.9", "0.99", "0.999"), each=1000)
dat <- data.frame(results.vec, h2.vec)
# Create side-by-side violin plots
p2 <- ggplot(dat, aes(x=h2.vec, y=results.vec)) + 
  geom_violin(trim = FALSE) +  # Create violin plot
  stat_summary(fun = mean, geom = "point", shape = 20, size = 3, color = "black", position = position_dodge(width = 0.9)) +  # Mean
  stat_summary(fun = median, geom = "point", shape = 18, size = 3, color = "grey", position = position_dodge(width = 0.9)) +  # Median
  coord_cartesian(ylim = c(0, 3)) + ggtitle("n = 500") + xlab(expression(h^2)) + ylab("") + 
  theme_bw() +   theme(plot.title = element_text(hjust = 0.5)) +   geom_hline(yintercept = 1, color = "red", linetype = "solid", size = 0.3) 





uu <- 1
results <- matrix(0, 1000, 10)
for(kk in which(params[,1] == 1000)){
  tmp <- readRDS(paste("/WidthResults/Width_", kk, ".RDS", sep=""))
  cat(uu, "\n")
  results[,uu] <- tmp$ALBIwidth/tmp$CIwidth
  uu <- uu + 1
}
colnames(results) <- c(0.0001, .001, .01, .1, 0.25, 0.5, 0.75, 0.9, 0.99, 0.999)

results.vec <- c(results)
h2.vec <- rep(c("0.0001", "0.001", "0.01", "0.1", "0.25", "0.5", "0.75", "0.9", "0.99", "0.999"), each=1000)
dat <- data.frame(results.vec, h2.vec)
# Create side-by-side violin plots
p3 <- ggplot(dat, aes(x=h2.vec, y=results.vec)) + 
  geom_violin(trim = FALSE) +  # Create violin plot
  stat_summary(fun = mean, geom = "point", shape = 20, size = 3, color = "black", position = position_dodge(width = 0.9)) +  # Mean
  stat_summary(fun = median, geom = "point", shape = 18, size = 3, color = "grey", position = position_dodge(width = 0.9)) +  # Median
  coord_cartesian(ylim = c(0, 3)) + ggtitle("n = 1000") + xlab(expression(h^2)) + ylab("ALBI to Score width ratio") + 
  theme_bw() +   theme(plot.title = element_text(hjust = 0.5)) +   geom_hline(yintercept = 1, color = "red", linetype = "solid", size = 0.3) 


uu <- 1
results <- matrix(0, 1000, 10)
for(kk in which(params[,1] == 2000)){
  tmp <- readRDS(paste("/WidthResults/Width_", kk, ".RDS", sep=""))
  cat(uu, "\n")
  results[,uu] <- tmp$ALBIwidth/tmp$CIwidth
  uu <- uu + 1
}
colnames(results) <- c(0.0001, .001, .01, .1, 0.25, 0.5, 0.75, 0.9, 0.99, 0.999)

results.vec <- c(results)
h2.vec <- rep(c("0.0001", "0.001", "0.01", "0.1", "0.25", "0.5", "0.75", "0.9", "0.99", "0.999"), each=1000)
dat <- data.frame(results.vec, h2.vec)
# Create side-by-side violin plots
p4 <- ggplot(dat, aes(x=h2.vec, y=results.vec)) + 
  geom_violin(trim = FALSE) +  # Create violin plot
  stat_summary(fun = mean, geom = "point", shape = 20, size = 3, color = "black", position = position_dodge(width = 0.9)) +  # Mean
  stat_summary(fun = median, geom = "point", shape = 18, size = 3, color = "grey", position = position_dodge(width = 0.9)) +  # Median
  coord_cartesian(ylim = c(0, 3)) + ggtitle("n = 2000") + xlab(expression(h^2)) + ylab("") + 
  theme_bw() +   theme(plot.title = element_text(hjust = 0.5)) +   geom_hline(yintercept = 1, color = "red", linetype = "solid", size = 0.3) 



library(cowplot)
pdf("/WidthComparison.pdf",  width=9, height=6.2)
plot_grid(p1+ theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)), p2+ theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)),
          p3+ theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)), p4 + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)), nrow=2)
dev.off()