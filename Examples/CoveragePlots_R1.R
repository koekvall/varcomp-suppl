
tmp <- readRDS("/CoveragePlots/Coverage_n300_R1.RDS")
h2vec00 <- c("0.0001", "0.001", "0.01", "0.1", "0.25", "0.5", "0.75", "0.9", "0.99", "0.999", "0.9999")
h2vec <- c(10^seq(-4, -1, length=19), 0.25, 0.5, 0.75, rev(1-10^seq(-4, -1, length=19)))
h2vectmp <- c(1:19, 25, 31, 37, 43, 44:61)
waldCIs <- tmp$waldCIs
lrtCIs <- tmp$lrtCIs
scoreCIs <- tmp$scoreCIs
lb <- c(waldCIs[2,], lrtCIs[2,], scoreCIs[2,])
ub <- c(waldCIs[1,], lrtCIs[1,], scoreCIs[1,])
mp <- c(apply(waldCIs, 2, mean), apply(lrtCIs, 2, mean),apply(scoreCIs, 2, mean))
Method <- c(rep("Wald", length(h2vectmp)), rep("LRT", length(h2vectmp)), rep("Proposed", length(h2vectmp)))
h2 <- rep(h2vectmp, 3)
dat <- data.frame(h2, lb, ub, mp,  Method)
library(ggplot2)
p1 <- ggplot(dat, aes(x = h2, y = mp, color=Method)) + ggtitle("n = 300")+   xlab(expression(h^2)) + 
  ylab("Coverage probability") + coord_cartesian(ylim = c(0.845, 1)) + 
  geom_ribbon(aes(x = h2, ymin = lb, ymax = ub, fill = Method), alpha= 0.8, colour = NA) +
  #geom_line() +  # Add line plot
  theme_bw() +   # Black and white theme
  # geom_vline(xintercept = 8, col="red", linewidth = 0.25) + 
  # geom_vline(xintercept = 16, col="red", linewidth = 0.25) + 
  scale_fill_grey(start = 0.0, end = 0.85) +  # Greyscale for ribbon fill
  theme_bw() + geom_hline(yintercept = 0.95) + scale_x_continuous(breaks = seq(1, 61, by = 6), labels = h2vec00) +
  theme(plot.title = element_text(hjust = 0.5))


tmp <- readRDS("/CoveragePlots/Coverage_n1000_R1.RDS")

h2vec00 <- c("0.0001", "0.001", "0.01", "0.1", "0.25", "0.5", "0.75", "0.9", "0.99", "0.999", "0.9999")
h2vec <- c(10^seq(-4, -1, length=19), 0.25, 0.5, 0.75, rev(1-10^seq(-4, -1, length=19)))
h2vectmp <- c(1:19, 25, 31, 37, 43, 44:61)
waldCIs <- tmp$waldCIs
lrtCIs <- tmp$lrtCIs
scoreCIs <- tmp$scoreCIs
lb <- c(waldCIs[2,], lrtCIs[2,], scoreCIs[2,])
ub <- c(waldCIs[1,], lrtCIs[1,], scoreCIs[1,])
mp <- c(apply(waldCIs, 2, mean), apply(lrtCIs, 2, mean),apply(scoreCIs, 2, mean))
Method <- c(rep("Wald", length(h2vectmp)), rep("LRT", length(h2vectmp)), rep("Proposed", length(h2vectmp)))
h2 <- rep(h2vectmp, 3)
dat <- data.frame(h2, lb, ub, mp,  Method)
library(ggplot2)
p2 <- ggplot(dat, aes(x = h2, y = mp, color=Method)) + ggtitle("n = 1000")+   xlab(expression(h^2)) + 
  ylab("") + coord_cartesian(ylim = c(0.90, 1)) + 
  geom_ribbon(aes(x = h2, ymin = lb, ymax = ub, fill = Method), alpha= 0.8, colour = NA) +
  #geom_line() +  # Add line plot
  theme_bw() +   # Black and white theme
  # geom_vline(xintercept = 8, col="red", linewidth = 0.25) + 
  # geom_vline(xintercept = 16, col="red", linewidth = 0.25) + 
  scale_fill_grey(start = 0.0, end = 0.85) +  # Greyscale for ribbon fill
  theme_bw() + geom_hline(yintercept = 0.95) + scale_x_continuous(breaks = seq(1, 61, by = 6), labels = h2vec00) +
  theme(plot.title = element_text(hjust = 0.5))


# 
library(cowplot)
pdf("/CoveragePlots/CoveragePlot.pdf", width=11, height=3)
plot_grid(p1 + theme(legend.position="", axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)), p2+ theme(legend.position="", axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)), nrow=1, rel_widths=c(1, 1), align="h")
dev.off()

pdf("/CoveragePlots/CoveragePlot_Legend.pdf", width=5, height=5)
plot_grid(p1 + theme(legend.position="bottom"))
dev.off()







tmp <- readRDS("/CoveragePlots/Coverage_n300_Zoomed_R1.RDS")

h2vec <- log10(c(10^seq(-2.2, -1, length=21)))
h2vectmp <- 1:21
waldCIs <- tmp$waldCIs
lrtCIs <- tmp$lrtCIs
scoreCIs <- tmp$scoreCIs
lb <- c(waldCIs[2,], lrtCIs[2,], scoreCIs[2,])
ub <- c(waldCIs[1,], lrtCIs[1,], scoreCIs[1,])
mp <- c(apply(waldCIs, 2, mean), apply(lrtCIs, 2, mean),apply(scoreCIs, 2, mean))
Method <- c(rep("Wald", length(h2vec)), rep("LRT", length(h2vec)), rep("Proposed", length(h2vec)))
h2 <- rep(h2vectmp, 3)
dat <- data.frame(h2, lb, ub, mp, h2, Method)
p3 <- ggplot(dat, aes(x = h2, y = mp, color=Method)) + ggtitle("") +  
  xlab(expression(log[10](h^2))) + ylab("Coverage probability")  + 
  coord_cartesian(ylim = c(0.845, 1)) +    ggtitle("n = 300") + 
  geom_ribbon(aes(x = h2, ymin = lb, ymax = ub, fill = Method), alpha= 0.8,colour = NA)+
  #geom_line() +  # Add line plot
  theme_bw() +   # Black and white theme
  scale_fill_grey(start = 0.0, end = 0.85) +  # Greyscale for ribbon fill
  geom_hline(yintercept = 0.95, linetype = "solid") +
  scale_x_continuous(breaks = seq(1,21, by=4), labels = round(h2vec[seq(1,21, by=4)], 2)) +
  theme(plot.title = element_text(hjust = 0.5)) # Center the title



tmp <- readRDS("/CoveragePlots/Coverage_n1000_Zoomed_R1.RDS")


h2vec <- log10(c(10^seq(-2.2, -1, length=21)))
h2vectmp <- 1:21
waldCIs <- tmp$waldCIs
lrtCIs <- tmp$lrtCIs
scoreCIs <- tmp$scoreCIs
lb <- c(waldCIs[2,], lrtCIs[2,], scoreCIs[2,])
ub <- c(waldCIs[1,], lrtCIs[1,], scoreCIs[1,])
mp <- c(apply(waldCIs, 2, mean), apply(lrtCIs, 2, mean),apply(scoreCIs, 2, mean))
Method <- c(rep("Wald", length(h2vec)), rep("LRT", length(h2vec)), rep("Proposed", length(h2vec)))
h2 <- rep(h2vectmp, 3)
dat <- data.frame(h2, lb, ub, mp, h2, Method)
p4 <- ggplot(dat, aes(x = h2, y = mp, linetype = Method)) +
  ggtitle("n = 1000") +
  xlab(expression(log[10](h^2)))  +
 ylab("")  +
  coord_cartesian(ylim = c(0.845, 1)) +
  geom_ribbon(aes(x = h2, ymin = lb, ymax = ub, fill = Method), alpha= 0.8, colour = NA) +
  #geom_line() +  # Add line plot
  theme_bw() +   # Black and white theme
  scale_fill_grey(start = 0.0, end = 0.85) +  # Greyscale for ribbon fill
  geom_hline(yintercept = 0.95, linetype = "solid") +
  scale_x_continuous(breaks = seq(1,21, by=4), labels = round(h2vec[seq(1,21, by=4)], 2)) +
  theme(plot.title = element_text(hjust = 0.5)) # Center the title

library(cowplot)
pdf("/CoveragePlots/CoverageZoomed.pdf", width=11, height=3)
plot_grid(p3 + theme(legend.position=""), p4+ theme(legend.position=""), nrow=1, rel_widths=c(1, 1), align="h")
dev.off()


