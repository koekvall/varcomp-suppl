# setwd("TimingComparison")
temp <- readRDS("TimingResults.RDS")
nseq <- c(seq(200, 2000, length=9), c(2250, 2500, 2750, 3000)) #

for (i in 1:12) {
  temp[[i]] <- temp[[i]][,1:9] # Timing from 200 to 2000
}


# 01
mp <- c(apply(temp$ours.times01, 2, mean), apply(temp$fiesta.times01, 2, mean), apply(temp$albi.times01, 2, mean), apply(temp$rlrt.times01, 2, mean))
lb <- c(apply(temp$ours.times01, 2, mean) - 2*apply(temp$ours.times01, 2, sd)/sqrt(100),
        apply(temp$fiesta.times01, 2, mean)- 2*apply(temp$fiesta.times01, 2, sd)/sqrt(100),
        apply(temp$albi.times01, 2, mean) - 2*apply(temp$albi.times01, 2, sd)/sqrt(100),
        apply(temp$rlrt.times01, 2, mean) - 2*apply(temp$rlrt.times01, 2, sd)/sqrt(100))
ub <- c(apply(temp$ours.times01, 2, mean) + 2*apply(temp$ours.times01, 2, sd)/sqrt(100),
        apply(temp$fiesta.times01, 2, mean) + 2*apply(temp$fiesta.times01, 2, sd)/sqrt(100),
        apply(temp$albi.times01, 2, mean) + 2*apply(temp$albi.times01, 2, sd)/sqrt(100),
        apply(temp$rlrt.times01, 2, mean) + 2*apply(temp$rlrt.times01, 2, sd)/sqrt(100))
lb <- log(lb,base = 10)
ub <- log(ub,base = 10)
mp <- log(mp,base = 10)
Method <- c(rep("Proposed", length(nseq)), rep("FIESTA", length(nseq)), rep("ALBI", length(nseq)), rep("RLRT", length(nseq)))
Nseq <- rep(nseq, 4)
dat <- data.frame(Nseq, lb, ub, mp, Method)



library(ggplot2)
p1 <- ggplot(dat, aes(x = Nseq, y = mp, linetype = Method)) + ggtitle(expression(paste(h^2, " = " , 0.01))) +   xlab(expression(n)) +
  ylab(expression(paste("Computing time in ", log[10], "(secs)"))) + geom_line(linewidth=0.5)+
  geom_ribbon(aes(x = Nseq, ymin = lb, ymax = ub), alpha = 0.2) +
  scale_linetype_manual(values = c(3,4,1,2)) +
  theme_bw() + theme(legend.position="") +  ylim(c(-4, 4)) +
  theme(plot.title = element_text(hjust = 0.5))  #+ theme(legend.position="bottom")+
  # theme(axis.title=element_text(size=16), plot.title = element_text(size = 18), axis.text = element_text(size=12))


# 1
# nseq <- c(seq(200, 2000, length=9))
mp <- c(apply(temp$ours.times1, 2, mean), apply(temp$fiesta.times1, 2, mean), apply(temp$albi.times1, 2, mean), apply(temp$rlrt.times1, 2, mean))
lb <- c(apply(temp$ours.times1, 2, mean) - 2*apply(temp$ours.times1, 2, sd)/sqrt(100),
        apply(temp$fiesta.times1, 2, mean)- 2*apply(temp$fiesta.times1, 2, sd)/sqrt(100),
        apply(temp$albi.times1, 2, mean) - 2*apply(temp$albi.times1, 2, sd)/sqrt(100),
        apply(temp$rlrt.times1, 2, mean) - 2*apply(temp$rlrt.times1, 2, sd)/sqrt(100))
ub <- c(apply(temp$ours.times1, 2, mean) + 2*apply(temp$ours.times1, 2, sd)/sqrt(100),
        apply(temp$fiesta.times1, 2, mean) + 2*apply(temp$fiesta.times1, 2, sd)/sqrt(100),
        apply(temp$albi.times1, 2, mean) + 2*apply(temp$albi.times1, 2, sd)/sqrt(100),
        apply(temp$rlrt.times1, 2, mean) + 2*apply(temp$rlrt.times1, 2, sd)/sqrt(100))
lb <- log(lb,base = 10)
ub <- log(ub,base = 10)
mp <- log(mp,base = 10)
Method <- c(rep("Proposed", length(nseq)), rep("FIESTA", length(nseq)), rep("ALBI", length(nseq)), rep("RLRT", length(nseq)))
Nseq <- rep(nseq, 4)
dat <- data.frame(Nseq, lb, ub, mp, Method)
#Method <- factor(Method, ordered=TRUE, levels=c("ALBI", "Proposed", "FIESTA"))
library(ggplot2)
p2 <- ggplot(dat, aes(x = Nseq, y = mp, linetype = Method)) + ggtitle(expression(paste(h^2, " = " , 0.1))) +   xlab(expression(n)) +
  ylab("") + geom_line(linewidth=0.5)+
  geom_ribbon(aes(x = Nseq, ymin = lb, ymax = ub), alpha = 0.2) +
  scale_linetype_manual(values = c(3,4,1,2)) +
  theme_bw() + theme(legend.position="") +  ylim(c(-4, 4)) +
  theme(plot.title = element_text(hjust = 0.5))

  # ggplot(dat, aes(x = Nseq, y = mp, linetype = Method)) + ggtitle(expression(paste(h^2, " = " , 0.1))) +   xlab(expression(n)) +
  # ylab("") + geom_line(linewidth=0.5)+ geom_ribbon(aes(x = Nseq, ymin = lb, ymax = ub, fill = Method), alpha= 0.3, colour = NA) +
  # theme_bw() + ylim(c(-4, 4)) + theme(legend.position="") +
  # theme(plot.title = element_text(hjust = 0.5))




# 5
# nseq <- c(seq(200, 2000, length=9))
mp <- c(apply(temp$ours.times5, 2, mean), apply(temp$fiesta.times5, 2, mean), apply(temp$albi.times5, 2, mean), apply(temp$rlrt.times5, 2, mean))
lb <- c(apply(temp$ours.times5, 2, mean) - 2*apply(temp$ours.times5, 2, sd)/sqrt(100),
        apply(temp$fiesta.times5, 2, mean)- 2*apply(temp$fiesta.times5, 2, sd)/sqrt(100),
        apply(temp$albi.times5, 2, mean) - 2*apply(temp$albi.times5, 2, sd)/sqrt(100),
        apply(temp$rlrt.times5, 2, mean) - 2*apply(temp$rlrt.times5, 2, sd)/sqrt(100))
ub <- c(apply(temp$ours.times5, 2, mean) + 2*apply(temp$ours.times5, 2, sd)/sqrt(100),
        apply(temp$fiesta.times5, 2, mean) + 2*apply(temp$fiesta.times5, 2, sd)/sqrt(100),
        apply(temp$albi.times5, 2, mean) + 2*apply(temp$albi.times5, 2, sd)/sqrt(100),
        apply(temp$rlrt.times5, 2, mean) + 2*apply(temp$rlrt.times5, 2, sd)/sqrt(100))
lb <- log(lb,base = 10)
ub <- log(ub,base = 10)
mp <- log(mp,base = 10)
Method <- c(rep("Proposed", length(nseq)), rep("FIESTA", length(nseq)), rep("ALBI", length(nseq)), rep("RLRT", length(nseq)))
#Method <- factor(Method, ordered=TRUE, levels=c("ALBI", "Proposed", "FIESTA"))
Nseq <- rep(nseq, 4)
dat <- data.frame(Nseq, lb, ub, mp, Method)

library(ggplot2)
p3 <- ggplot(dat, aes(x = Nseq, y = mp, linetype = Method)) + ggtitle(expression(paste(h^2, " = " , 0.5))) +   xlab(expression(n)) +
  ylab("") + geom_line(linewidth=0.5)+
  geom_ribbon(aes(x = Nseq, ymin = lb, ymax = ub), alpha = 0.2) +
  scale_linetype_manual(values = c(3,4,1,2)) +
  theme_bw() + theme(legend.position="") +  ylim(c(-4, 4)) +
  theme(plot.title = element_text(hjust = 0.5))
# ggplot(dat, aes(x = Nseq, y = mp, color=Method)) + ggtitle(expression(paste(h^2, " = " , 0.5))) +   xlab(expression(n)) +
#   ylab("") + geom_line(linewidth=0.5)+ geom_ribbon(aes(x = Nseq, ymin = lb, ymax = ub, fill = Method), alpha= 0.3, colour = NA) +
#   theme_bw() + ylim(c(-4, 4)) +
#   theme(plot.title = element_text(hjust = 0.5))


library(cowplot)
pdf("time_3000.pdf", width=9, height=2.5)
plot_grid(p1, p2, p3 + theme(legend.position=""), nrow = 1, rel_widths = c(1, 1, 1))
dev.off()
pdf("time_legend3000.pdf", height=0.3)
# plot_grid(p3 + theme(legend.position="bottom"))
# plot_grid(get_legend(p3 + theme(legend.position = 'bottom') +
#                        theme(legend.title = element_text(size=15))+
#                        theme(legend.text = element_text(size=15))))
leg <- cowplot::get_plot_component(p3 + theme(legend.position = 'bottom') +
                                     theme(legend.key.width = unit(1, 'cm')) +
                                     theme(legend.title = element_text(size=15)) +
                                     theme(legend.text = element_text(size=15)), 'guide-box-bottom', return_all = TRUE)
cowplot::ggdraw(leg)
dev.off()






# -----------------------------------
# Add width density plots
# -----------------------------------
params <-  expand.grid(n = c(500, 1000),
  h2 = c(0.0001, .001, .01, .1, 0.25, 0.5, 0.75, 0.9, 0.99, 0.999))
library(ggplot2)
path <- "/Users/amolstad/Desktop/TimingPlots/Results/"

t0 <- readRDS(paste(path, "Width_2.RDS", sep=""))
data <- data.frame(widths = c(t0$CIwidth, t0$ALBIwidth),
        CI = rep(c("Score", "FIESTA"), each = length(t0$CIwidth)))
s1 <- ggplot(data, aes(x=widths, fill=CI)) + geom_density(alpha=0.25) + ggtitle(expression(paste(h^2, " = 0.0001, n = 1000")))+
 theme_bw() + xlab("CI width") + ylab("Density")+ theme(legend.position="") +
  theme(plot.title = element_text(hjust = 0.5))
  mean(t0$ALBIcovered)

t0 <- readRDS(paste(path, "Width_12.RDS", sep=""))
data <- data.frame(widths = c(t0$CIwidth, t0$ALBIwidth),
        CI = rep(c("Score", "FIESTA"), each = length(t0$CIwidth)))
s2 <- ggplot(data, aes(x=widths, fill=CI)) + geom_density(alpha=0.25) + ggtitle(expression(paste(h^2, " = " , 0.5, ", n = 1000")))+
 theme_bw() + xlab("CI width")+ ylab("") + theme(legend.position="") +
  theme(plot.title = element_text(hjust = 0.5))

t0 <- readRDS(paste(path, "Width_16.RDS", sep=""))
data <- data.frame(widths = c(t0$CIwidth, t0$ALBIwidth),
        Method = rep(c("Proposed", "FIESTA"), each = length(t0$CIwidth)))
s3 <- ggplot(data, aes(x=widths, fill=Method)) + geom_density(alpha=0.25) + ggtitle(expression(paste(h^2, " = " , 0.9, ", n = 1000")))+
 theme_bw() + xlab("CI width") + ylab("")  +
  theme(plot.title = element_text(hjust = 0.5))


dev.new()
plot_grid(s1, s2, s3 + theme(legend.position=""), nrow = 1)

plot_grid(s3 + theme(legend.position="bottom"), nrow = 1)





