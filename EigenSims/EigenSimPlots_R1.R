# -----------------------------------------------
# Create plots for simulation studies
# ----------------------------------------------
params <-  expand.grid(n=c(seq(200, 2000, length=9)), 
	h2 = c(.0001, .01, 0.1, 0.5), rho = c(0.1, 0.5, 0.95, 0.999))
fullResult <- NULL
fullResultSE <- NULL
fulln <- NULL
h2 <- NULL

path <- "/Results_R1/"

makeWidthPlot <- function(rho){
	fullResult <- NULL
	fullResultSE <- NULL
	fulln <- NULL
	h2 <- NULL
	# plot 1: h2 = 0.001, rho = 0.5, n varying
	inds <- which(params[,2] == 0.0001 & params[,3] == rho)
	n <- rep(0, length(inds))
	result <- rep(0, length(inds))
	result_se <- rep(0, length(inds))
	for(kk in 1:length(inds)){
		t0 <- readRDS(paste(path, "Coverage_", inds[kk], "_R1.RDS", sep=""))
		n[kk] <- params[inds[kk], 1]
		result[kk] <- mean(t0$CIwidth, na.rm=TRUE)
		result_se[kk] <- sd(t0$CIwidth, na.rm=TRUE)/sqrt(sum(!is.na(t0$CIwidth)))
		# cat(kk, "\n")
	}

	fullResult <- c(fullResult, result)
	fullResultSE <- c(fullResultSE, result_se)
	fulln <- c(fulln,n)
	h2 <- c(h2, rep(0.0001, length(result)))


	# plot 1: h2 = 0.1, rho = 0.5, n varying
	inds <- which(params[,2] == 0.01 & params[,3] == rho)
	n <- rep(0, length(inds))
	result <- rep(0, length(inds))
	result_se <- rep(0, length(inds))
	for(kk in 1:length(inds)){
		t0 <- readRDS(paste(path, "Coverage_", inds[kk], "_R1.RDS", sep=""))
		n[kk] <- params[inds[kk], 1]
		result[kk] <- mean(t0$CIwidth, na.rm=TRUE)
		result_se[kk] <- sd(t0$CIwidth, na.rm=TRUE)/sqrt(sum(!is.na(t0$CIwidth)))
		#cat(kk, "\n")
	}

	fullResult <- c(fullResult, result)
	fullResultSE <- c(fullResultSE, result_se)
	fulln <- c(fulln,n)
	h2 <- c(h2, rep(0.01, length(result)))

	inds <- which(params[,2] == 0.1 & params[,3] == rho)
	n <- rep(0, length(inds))
	result <- rep(0, length(inds))
	result_se <- rep(0, length(inds))
	for(kk in 1:length(inds)){
		t0 <- readRDS(paste(path, "Coverage_", inds[kk], "_R1.RDS", sep=""))
		n[kk] <- params[inds[kk], 1]
		result[kk] <- mean(t0$CIwidth, na.rm=TRUE)
		result_se[kk] <- sd(t0$CIwidth, na.rm=TRUE)/sqrt(sum(!is.na(t0$CIwidth)))
		#cat(kk, "\n")
	}

	fullResult <- c(fullResult, result)
	fullResultSE <- c(fullResultSE, result_se)
	fulln <- c(fulln,n)
	h2 <- c(h2, rep(0.1, length(result)))


	inds <- which(params[,2] == 0.5 & params[,3] == rho)
	n <- rep(0, length(inds))
	result <- rep(0, length(inds))
	result_se <- rep(0, length(inds))
	for(kk in 1:length(inds)){
		t0 <- readRDS(paste(path, "Coverage_", inds[kk], "_R1.RDS", sep=""))
		n[kk] <- params[inds[kk], 1]
		result[kk] <- mean(t0$CIwidth, na.rm=TRUE)
		result_se[kk] <- sd(t0$CIwidth, na.rm=TRUE)/sqrt(sum(!is.na(t0$CIwidth)))
		#cat(kk, "\n")
	}


	fullResult <- c(fullResult, result)
	fullResultSE <- c(fullResultSE, result_se)
	fulln <- c(fulln,n)
	h2 <- c(h2, rep(0.5, length(result)))

	dat <- data.frame("width" = fullResult, 
	"lower" = fullResult - 2*fullResultSE,
	"upper" = fullResult + 2*fullResultSE,
	 "n" = fulln, 
	"h2" = factor(h2))


	library(ggplot2)
	aOut <- ggplot(dat, aes(x = n, y = width, color = h2)) + geom_ribbon(aes(ymin=lower, ymax=upper, fill=h2)) + 
	  xlab("") + geom_line() + theme_bw() + 
	  scale_color_manual(values = c("grey0", "grey30",  "grey65", "grey90"))  +
	  scale_fill_manual(values = c("grey0", "grey30", "grey65", "grey90")) + ylim(c(0, 1)) +  theme(panel.grid.major = element_blank(),   # Remove major gridlines
	                                                                                                  panel.grid.minor = element_blank()) +
	  labs(fill = expression(h^2), color = expression(h^2)) + 
	  guides(fill = guide_legend(title = expression(h^2)),
	         color = guide_legend(title = expression(h^2)))
	return(aOut)

}


a1 <- makeWidthPlot(rho = 0.1)  + ylab("Average CI width") + ggtitle(expression(paste(rho, " = ", 0.1))) + theme(legend.position="",plot.title = element_text(hjust = 0.5),)
a1 + theme(legend.position="bottom")
a2 <- makeWidthPlot(rho = 0.5) + ylab("") + ggtitle(expression(paste(rho, " = ", 0.5))) + theme(legend.position="",plot.title = element_text(hjust = 0.5))+ theme(axis.title.y=element_blank())
a3 <- makeWidthPlot(rho = 0.95) + ylab("") + ggtitle(expression(paste(rho, " = ", 0.95))) + theme(legend.position="",plot.title = element_text(hjust = 0.5))+ theme(axis.title.y=element_blank())
a4 <- makeWidthPlot(rho = 0.999) + ylab("") + ggtitle(expression(paste(rho, " = ", 0.999))) + theme(legend.position="",plot.title = element_text(hjust = 0.5))+ theme(axis.title.y=element_blank())








makeCoveragePlot <- function(rho){
	fullResult <- NULL
	fullResultSE <- NULL

	fulln <- NULL
	h2 <- NULL


	# plot 1: h2 = 0.1, rho = 0.5, n varying
	inds <- which(params[,2] == 0.0001 & params[,3] == rho)
	n <- rep(0, length(inds))
	result <- rep(0, length(inds))
	result_se <- rep(0, length(inds))
	for(kk in 1:length(inds)){
	t0 <- readRDS(paste(path, "Coverage_", inds[kk], "_R1.RDS", sep=""))
	n[kk] <- params[inds[kk], 1]
	result[kk] <- mean(t0$scoreCoverage[,1])
	result_se[kk] <- sd(t0$scoreCoverage[,1], na.rm=TRUE)/sqrt(sum(!is.na(t0$scoreCoverage[,1])))
	cat(kk, "\n")
	}

	fullResult <- c(fullResult, result)
	fullResultSE <- c(fullResultSE, result_se)
	fulln <- c(fulln,n)
	h2 <- c(h2, rep(0.0001, length(result)))



	# plot 1: h2 = 0.1, rho = 0.5, n varying
	inds <- which(params[,2] == 0.01 & params[,3] == rho)
	n <- rep(0, length(inds))
	result <- rep(0, length(inds))
	result_se <- rep(0, length(inds))
	for(kk in 1:length(inds)){
	t0 <- readRDS(paste(path, "Coverage_", inds[kk], "_R1.RDS", sep=""))
	n[kk] <- params[inds[kk], 1]
	result[kk] <- mean(t0$scoreCoverage[,1])
	result_se[kk] <- sd(t0$scoreCoverage[,1], na.rm=TRUE)/sqrt(sum(!is.na(t0$scoreCoverage[,1])))
	cat(kk, "\n")
	}

	fullResult <- c(fullResult, result)
	fullResultSE <- c(fullResultSE, result_se)
	fulln <- c(fulln,n)
	h2 <- c(h2, rep(0.01, length(result)))

	inds <- which(params[,2] == 0.1 & params[,3] == rho)
	n <- rep(0, length(inds))
	result <- rep(0, length(inds))
	result_se <- rep(0, length(inds))
	for(kk in 1:length(inds)){
	t0 <- readRDS(paste(path, "Coverage_", inds[kk], "_R1.RDS", sep=""))
	n[kk] <- params[inds[kk], 1]
	result[kk] <- mean(t0$scoreCoverage[,1])
	result_se[kk] <- sd(t0$scoreCoverage[,1], na.rm=TRUE)/sqrt(sum(!is.na(t0$scoreCoverage[,1])))
	cat(kk, "\n")
	}

	fullResult <- c(fullResult, result)
	fullResultSE <- c(fullResultSE, result_se)
	fulln <- c(fulln,n)
	h2 <- c(h2, rep(0.1, length(result)))


	inds <- which(params[,2] == 0.5 & params[,3] == rho)
	n <- rep(0, length(inds))
	result <- rep(0, length(inds))
	result_se <- rep(0, length(inds))
	for(kk in 1:length(inds)){
	t0 <- readRDS(paste(path, "Coverage_", inds[kk], "_R1.RDS", sep=""))
	n[kk] <- params[inds[kk], 1]
	result[kk] <- mean(t0$scoreCoverage[,1])
	result_se[kk] <- sd(t0$scoreCoverage[,1], na.rm=TRUE)/sqrt(sum(!is.na(t0$scoreCoverage[,1])))
	cat(kk, "\n")
	}


	fullResult <- c(fullResult, result)
	fullResultSE <- c(fullResultSE, result_se)
	fulln <- c(fulln,n)
	h2 <- c(h2, rep(0.5, length(result)))

	dat <- data.frame("width" = fullResult, 
	"lower" = fullResult - 2*fullResultSE,
	"upper" = fullResult + 2*fullResultSE,
	"n" = fulln, 
	"h2" = factor(h2))


	library(ggplot2)
	sOut <- ggplot(dat, aes(x = n, y = width, color = h2)) + 
	geom_errorbar(aes(x = n, ymin=lower, ymax=upper), size= 0.6, width = 0, position = position_dodge(120)) + 
	#geom_ribbon(aes(x = n, ymin=lower, ymax=upper, fill=h2), alpha = 0.3, color=NA) + 
	#geom_line(alpha = 0) + 
	theme_bw() + scale_x_continuous(breaks = unique(fulln)[c(1,3,5,7,9, 11, 13, 15)]) + 
	  scale_color_manual(values = c("grey0", "grey30",  "grey65", "grey80")) +  theme(panel.grid.major = element_blank(),   # Remove major gridlines
	                                                                                  panel.grid.minor = element_blank()) +
	ylim(c(0.938, 0.962))+ geom_hline(yintercept = 0.95) #+ theme(axis.text.x = element_text(angle = 25, vjust = 1, hjust=1))
	return(sOut)


}


s1 <- makeCoveragePlot(rho = 0.1)  + ylab("Coverage probability") + theme(legend.position="", plot.title = element_text(hjust = 0.5))
s2 <- makeCoveragePlot(rho = 0.5)  + ylab("") + theme(legend.position="",plot.title = element_text(hjust = 0.5)) + theme(axis.title.y=element_blank())
s3 <- makeCoveragePlot(rho = 0.95)  + ylab("") + theme(legend.position="",plot.title = element_text(hjust = 0.5)) + theme(axis.title.y=element_blank())
s4 <- makeCoveragePlot(rho = 0.999)  + ylab("") + theme(legend.position="",plot.title = element_text(hjust = 0.5)) + theme(axis.title.y=element_blank())



library(cowplot)
pdf(file="EigenSims.pdf", width=9, height=4.3)
plot_grid(a1, a2,   a3,  a4,  NULL , NULL, NULL, NULL,  s1, s2, s3, s4, nrow=3, 
          rel_heights=c(1, -0.1, 1), rel_widths=c(1, 1,  1, 1), align="h")
dev.off()

t00 <- plot_grid(a1, a2,   a3,  a4,  NULL , NULL, NULL, NULL,  s1, s2, s3, s4, nrow=3, rel_heights=c(1, -0.1, 1), rel_widths=c(1, 1,  1, 1))
t000 <- get_legend(t00  + theme(legend.position="bottom"))