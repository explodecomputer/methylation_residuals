library(plyr)
library(ggplot2)
library(matrixStats)

tp <- c("cord", "F7", "15up", "antenatal", "FOM")

l <- list()
for(i in 1:length(tp))
{
	message(tp[i])
	load(paste0("~/data/deprecated/aries/releasev2_apr2014/ALN.tost.betas.", tp[i], ".releasev2_apr2014.Rdata"))
	a <- rowQuantiles(tbetas, prob=c(0.25, 0.75))
	cpg <- read.csv(paste0("~/repo/mQTL-partitioning/filter_run5_gwas/data/", tp[i], "_sorted.gwcannotated"))
	index <- rownames(tbetas) %in% cpg$CPG
	b <- a[,2] - a[,1]
	l[[i]] <- data.frame(a, genetic=index, timepoint=tp[i])
}

dat <- rbind.fill(l)
levels(dat$timepoint) <- c("Birth", "Childhood", "Adolescence", "Pregnancy", "Middle age")

with(dat, tapply(dat$X75-dat$X25, list(timepoint, genetic), mean))
save(dat, file="~/repo/methylation_residuals/iqr.RData")

ggplot(d, aes(x=)) +
geom_density(aes(fill=index), alpha=0.2) +
facet_wrap(~ timepoint)



