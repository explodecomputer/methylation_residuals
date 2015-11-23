load("../methylation_data/data/correlation_F7_15up.RData")
cpgs <- as.character(subset(dat, variable=="QTLs")$cpg)
load("data/variances.RData")

dat$qtl <- "All"
dat$qtl[dat$cpg %in% cpgs] <- "With QTL"
with(dat, tapply(vars, list(timepoint, qtl), median))
with(dat, tapply(vars, list(timepoint, qtl), mean))

library(ggplot2)
ggplot(dat, aes(x=timepoint,y=vars)) +
geom_boxplot(aes(fill=qtl)) + 
labs(y="Probe variance (-log10)", x="Time point", fill="CpG site") +
scale_y_continuous(limits=c(1e-7, 1), trans="log10")
ggsave(file="images/variances_mqtl.pdf")

