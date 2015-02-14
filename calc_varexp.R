library(plyr)

wd <- "~/repo/mQTL-partitioning/filter_run5_gwas/data"
timepoint <- c("15up", "antenatal", "cord", "F7", "FOM")
filename <- file.path(wd, paste(timepoint, ".txt", sep=""))

calcHsq <- function(pval, n)
{
	fval <- qf(pval, 1, n-1, low=FALSE)
	h2 <- fval / (fval + n - 1)
	return(h2)
}

interpolateMod <- function(n)
{
	p <- 1:300
	pvals <- 10^-p
	hs <- calcHsq(pvals, n)
	mod <- lm(hs ~ p + I(p^2) + I(p^3))
	return(mod)
}

inferInterpolate <- function(pval, mod)
{
	new <- data.frame(p = -log10(pval))
	predict(mod, new)
}


b <- list()
for(i in 1:length(timepoint))
{
	cat(i, "\n")
	a <- read.table(filename[i], he=T, sep=",")
	a <- subset(a, select=c(SNP, CPG, pval, Trans))
	a$timepoint <- timepoint[i]
	a$pval[a$pval==0] <- min(a$pval[a$pval!=0], na.rm=T)
	a$varexp <- calcHsq(a$pval, 800)
	b[[i]] <- a
}

qtldat <- rbind.fill(b)
index <- which(is.na(qtldat$varexp))
mod <- interpolateMod(800)
qtldat$varexp[index] <- inferInterpolate(qtldat$pval[index], mod)

save(qtldat, file="~/repo/methylation_residuals/data/varexp.RData")