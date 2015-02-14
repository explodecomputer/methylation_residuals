library(dplyr)

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

qtldat <- rbind_all(b)
index <- which(is.na(qtldat$varexp))
mod <- interpolateMod(800)
qtldat$varexp[index] <- inferInterpolate(qtldat$pval[index], mod)

qtldat <- subset(qtldat, select=c(CPG, timepoint, Trans, varexp))
levels(qtldat$Trans) <- c("hsq1", "hsq2")

qd <- group_by(qtldat, CPG, timepoint, Trans)
qd <- summarise(qd, varexp=sum(varexp))

qd2 <- group_by(qd, CPG, timepoint)
qd2 <- summarise(qd2, varexp=sum(varexp))

qd2$Trans <- "hsq"

qd2 <- subset(qd2, select=c(CPG, timepoint, Trans, varexp))
qd <- rbind(qd, qd2)


save(qd, file="~/repo/methylation_residuals/data/varexp.RData")