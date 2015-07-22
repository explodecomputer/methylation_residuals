library(dplyr)
library(reshape2)


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

harmoniseDat <- function(qd, dat, annot)
{
	levels(dat$timepoint) <- c("Birth", "Childhood", "Pregnancy", "Middle age", "Adolescence")
	qd$timepoint <- as.factor(qd$timepoint)
	levels(qd$timepoint) <- c("Adolescence", "Pregnancy", "Birth", "Childhood", "Middle age")
	dat$hsq <- dat$hsq1 + dat$hsq2
	dat2 <- melt(subset(dat, select=c(CPG, timepoint, hsq, hsq1, hsq2)), id.vars=c("CPG", "timepoint"))
	names(qd) <- names(dat2)

	qd$experiment <- "mQTL"
	dat2$experiment <- "SNP heritability"

	datl <- rbind(dat2, qd)
	levels(datl$variable) <- c("Total", "Cis", "Trans")
	datl$timepoint <- factor(datl$timepoint, levels=c("Birth", "Childhood", "Adolescence", "Pregnancy", "Middle age"))

	# Put 0s for cpgs that have no mqtl

	datl$code <- with(datl, paste(CPG, timepoint, variable))
	a <- subset(datl, experiment=="SNP heritability")
	b <- subset(datl, experiment=="mQTL", select=c(value, experiment, code))
	ab <- merge(a, b, by="code", all.x=TRUE)
	ab$value.y[is.na(ab$value.y)] <- 0

	b <- subset(ab, select=c(CPG, timepoint, variable, value.y))
	a <- subset(a, select=-c(code))
	b$experiment <- "mQTL"
	names(b) <- names(a)
	datl <- rbind(a, b)

	datl <- merge(datl, annot, by="CPG", all.x=TRUE)
	levels(datl$Feature) <- c("3' UTR", "5' UTR", "Coding", "Downstream", "Essential splice site", "Intergenic", "Intronic", "Regulatory region", "Splice site", "Upstream", "Within mature miRNA", "Within non coding gene")
	datl$Feature[datl$Feature=="Essential splice site"] <- "Splice site"
	datl$Feature[datl$Feature=="Within mature miRNA"] <- "Within non coding gene"
	datl$Feature <- factor(datl$Feature)

	return(datl)
}


wd <- "~/repo/mQTL-partitioning/filter_run5_gwas/data"
timepoint <- c("15up", "antenatal", "cord", "F7", "FOM")
filename <- file.path(wd, paste(timepoint, ".txt", sep=""))


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
qtldat <- subset(qtldat, pval < 1.25e-13)
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


annot <- read.table("~/repo/methylation_residuals/data/annot.txt", sep="\t")
annot <- subset(annot, select=c(V1, V8))
names(annot) <- c("CPG", "Feature")

load("~/repo/methylation_residuals/hsq_ct_aries.RData")
dat_ct <- dat
load("~/repo/methylation_residuals/hsq_ct_chr_aries.RData")
dat_ct_chr <- dat


datl <- harmoniseDat(qd, dat_ct, annot)
datl_chr <- harmoniseDat(qd, dat_ct_chr, annot)

save(datl, file="~/repo/methylation_residuals/data/datl_v2.RData")
save(datl_chr, file="~/repo/methylation_residuals/data/datl_chr_v2.RData")
save(qd, file="~/repo/methylation_residuals/data/varexp.RData")