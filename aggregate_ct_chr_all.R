setwd("~/repo/methylation_residuals/res_ct_chr")

a1 <- list()
a2 <- list()
a3 <- list()
a4 <- list()
a5 <- list()
a6 <- list()
a7 <- list()
for(i in 1:91)
{
	cat(i, "\n")
	fn <- paste("allres", i, ".RData", sep="")
	if(file.exists(fn))
	{
		load(fn)
		a1[[i]] <- unlist(l1)
		a2[[i]] <- unlist(l2)
		a3[[i]] <- unlist(l3)
		a4[[i]] <- unlist(l4)
		a5[[i]] <- unlist(l5)
		a6[[i]] <- unlist(l6)
		a7[[i]] <- do.call(rbind, l7)
	} else {
		message(fn, "missing")
	}
}

load("../data/parameters.RData")


f1 <- unlist(a1)
f2 <- unlist(a2)
f3 <- unlist(a3)
f4 <- unlist(a4)
f5 <- unlist(a5)
f6 <- unlist(a6)
f7 <- data.frame(do.call(rbind, a7))

names(f7) <- c("CPG", "timepoint")

dat <- data.frame(hsq1=f1, se1=f2, hsq2=f3, se2=f4, n=f5, f7)
dat$timepoint <- factor(dat$timepoint, levels=c("cord", "F7", "15up", "antenatal", "FOM"))

save(dat, file="../hsq_ct_chr_aries.RData")


