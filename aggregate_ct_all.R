setwd("~/repo/methylation_residuals/res")

a1 <- list()
a2 <- list()
a3 <- list()
a4 <- list()
for(i in 1:91)
{
	cat(i, "\n")
	load(paste("allres", i, ".RData", sep=""))
	a1[[i]] <- unlist(l1)
	a2[[i]] <- unlist(l2)
	a3[[i]] <- unlist(l3)
	a4[[i]] <- do.call(rbind, l4)
}

load("../data/parameters.RData")


f1 <- unlist(a1)
f2 <- unlist(a2)
f3 <- unlist(a3)
f4 <- data.frame(do.call(rbind, a4))

names(f4) <- c("CPG", "timepoint")

dat <- data.frame(hsq=f1, se=f2, n=f3, f4)
dat$timepoint <- factor(dat$timepoint, levels=c("cord", "F7", "15up", "antenatal", "FOM"))

save(dat, file="../hsq_ct_aries.RData")


