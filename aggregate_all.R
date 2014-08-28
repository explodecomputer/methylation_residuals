setwd("~/repo/methylation_residuals/res")

a1 <- list()
a2 <- list()
a3 <- list()
a4 <- list()
for(i in 1:122)
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

save(dat, file="../hsq_aries.RData")


with(dat, tapply(hsq, timepoint, mean))
with(dat, tapply(hsq, timepoint, sd))

library(ggplot2)

ggplot(dat, aes(x=timepoint, y=hsq)) +
geom_boxplot()


summary(lm(hsq ~ timepoint, dat))
summary(lm(hsq ~ as.numeric(timepoint), dat))
summary(lm(hsq ~ as.numeric(timepoint) + I(as.numeric(timepoint)^2), dat))

x <- 0:6
y <- 0.05 + x*0.0986 + x^2*-0.015
plot(y ~ x)


pdf("~/repo/methylation_residuals/images/h2_estimates.pdf")
boxplot(hsq ~ as.numeric(timepoint), dat)
abline(lm(hsq ~ as.numeric(timepoint), dat), col="blue")
lines(y ~ x, col="red")
dev.off()



geom_density(aes(fill=timepoint)) +
scale_x_log10()

stat_density(aes(ymax = ..density..,  ymin = -..density..),
    fill = "grey50", colour = "grey50",
    geom = "ribbon", position = "identity")
