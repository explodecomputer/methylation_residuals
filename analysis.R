library(ggplot2)
library(plyr)

load("~/repo/methylation_residuals/hsq_aries.RData")

levels(dat$timepoint) <- c("Birth", "Childhood", "Adolescence", "Pregnancy", "Middle age")

with(dat, tapply(hsq, timepoint, mean))
with(dat, tapply(hsq, timepoint, sd))

summary(lm(hsq ~ timepoint, dat))
summary(lm(hsq ~ as.numeric(timepoint), dat))
summary(lm(hsq ~ as.numeric(timepoint) + I(as.numeric(timepoint)^2), dat))


s <- ddply(dat, .(timepoint), summarise, ave = mean(hsq), se = sd(hsq) / sqrt(mean(n)))

ggplot(dat, aes(x=hsq)) +
geom_density() +
facet_grid(timepoint ~ .) +
geom_vline(data=s, aes(xintercept=ave), colour="red") +
labs(y = "Density", x="Heritability estimate")
ggsave(file="~/repo/methylation_residuals/images/distributions.pdf", width=5, height=7)

dat2 <- subset(dat, hsq != 1e-6)
s2 <- ddply(dat2, .(timepoint), summarise, ave = mean(hsq))

ggplot(dat2, aes(x=hsq)) +
geom_density() +
facet_grid(timepoint ~ .) +
geom_vline(data=s2, aes(xintercept=ave), colour="red")



ggplot(dat, aes(x=se)) +
geom_histogram() +
facet_grid(timepoint ~ .)


x <- 0:6
y <- 0.05 + x*0.0986 + x^2*-0.015
plot(y ~ x)


pdf("~/repo/methylation_residuals/images/h2_estimates.pdf")
boxplot(hsq ~ as.numeric(timepoint), dat)
abline(lm(hsq ~ as.numeric(timepoint), dat), col="blue")
lines(y ~ x, col="red")
dev.off()


a <- ddply(dat, .(timepoint), function(x) quantile(x$hsq))
b <- ddply(dat2, .(timepoint), function(x) quantile(x$hsq))
write.table(a, file="~/repo/methylation_residuals/images/h2_estimates.csv", row=F, col=T, qu=F, sep=",")
write.table(b, file="~/repo/methylation_residuals/images/h2_estimates_non0.csv", row=F, col=T, qu=F, sep=",")


