# objectives:

# what is the average snp heritability
# - cis / trans
# how does this change with age
# what proportion of estimated genetic variation is captured by SNPs

# is there a difference in SNP heritability for different genomic features




library(ggplot2)
library(dplyr)
library(reshape2)


load("~/repo/methylation_residuals/data/datl_v2.RData")
load("~/repo/mQTL-partitioning/filter_run5_gwas/data/condanal_results.RData")
naeem <- subset(condres, select=c(CPG, keep))
datl <- merge(datl, naeem, by="CPG")

# total h2 divided into cis and trans for each time point

with(datl, tapply(value, list(timepoint, variable, experiment), mean))

p1 <- group_by(datl, timepoint, variable, experiment)
p1 <- summarise(p1, ave = mean(value), med=median(value), se = sqrt(var(value) / length(value)), c05 = quantile(value, c(0.25)), c95 = quantile(value, 0.75))

# ggplot(subset(p1, experiment=="SNP heritability" & variable != "Total"), aes(y=ave, x=timepoint)) +
# geom_bar(stat="identity", aes(fill=variable)) +
# labs(y="SNP heritability", fill="Variance\ncomponent", x="")


# proportion explained h2 divided into cis and trans

fr <- subset(p1, experiment=="mQTL")
fr$prop <- subset(p1, experiment=="mQTL")$ave / subset(p1, experiment=="SNP heritability")$ave
fr$ave <- subset(p1, experiment=="SNP heritability")$ave - subset(p1, experiment=="mQTL")$ave
temp <- subset(p1, variable != "Total" & experiment == "mQTL", select=c(timepoint, variable, ave))
fr <- subset(fr, variable != "Total", select=c(timepoint, variable, ave))
fr$expl <- "Unexplained"
temp$expl <- "Explained"

fr <- rbind(fr, temp)
fr$lab <- factor(paste(fr$variable, fr$expl))
fr$lab <- factor(fr$lab, levels=levels(fr$lab)[c(2,1,4,3)])


ggplot(fr, aes(y=ave, x=timepoint, order=order(timepoint))) +
geom_bar(stat="identity", aes(fill=lab, order=order(lab))) +
scale_fill_brewer(type="qual", palette=3) +
labs(y="SNP heritability", x="", fill="")
ggsave(file="~/repo/methylation_residuals/images/h2_partition_cistrans.pdf")



# is there a difference in SNP heritability for different genomic features

summary(lm(value ~ Feature + timepoint, subset(datl, variable=="Total"&experiment=="SNP heritability")))
summary(lm(value ~ Feature + timepoint, subset(datl, variable=="Cis"&experiment=="SNP heritability")))
summary(lm(value ~ Feature + timepoint, subset(datl, variable=="Trans"&experiment=="SNP heritability")))
summary(lm(value ~ Feature * as.numeric(timepoint), subset(datl, variable=="Cis"&experiment=="SNP heritability")))
summary(lm(value ~ Feature * timepoint, subset(datl, variable=="Cis"&experiment=="SNP heritability")))
# No interaction between h2 and time point, just take childhood for simplicity


p2 <- group_by(subset(datl, experiment == "SNP heritability" & timepoint == "Childhood" & variable == "Cis"), Feature)
p2 <- summarise(p2, ave = mean(value), med=median(value), se = sqrt(var(value) / length(value)), c05 = quantile(value, c(0.25)), c95 = quantile(value, 0.75))
p2$Feature <- factor(p2$Feature, levels=levels(p2$Feature)[order(p2$ave)])

p3 <- group_by(subset(datl, experiment == "SNP heritability" & timepoint == "Childhood"), Feature, variable)
p3 <- summarise(p3, ave = mean(value), med=median(value), se = sqrt(var(value) / length(value)), c05 = quantile(value, c(0.25)), c95 = quantile(value, 0.75))
p3$Feature <- factor(p3$Feature, levels=levels(p3$Feature)[order(p3$ave[p3$variable=="Cis"])])

ggplot(p2, aes(y=ave, x=Feature)) +
geom_bar(stat="identity") +
geom_errorbar(aes(ymax=ave+se, ymin=ave-se, width=0.2)) +
geom_hline(data=subset(p1, experiment=="SNP heritability" & timepoint == "Childhood" & variable == "Cis"), aes(yintercept=ave)) +
theme(axis.text.x=element_text(angle=90, hjust=1, vjust=0.5), axis.ticks.x=element_blank()) +
labs(x="", y="Cis SNP heritability")
ggsave(file="~/repo/methylation_residuals/images/cis_h2_features.pdf")

ggplot(subset(p3, variable!="Total"), aes(y=ave, x=Feature)) +
geom_bar(stat="identity") +
geom_errorbar(aes(ymax=ave+se, ymin=ave-se, width=0.2)) +
geom_hline(data=subset(p1, experiment=="SNP heritability" & timepoint == "Childhood" & variable != "Total"), aes(yintercept=ave)) +
facet_grid(variable ~ ., scale="free_y") +
theme(axis.text.x=element_text(angle=90, hjust=1, vjust=0.5), axis.ticks.x=element_blank()) +
labs(x="", y="Mean SNP heritability")
ggsave(file="~/repo/methylation_residuals/images/h2_features.pdf")



temp1 <- subset(datl, variable=="Cis", select=c(CPG, timepoint, value, experiment, keep))
temp1 <- temp1[order(temp1$CPG, temp1$timepoint, temp1$experiment),]

temp1a <- subset(temp1, experiment=="mQTL")
temp1b <- subset(temp1, experiment!="mQTL")

plot(temp1a$value ~ temp1b$value, alpha=0.1)
temp1a$value2 <- temp1b$value

ggplot(temp1a, aes(x=value2, y=value)) +
geom_point(alpha=0.2) +
facet_grid(keep ~ .) +
labs(x="Cis SNP h2", y="Cis meQTL h2")
ggsave(file="~/repo/methylation_residuals/images/cis_naeem.png")
