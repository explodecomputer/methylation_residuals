library(ggplot2)
library(reshape2)
library(plyr)

check_pc_sample_type <- function(timepoint, npc)
{
	load(paste0("~/data/aries/adjusted_all.v1.", timepoint, ".RData"))
	index <- seq(1, ncol(M2), len=20000)
	pc <- prcomp(M2[,index])

	ids <- read.csv("~/data/deprecated/aries/releasev2_apr2014/ARIES.manifest.releasev2_apr2014.txt")

	table(ids$sample_type)
	levels(ids$sample_type) <- c("Blood spots", "White cells", "PBL", "White cells", "White cells", "Whole blood", "Whole blood")

	id <- subset(ids, time_point==timepoint)
	id <- id[match(rownames(M2), id$ALN), ]


	x1 <- melt(pc$x[,c(1,3,5,7)])
	x2 <- melt(pc$x[,c(2,4,6,8)])
	names(x1) <- paste0(names(x1), ".x")
	names(x2) <- paste0(names(x2), ".y")
	x <- cbind(x1, x2)
	x$pc <- paste(x$Var2.x, "vs", x$Var2.y)
	id <- subset(id, select=c(ALN, sample_type))

	x <- merge(x, id, by.x="Var1.x", by.y="ALN")
	x$timepoint <- timepoint


	id <- data.frame(id, pc$x[,1:npcs])
	id$sample_type <- factor(id$sample_type)


	pval <- array(0, npc)
	for(i in 1:npc)
	{
		form <- as.formula(paste0("PC", i, " ~ sample_type"))
		pval[i] <- anova(lm(form, id))$P[1]
	}

	d <- data.frame(Timepoint = timepoint, PC = paste0("PC",1:npc), P.value=pval)
	return(list(d=d, x=x))
}

tp <- c("cord", "F7", "antenatal", "FOM")
l <- list()
for(i in 1:length(tp))
{
	message(tp[i])
	l[[i]] <- check_pc_sample_type(tp[i], 10)
}



x <- rbind.fill(lapply(l, function(x) x$x))

x$timepoint <- as.factor(x$timepoint)
levels(x$timepoint) <- c("Pregnancy", "Birth", "Childhood", "Middle age")
x$timepoint <- factor(x$timepoint, levels=c("Birth", "Childhood", "Pregnancy", "Middle age"))

ggplot(x, aes(value.x, value.y)) +
geom_point(aes(colour=sample_type)) +
facet_grid(timepoint ~ pc, scale="free") +
scale_colour_brewer(type="qual") +
labs(x="PC value", y="PC value", colour="Sample type")
ggsave(file="~/repo/methylation_residuals/images/pca_methylation_sampletype.pdf")

