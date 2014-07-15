source("~/repo/explodecomputr/R/functions.R")

tp <- c("cord", "F7", "15up", "antenatal", "FOM")
origfile <- paste("~/repo/mQTL-partitioning/filter_run1_gwas/normalise/ALN.tost.betas.", tp, ".OCT.ranktransformed.Rdata", sep="")
textfile <- paste("~/repo/methylation_residuals/data/ALN.tost.betas.", tp, ".OCT.ranktransformed.txt", sep="")

file.exists(origfile)


# Write files to text format - faster to read

for(i in 1:5)
{
	cat(i, "\n")
	load(origfile[i])
	write.table(tbetas, file=textfile[i], row=F, col=F, qu=F)
}



# Correspond aln to IDs in GRM

kidsfam <- "~/repo/mQTL-partitioning/grm/alspac_hm3_kids.fam"
mothersfam <- "~/repo/mQTL-partitioning/grm/alspac_hm3_mothers.fam"

load("ALN.tost.betas.FOM.OCT.ranktransformed.Rdata")
mt <- tbetas
mothers <- readFam(mothersfam)[,1:2]
mothers$aln <- gsub("[A-Z]", "", mothers$FID)
mothers <- subset(mothers, !duplicated(aln))
table(mothers$aln %in% colnames(mt))
table(colnames(mt) %in% mothers$aln)
mothersaln <- data.frame(aln=colnames(mt), index=1:ncol(mt), stringsAsFactors=FALSE)
mothers <- merge(mothersaln, mothers, by="aln", all.x=TRUE)
mothers <- mothers[order(mothers$index), ]
mindex <- is.na(mothers$FID)
table(mothers$aln == colnames(mt))
mothers$FID[mindex] <- paste(mothers$aln[mindex], "missing", sep="")
mothers$IID[mindex] <- paste(mothers$aln[mindex], "missing", sep="")
mothers <- subset(mothers, select=-c(aln, index))
write.table(mothers, file="~/repo/methylation_residuals/data/mothers.id", row=F, col=F, qu=F)

load("ALN.tost.betas.F7.OCT.ranktransformed.Rdata")
kt <- tbetas
kids <- readFam(kidsfam)[,1:2]
kids$aln <- gsub("[A-Z]", "", kids$FID)
kids <- subset(kids, !duplicated(aln))
table(kids$aln %in% colnames(mt))
table(colnames(mt) %in% kids$aln)
kidsaln <- data.frame(aln=colnames(mt), index=1:ncol(mt), stringsAsFactors=FALSE)
kids <- merge(kidsaln, kids, by="aln", all.x=TRUE)
kids <- kids[order(kids$index), ]
kindex <- is.na(kids$FID)
table(kids$aln == colnames(mt))
kids$FID[kindex] <- paste(kids$aln[kindex], "missing", sep="")
kids$IID[kindex] <- paste(kids$aln[kindex], "missing", sep="")
kids <- subset(kids, select=-c(aln, index))
write.table(kids, file="~/repo/methylation_residuals/data/kids.id", row=F, col=F, qu=F)


# Make parameter file

all(rownames(mt) == rownames(kt))

temp <- data.frame(cpg=rownames(kt), index=1:nrow(kt))
params <- expand.grid(cpg = rownames(kt), timepoint=tp)
params <- merge(params, temp, by="cpg", all=TRUE)
params <- params[order(params$timepoint, params$index),]
params$group <- "mothers"
params$group[params$timepoint %in% tp[1:3]] <- "kids"
params$cpg <- as.character(params$cpg)
params$timepoint <- as.character(params$timepoint)
params$group <- as.character(params$group)

save(params, file="~/repo/methylation_residuals/data/parameters.RData")


