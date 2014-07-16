source("~/repo/explodecomputr/R/functions.R")
makeIds <- function(famfile, timepoint, tbetas)
{
	people <- readFam(famfile)[,1:2]
	people$aln <- gsub("[A-Z]", "", people$FID)
	people <- subset(people, !duplicated(aln))
	print(table(people$aln %in% colnames(tbetas)))
	print(table(colnames(tbetas) %in% people$aln))
	peoplealn <- data.frame(aln=colnames(tbetas), index=1:ncol(tbetas), stringsAsFactors=FALSE)
	people <- merge(peoplealn, people, by="aln", all.x=TRUE)
	people <- people[order(people$index), ]
	pindex <- is.na(people$FID)
	print(table(people$aln == colnames(tbetas)))
	people$FID[pindex] <- paste(people$aln[pindex], "missing", sep="")
	people$IID[pindex] <- paste(people$aln[pindex], "missing", sep="")
	people <- subset(people, select=-c(aln, index))
	write.table(people, file=paste("~/repo/methylation_residuals/data/", timepoint, ".id", sep=""), row=F, col=F, qu=F)
}

tp <- c("cord", "F7", "15up", "antenatal", "FOM")
famfile <- c("~/repo/mQTL-partitioning/grm/alspac_hm3_kids.fam", "~/repo/mQTL-partitioning/grm/alspac_hm3_kids.fam", "~/repo/mQTL-partitioning/grm/alspac_hm3_kids.fam", "~/repo/mQTL-partitioning/grm/alspac_hm3_mothers.fam", "~/repo/mQTL-partitioning/grm/alspac_hm3_mothers.fam")
origfile <- paste("~/repo/mQTL-partitioning/filter_run1_gwas/normalise/ALN.tost.betas.", tp, ".OCT.ranktransformed.Rdata", sep="")
textfile <- paste("~/repo/methylation_residuals/data/", tp, ".phen", sep="")
grmfile <- c("~/repo/mQTL-partitioning/grm/alspac_hm3_kids_unrelated", "~/repo/mQTL-partitioning/grm/alspac_hm3_kids_unrelated", "~/repo/mQTL-partitioning/grm/alspac_hm3_kids_unrelated", "~/repo/mQTL-partitioning/grm/alspac_hm3_mothers_unrelated", "~/repo/mQTL-partitioning/grm/alspac_hm3_mothers_unrelated")

file.exists(origfile)
file.exists(famfile)

# Write files to text format - faster to parse
# Make new ID files

for(i in 2:5)
{
	cat(i, "\n")
	load(origfile[i])
	makeIds(famfile[i], tp[i], tbetas)
	cmd <- paste("gcta64 --grm ", grmfile[i], " --keep ", paste("~/repo/methylation_residuals/data/", tp[i], ".id", sep=""), " --make-grm --out ", paste("~/repo/methylation_residuals/data/", tp[i], sep=""), sep="")
	system(cmd)
	# write.table(tbetas, file=textfile[i], row=F, col=F, qu=F)
}


# Make parameter file

temp <- data.frame(cpg=rownames(tbetas), index=1:nrow(tbetas))
params <- expand.grid(cpg = rownames(tbetas), timepoint=tp)
params <- merge(params, temp, by="cpg", all=TRUE)
params <- params[order(params$timepoint, params$index),]
params$group <- "mothers"
params$group[params$timepoint %in% tp[1:3]] <- "kids"
params$cpg <- as.character(params$cpg)
params$timepoint <- as.character(params$timepoint)
params$group <- as.character(params$group)

save(params, file="~/repo/methylation_residuals/data/parameters.RData")


