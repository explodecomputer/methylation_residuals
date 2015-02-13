source("~/repo/explodecomputr/R/functions.R")
makeIds <- function(famfile, timepoint, M2)
{
	cat("Making IDs\n")
	people <- readFam(famfile)[,1:2]
	people$aln <- gsub("[A-Z]", "", people$FID)
	people <- subset(people, !duplicated(aln))
	print(table(people$aln %in% colnames(M2)))
	print(table(colnames(M2) %in% people$aln))
	peoplealn <- data.frame(aln=colnames(M2), index=1:ncol(M2), stringsAsFactors=FALSE)
	people <- merge(peoplealn, people, by="aln", all.x=TRUE)
	people <- people[order(people$index), ]
	pindex <- is.na(people$FID)
	print(table(people$aln == colnames(M2)))
	people$FID[pindex] <- paste(people$aln[pindex], "missing", sep="")
	people$IID[pindex] <- paste(people$aln[pindex], "missing", sep="")
	people <- subset(people, select=-c(aln, index))
	write.table(people, file=paste("~/repo/methylation_residuals/data/", timepoint, ".id", sep=""), row=F, col=F, qu=F)
}

tp <- c("cord", "F7", "15up", "antenatal", "FOM")
famfile <- c("~/repo/mQTL-partitioning/grm/alspac_hm3_kids.fam", "~/repo/mQTL-partitioning/grm/alspac_hm3_kids.fam", "~/repo/mQTL-partitioning/grm/alspac_hm3_kids.fam", "~/repo/mQTL-partitioning/grm/alspac_hm3_mothers.fam", "~/repo/mQTL-partitioning/grm/alspac_hm3_mothers.fam")
origfile <- paste("~/data/aries/adjusted_all.v1.", tp, ".RData", sep="")
textfile <- paste("~/repo/methylation_residuals/data/", tp, "norm.phen", sep="")
grmfile <- c("~/repo/mQTL-partitioning/grm/alspac_hm3_kids_unrelated", "~/repo/mQTL-partitioning/grm/alspac_hm3_kids_unrelated", "~/repo/mQTL-partitioning/grm/alspac_hm3_kids_unrelated", "~/repo/mQTL-partitioning/grm/alspac_hm3_mothers_unrelated", "~/repo/mQTL-partitioning/grm/alspac_hm3_mothers_unrelated")

file.exists(origfile)
file.exists(famfile)

# Write files to text format - faster to parse
# Make new ID files

for(i in 1:5)
{
	cat(i, "\n")
	load(origfile[i])
	M2 <- t(M2)
	makeIds(famfile[i], tp[i], M2)
	cmd <- paste("gcta64 --grm ", grmfile[i], " --keep ", paste("~/repo/methylation_residuals/data/", tp[i], ".id", sep=""), " --make-grm --grm-cutoff 0.05 --out ", paste("~/repo/methylation_residuals/data/", tp[i], sep=""), sep="")
	system(cmd)
	write.table(M2, file=textfile[i], row=F, col=F, qu=F)
}


# Make parameter file
temp <- data.frame(cpg=rownames(M2), index=1:nrow(M2))
params <- expand.grid(cpg = rownames(M2), timepoint=tp)
params <- merge(params, temp, by="cpg", all=TRUE)
params <- params[order(params$timepoint, params$index),]
params$group <- "mothers"
params$group[params$timepoint %in% tp[1:3]] <- "kids"
params$cpg <- as.character(params$cpg)
params$timepoint <- as.character(params$timepoint)
params$group <- as.character(params$group)


# Make SNP info file
# Mothers bimfile is identical to kids bimfile

snpinfo <- readBim("~/repo/mQTL-partitioning/grm/alspac_hm3_kids.bim")[,c(1,2,4)]


load("~/repo/mQTL-partitioning/filter_run1_gwas/probe_info/filtered_probe_list.RData")
probeinfo <- subset(probeinfo, CHR %in% 1:22)

tab <- list()
for(i in 1:22)
{
	cat(i, "\n")
	r <- range(subset(snpinfo, CHR == i)$BP)
	tab[[i]] <- subset(probeinfo, CHR37 == i & (COORDINATE_36 < r[1] | COORDINATE_36 > r[2]))$TargetID
}
probeinfo <- subset(probeinfo, ! TargetID %in% unlist(tab))
params <- subset(params, cpg %in% probeinfo$TargetID)

save(params, probeinfo, snpinfo, file="~/repo/methylation_residuals/data/parameters.RData")
