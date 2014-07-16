load("~/repo/methylation_residuals/data/parameters.RData")

makePhen <- function(allphenfile, idfile, row, phenfile)
{
	cmd <- paste("head -n ", row, " ", allphenfile, " | tail -n 1 | tr ' ' '\n' > ", phenfile, ".temp", sep="")
	system(cmd)
	cmd <- paste("paste -d ' ' ", idfile, " ", phenfile, ".temp > ", phenfile, sep="")
	system(cmd)
	system(paste("rm ", phenfile, ".temp", sep=""))
}

runGcta <- function(grmfile, phenfile, outfile, flags="")
{
	cmd <- paste("gcta64 ", flags, " --grm ", grmfile, " --reml --reml-no-lrt --pheno ", phenfile, " --reml-pred-rand --out ", outfile, sep="")
	system(cmd)
}

arguments <- commandArgs(T)
jid <- as.numeric(arguments[1])
nrun <- as.numeric(arguments[2])

first <- (jid - 1) * nrun + 1
last <- min(nrow(params), jid * nrun)

print(c(first, last))

# grmroot <- "~/repo/methylation_residuals/data/alspac_hm3_GROUP_unrelated_m"
# allphenroot <- "~/repo/methylation_residuals/data/ALN.tost.betas.TIMEPOINT.OCT.ranktransformed.txt"
# idroot <- "~/repo/methylation_residuals/data/GROUP.id"
# phenroot <- "~/repo/methylation_residuals/results/TIMEPOINT/CPG.phen"
# outroot <- "~/repo/methylation_residuals/results/TIMEPOINT/CPG"

for(i in first:last)
{
	grmfile <- paste("~/repo/methylation_residuals/data/", params$timepoint[i], sep="")
	allphenfile <- paste("~/repo/methylation_residuals/data/", params$timepoint[i], ".phen", sep="")
	idfile <- paste("~/repo/methylation_residuals/data/", params$timepoint[i], ".id", sep="")
	outfile <- paste("~/repo/methylation_residuals/results/", params$timepoint[i], "/", params$cpg[i], sep="")
	phenfile <- paste("~/repo/methylation_residuals/results/", params$timepoint[i], "/", params$cpg[i], ".phen", sep="")

	makePhen(allphenfile, idfile, params$index[i], phenfile)
	runGcta(grmfile, phenfile, outfile)
}