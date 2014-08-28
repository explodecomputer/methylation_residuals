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

readPreds <- function(rootname)
{
	filename <- paste(rootname, ".indi.blp", sep="")
	if(file.exists(filename))
	{
		a <- read.table(filename, colClass=c("character", "character", "numeric", "numeric", "numeric", "numeric"))
	}
}

readHsqs <- function(rootname)
{
	filename <- paste(rootname, ".hsq", sep="")
	if(file.exists(filename))
	{
		a <- read.table(filename, he=T, fill=TRUE)
		return(a)
	} else {
		return(NULL)
	}
}

removeFiles <- function(rootname)
{
	a <- paste(rootname, c(".phen", ".hsq", ".indi.blp"), sep="")
	unlink(a)
}



arguments <- commandArgs(T)
jid <- as.numeric(arguments[1])
nrun <- as.numeric(arguments[2])

savefile1 <- paste("~/repo/methylation_residuals/res/results_hsq", jid, ".RData", sep="")
savefile2 <- paste("~/repo/methylation_residuals/res/results_pred", jid, ".RData", sep="")

load("~/repo/methylation_residuals/data/parameters.RData")
first <- (jid - 1) * nrun + 1
last <- min(nrow(params), jid * nrun)


print(c(first, last))

# grmroot <- "~/repo/methylation_residuals/data/alspac_hm3_GROUP_unrelated_m"
# allphenroot <- "~/repo/methylation_residuals/data/ALN.tost.betas.TIMEPOINT.OCT.ranktransformed.txt"
# idroot <- "~/repo/methylation_residuals/data/GROUP.id"
# phenroot <- "~/repo/methylation_residuals/results/TIMEPOINT/CPG.phen"
# outroot <- "~/repo/methylation_residuals/results/TIMEPOINT/CPG"

preds <- list()
hsqs <- list()
nom <- list()
j <- 1
for(i in first:last)
{
	cat(i, "\n")
	outfile <- paste("~/repo/methylation_residuals/res/", params$timepoint[i], "/", params$cpg[i], sep="")
	if(!file.exists(paste(outfile, ".hsq", sep="")))
	{
		grmfile <- paste("~/repo/methylation_residuals/data/", params$timepoint[i], sep="")
		allphenfile <- paste("~/repo/methylation_residuals/data/", params$timepoint[i], ".phen", sep="")
		idfile <- paste("~/repo/methylation_residuals/data/", params$timepoint[i], ".id", sep="")
		phenfile <- paste("~/repo/methylation_residuals/res/", params$timepoint[i], "/", params$cpg[i], ".phen", sep="")
		makePhen(allphenfile, idfile, params$index[i], phenfile)
		runGcta(grmfile, phenfile, outfile)
		nom[[j]] <- c(params$cpg[i], params$timepoint[i])
		preds[[j]] <- readPreds(outfile)
		hsqs[[j]] <- readHsqs(outfile)
		j <- j + 1
		removeFiles(outfile)
	}
}

save(nom, hsqs, file=savefile1)
save(nom, preds, file=savefile2)
