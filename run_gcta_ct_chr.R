makePhen <- function(allphenfile, idfile, row, phenfile)
{
	cmd <- paste("head -n ", row, " ", allphenfile, " | tail -n 1 | tr ' ' '\n' > ", phenfile, ".temp", sep="")
	system(cmd)
	cmd <- paste("paste -d ' ' ", idfile, " ", phenfile, ".temp > ", phenfile, sep="")
	system(cmd)
	system(paste("rm ", phenfile, ".temp", sep=""))
}

getCisChr <- function(probeinfo, cpg)
{
	return(probeinfo$CHR[probeinfo$TargetID == cpg][1])
}

runGcta <- function(cpg, grmroot, phenfile, outfile, flags="")
{
	chr <- getCisChr(probeinfo, cpg)
	mgrmfile <- paste(grmroot, chr, ".mgrm", sep="")
	cmd <- paste("gcta64 ", flags, " --mgrm ", mgrmfile, " --reml --reml-no-lrt --pheno ", phenfile, " --reml-pred-rand --out ", outfile, sep="")
	system(cmd)
}

readPreds <- function(rootname)
{
	filename <- paste(rootname, ".indi.blp", sep="")
	if(file.exists(filename))
	{
		a <- read.table(filename, colClass=c("character", "character", "numeric", "numeric", "numeric", "numeric"))
		a$probe <- basename(rootname)
		return(a)
	}
}

readHsqs <- function(rootname)
{
	filename <- paste(rootname, ".hsq", sep="")
	if(file.exists(filename))
	{
		a <- read.table(filename, he=T, fill=TRUE)
		a$probe <- basename(rootname)
		return(a)
	} else {
		return(NULL)
	}
}

removeFiles <- function(rootname)
{
	a <- paste(rootname, c(".phen", ".hsq", ".indi.blp", ".log", ".grm.id", ".grm.bin", ".grm.N.bin", ".snplist", ".mgrm"), sep="")
	unlink(a)
}

arguments <- commandArgs(T)
jid <- as.numeric(arguments[1])
nrun <- as.numeric(arguments[2])

savefile1 <- paste("~/repo/methylation_residuals/res_ct_chr/results_ct_hsq", jid, ".RData", sep="")
savefile2 <- paste("~/repo/methylation_residuals/res_ct_chr/results_ct_pred", jid, ".RData", sep="")

if(all(file.exists(c(savefile1, savefile2)))) q()

load("~/repo/methylation_residuals/data/parameters.RData")

win <- 1000000
first <- (jid - 1) * nrun + 1
last <- min(nrow(params), jid * nrun)

print(c(first, last))
preds <- list()
hsqs <- list()
nom <- list()

j <- 1
for(i in first:last)
{
	cat("Running param", i, "\n")
	outfile <- paste("~/repo/methylation_residuals/res_ct_chr/", params$timepoint[i], "/", params$cpg[i], sep="")
	if(!file.exists(paste(outfile, ".hsq", sep="")) & params$cpg[i] %in% probeinfo$TargetID)
	{
		grmroot <- paste("~/repo/methylation_residuals/grms/", params$group[i], sep="")
		allphenfile <- paste("~/repo/methylation_residuals/data/", params$timepoint[i], "norm.phen", sep="")
		idfile <- paste("~/repo/methylation_residuals/data/", params$timepoint[i], ".id", sep="")

		cisgrmfile <- outfile
		phenfile <- paste(outfile, ".phen", sep="")

		makePhen(allphenfile, idfile, params$index[i], phenfile)
		runGcta(params$cpg[i], grmroot, phenfile, outfile)

		nom[[j]] <- c(params$cpg[i], params$timepoint[i])
		preds[[j]] <- readPreds(outfile)
		hsqs[[j]] <- readHsqs(outfile)
		j <- j + 1

		removeFiles(outfile)

	}
}

save(nom, hsqs, file=savefile1)
save(nom, preds, file=savefile2)
