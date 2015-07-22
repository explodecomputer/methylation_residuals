tp <- c("cord", "F7", "15up", "antenatal", "FOM")
phenfile <- paste("~/repo/methylation_residuals/data/", tp, "norm.phen", sep="")
outfile <- paste("~/repo/methylation_residuals/data/", tp, "norm.lowpres.phen", sep="")
stopifnot(all(file.exists(phenfile)))

for(i in 1:length(phenfile))
{
	message(tp[i])
	a <- read.table(phenfile[i], colClass="numeric")
	write.table(round(a, 3), outfile[i], row=F, col=F, qu=F)
}
