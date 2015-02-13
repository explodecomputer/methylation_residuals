jid <- as.numeric(commandArgs(T))
gap <- 20

first <- (jid - 1) * gap + 1
last <- min(jid * gap, 1806)
# 122

n <- length(first:last)

filename <- paste("~/repo/methylation_residuals/res/results_ct_hsq", first:last, ".RData", sep="")
print(filename)
l1 <- list()
l2 <- list()
l3 <- list()
l4 <- list()
for(i in 1:n)
{
	cat(i, "\n")
	load(filename[i])
	nom[sapply(nom, is.null)] <- NULL
	hsqs[sapply(hsqs, is.null)] <- NULL
	l1[[i]] <- sapply(hsqs, function(x) { x[4,2] })
	l2[[i]] <- sapply(hsqs, function(x) { x[4,3] })
	l3[[i]] <- sapply(hsqs, function(x) { x[6,2] })
	l4[[i]] <- do.call(rbind, nom)
}

save(l1, l2, l3, l4, file=paste("~/repo/methylation_residuals/res/allres", jid, ".RData", sep=""))

