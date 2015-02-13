jid <- as.numeric(commandArgs(T))
gap <- 20

first <- (jid - 1) * gap + 1
last <- min(jid * gap, 1806)
# 122

n <- length(first:last)

filename <- paste("~/repo/methylation_residuals/res/results_ct_hsq", first:last, ".RData", sep="")
print(filename)
l1 <- l2 <- l3 <- l4 <- l5 <- l6 <- l7 <- list()
for(i in 1:n)
{
	cat(i, "\n")
	load(filename[i])
	nom[sapply(hsqs, is.null)] <- NULL
	hsqs[sapply(hsqs, is.null)] <- NULL
	l1[[i]] <- sapply(hsqs, function(x) { x[5,2] })
	l2[[i]] <- sapply(hsqs, function(x) { x[5,3] })
	l3[[i]] <- sapply(hsqs, function(x) { x[6,2] })
	l4[[i]] <- sapply(hsqs, function(x) { x[6,3] })
	l5[[i]] <- sapply(hsqs, function(x) { x[8,2] })
	l6[[i]] <- sapply(hsqs, function(x) { x[1,4] })
	l7[[i]] <- do.call(rbind, nom)
	index <- l7[[i]][,1] %in% l6[[i]]
	l7[[i]] <- l7[[i]][index,]
	stopifnot(dim(l7) == length(l1))
}

save(l1, l2, l3, l4, l5, l6, l7, file=paste("~/repo/methylation_residuals/res/allres", jid, ".RData", sep=""))

