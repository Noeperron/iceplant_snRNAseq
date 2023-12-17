set.seed(3112)

####Identifying empty droplets using emptyDrops####
#D8DS
limit <- 230 #All cells with less than this number of UMIs are considered as "empty"
D8DS.out <- emptyDrops(counts_D8DS, lower=limit, test.ambient=TRUE)
summary(D8DS.out$FDR <= 0.001)
sum(D8DS.out$Total <= limit & D8DS.out$Total > 0, na.rm=TRUE) #Gives number of empty droplets detected

par(mfrow=c(1, 1), mar=c(5, 4, 4, 2) + 0.1)
hist(D8DS.out$PValue[D8DS.out$Total <= limit & D8DS.out$Total > 0],
     xlab="P-value", main="", col="grey80") #Plot pvalues for droplet below limit. Adjust limit number until you get a uniform distribution.

#D8DC
limit <- 460 #All cells with less than this number of UMIs are considered as "empty"
D8DC.out <- emptyDrops(counts_D8DC, lower=limit, test.ambient=TRUE)
summary(D8DC.out$FDR <= 0.001)
sum(D8DC.out$Total <= limit & D8DC.out$Total > 0, na.rm=TRUE) #Gives number of empty droplets detected

par(mfrow=c(1, 1), mar=c(5, 4, 4, 2) + 0.1)
hist(D8DC.out$PValue[D8DC.out$Total <= limit & D8DC.out$Total > 0],
     xlab="P-value", main="", col="grey80") #Plot pvalues for droplet below limit. Adjust limit number until you get a uniform distribution.

#D8LS
limit <- 330 #All cells with less than this number of UMIs are considered as "empty"
D8LS.out <- emptyDrops(counts_D8LS, lower=limit, test.ambient=TRUE)
summary(D8LS.out$FDR <= 0.001)
sum(D8LS.out$Total <= limit & D8LS.out$Total > 0, na.rm=TRUE) #Gives number of empty droplets detected

par(mfrow=c(1, 1), mar=c(5, 4, 4, 2) + 0.1)
hist(D8LS.out$PValue[D8LS.out$Total <= limit & D8LS.out$Total > 0],
     xlab="P-value", main="", col="grey80") #Plot pvalues for droplet below limit. Adjust limit number until you get a uniform distribution.

#D8LC
limit <- 500 #All cells with less than this number of UMIs are considered as "empty"
D8LC.out <- emptyDrops(counts_D8LC, lower=limit, test.ambient=TRUE)
summary(D8LC.out$FDR <= 0.001)
sum(D8LC.out$Total <= limit & D8LC.out$Total > 0, na.rm=TRUE) #Gives number of empty droplets detected

par(mfrow=c(1, 1), mar=c(5, 4, 4, 2) + 0.1)
hist(D8LC.out$PValue[D8LC.out$Total <= limit & D8LC.out$Total > 0],
     xlab="P-value", main="", col="grey80") #Plot pvalues for droplet below limit. Adjust limit number until you get a uniform distribution.

