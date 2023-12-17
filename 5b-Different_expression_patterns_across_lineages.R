set.seed(3112)

# For reproducibility
RNGversion("3.5.0")
palette(brewer.pal(8, "Dark2"))
data(countMatrix, package = "tradeSeq")
counts <- as.matrix(countMatrix)
rm(countMatrix)
data(crv, package = "tradeSeq")
data(celltype, package = "tradeSeq")

foldername <- "/Trajectory/Markers/"
system( paste0("mkdir -p ", foldername) )

#Discovering genes with different expression patterns across lineages
patternRes <- patternTest(sce_fitted)
write.table( patternRes,
             paste0(foldername, "Different_exp_patterns_across_lineages.tsv"),
             col.names = T,
             row.names = T,
             sep = "\t",
             quote = F
)


oPat <- order(patternRes$waldStat, decreasing = TRUE)
plotSmoothers(sce_fitted, counts, gene = rownames(patternRes)[oPat][1]) #Plot expression of gene with most variation across lineages

