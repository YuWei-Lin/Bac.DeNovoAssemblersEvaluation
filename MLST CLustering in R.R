### Dissimilarity Clustering Computation
# Import MLST summary table
STable <- read.csv(file.choose(), header = T, stringsAsFactors = F)
# Take out columns that have no allele call.
STable <- STable[ ,-(10:11)]
STable <- STable[ ,-2]
# Rename matrix
rownames(STable) <- STable$StrainNum.
STable <- STable[ ,-1]
# Calculate percentage similarity for each pair of isolates
SimPer <- matrix(NA, nrow = 8, ncol = 8)
for (i in 1:nrow(STable)) {
  for (j in 1:nrow(STable)) {
    SimPer[i,j] <- 100-(length(intersect(STable[i,], STable[j,]))/7*100)
  }
}
# Rename percentage similarity matrix 
rownames(SimPer) <- rownames(STable)
colnames(SimPer) <- rownames(STable)
# Convert to relative distances
SimPer <- as.dist(SimPer)
# Perform hierarchical clustering
CLU <- hclust(SimPer, method = "average")
CLU <- as.dendrogram(CLU)
# Define nodePar
nodePar <- list(lab.cex = 0.8, pch = c(NA, 21), cex = 1.5, col = "red", bg = "red")
# Customized plot; remove labels
plot(CLU, main = "Streptococcus dysgalactiae equisimilis ST", xlim = c(0, 100), xlab = "Percentage Disimilarity", nodePar = nodePar, horiz = TRUE)
