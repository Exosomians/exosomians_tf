library(Biostrings)
library(proxy)

label.mat = read.csv('Data/TertiaryLabelsMat.csv')
dist.mat = pairwiseAlignment(design.mat$seq, design.mat$seq[1], scoreOnly=TRUE, type='local')
single.dist = pairwiseAlignment(merDesignMat$seq[1], merDesignMat$seq[10], scoreOnly=FALSE, type='local')
single.dist


fn <- function(x, y) { return(x-y); }
dist.mat = proxy::dist(c("w", "a", "5", "3", "g", "5", "3"), method = fn)

# string sizes diff in clustering

design.mat = design.mat[design.mat$id %in% label.mat$X,]
design.mat = design.mat[1:100,]
design.mat = merge(x=design.mat, y=label.mat, by.x="id", by.y="X", all.x=TRUE, all.y=FALSE)
View(design.mat)
dist.mat <- matrix(nrow = 100, ncol = 100)
dist.mat[, 2] <- runif(100, min = 0, max = 60)
length(dist.mat)
dim(dist.mat)
typeof(dist.mat)

for (i in 1:100) {
  print(i)
  dist.mat[i:100, i] <- pairwiseAlignment(design.mat$seq[i:100], design.mat$seq[i], scoreOnly=TRUE, type='local')
}
for (i in 1:100) {
  print(i)
  dist.mat[i, i:100] <- dist.mat[i:100, i]
}

colnames(dist.mat) <- design.mat$id
rownames(dist.mat) <- design.mat$id
write.csv(dist.mat, file = 'Data/distanceMatrix.csv', row.names = TRUE)

max(dist.mat)
max.simile = max(dist.mat)
min(dist.mat)

for (i in 1:100) {
  for (j in 1:100) {
    if (i == j) {
      dist.mat[i, j] = 0
    } else {
      dist.mat[i, j] <- max.simile - dist.mat[i, j]
    }
  }
}

max(dist.mat)
min(dist.mat)

dist.mat = as.dist(dist.mat)

max(dist.mat)
min(dist.mat)
? hclust
hclust_avg <- hclust(dist.mat, method = 'centroid')
plot(hclust_avg)
cut_avg <- cutree(hclust_avg, k = 10)
design.mat$hclust = cut_avg
plot(cut_avg)
View(design.mat)


# DistanceMatrix saved on server

library(Biostrings)
library(ape)

? pairwiseAlignment

design.mat = read.csv('Data/PrimaryDesignMat.csv')
label.mat = read.csv('Data/TertiaryLabelsMat.csv')
design.mat = design.mat[ design.mat$id %in% label.mat$X,]

# dist.mat = dist.dna(design.mat$seq, as.matrix = TRUE, model='raw')

dist.mat <- matrix(nrow = 27858, ncol = 27858)
tmp = pairwiseAlignment(design.mat$seq[5:27858], design.mat$seq[5], scoreOnly=TRUE, type='local')

for (i in 1:27858) {
  print(i)
  dist.mat[i:27858, i] <- pairwiseAlignment(design.mat$seq[i:27858], design.mat$seq[i], scoreOnly=TRUE, type='local')
}

for (i in 1:27858) {
  print(i)
  dist.mat[i, i:27858] <- dist.mat[i:27858, i]
}

colnames(dist.mat) <- design.mat$id
rownames(dist.mat) <- design.mat$id
write.csv(dist.mat, file = 'Data/distanceMatrix.csv', row.names = TRUE)




