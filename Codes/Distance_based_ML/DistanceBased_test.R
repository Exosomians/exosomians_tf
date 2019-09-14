library(Biostrings)
library(ape)

? pairwiseAlignment
? dist.dna

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

