library(Biostrings)

merDesignMat = read.csv('Data/MergedDesignMatLabel_SecondStruct_LenFilter.csv', stringsAsFactors = FALSE)
length(merDesignMat[,merDesignMat$id != merDesignMat$X])
# remove kmers because len is 0
merDesignMat = merDesignMat[,c('id', 'chr', 'seq', 'length', 'strand', 'ic', 'ev', 'label')]
dim(merDesignMat)
merDesignMat = merDesignMat[merDesignMat$length > 30 & merDesignMat$length < 41 & merDesignMat$strand == '+',]
dim(merDesignMat)
yesMerDesignMat = merDesignMat[merDesignMat$label == 'YES', ]
dim(yesMerDesignMat)
noMerDesignMat = merDesignMat[merDesignMat$label == 'NO', ]
dim(noMerDesignMat)

noMerDesignMat = noMerDesignMat[sample(nrow(noMerDesignMat), size = 60, replace = FALSE),]
yesMerDesignMat = yesMerDesignMat[sample(nrow(yesMerDesignMat), size = 60, replace = FALSE),]
merDesignMat = rbind(yesMerDesignMat, noMerDesignMat)
# shuffle rows
merDesignMat = merDesignMat[sample(nrow(merDesignMat)),]
View(merDesignMat)

# hierarchical clustering
# test distance
# stringDist(c('blflflflflallflklkldlfl', 'adbffafkkdfjjfefdf'), method = 'levenshtein', ignoreCase = TRUE, diag = FALSE, upper = FALSE, substitutionMatrix=NULL)
distanceObject = stringDist(merDesignMat$seq, method = 'levenshtein', ignoreCase = TRUE, diag = FALSE, upper = FALSE)

hclust_single <- hclust(distanceObject, method = 'single')
plot(hclust_single, labels=merDesignMat$label)

hclust_wardD <- hclust(distanceObject, method = 'ward.D')
plot(hclust_wardD, labels=merDesignMat$label)

hclust_wardD2 <- hclust(distanceObject, method = 'ward.D2')
plot(hclust_wardD2, labels=merDesignMat$label)

hclust_complete <- hclust(distanceObject, method = 'complete')
plot(hclust_complete, labels=merDesignMat$label)

hclust_centroid <- hclust(distanceObject, method = 'centroid')
plot(hclust_centroid, labels=merDesignMat$label)

hclust_avg <- hclust(distanceObject, method = 'average')
plot(hclust_avg, labels=merDesignMat$label)

hclust_mcquitty <- hclust(distanceObject, method = 'mcquitty')
plot(hclust_mcquitty, labels=merDesignMat$label)

hclust_median <- hclust(distanceObject, method = 'median')
plot(hclust_median, labels=merDesignMat$label)


cut_avg <- cutree(hclust_avg, k = 10)
cut_avg
merDesignMat$hclust = cut_avg
plot(cut_avg)
notPredMat = merDesignMat[merDesignMat$hclust <= 2,]
dim(notPredMat)
predYesMat = merDesignMat[merDesignMat$hclust <= 311,]
predYesMat = predYesMat[!predYesMat$id %in% notPredMat$id,]
dim(predYesMat)
predNoMat = merDesignMat[merDesignMat$hclust > 311,]
predNoMat = predNoMat[!predNoMat$id %in% notPredMat$id,]
dim(predNoMat)

