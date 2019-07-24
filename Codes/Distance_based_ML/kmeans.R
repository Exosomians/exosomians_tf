library(Biostrings)
library(class)
library(cluster)

merDesignMat = read.csv('Data/MergedDesignMatLabel_LenFilter.csv', stringsAsFactors = FALSE)
length(merDesignMat[,merDesignMat$id != merDesignMat$X])
# remove kmers because len is 0
merDesignMat = merDesignMat[,c('id', 'chr', 'seq', 'length', 'strand', 'ic', 'ev', 'label')]
dim(merDesignMat)
merDesignMat = merDesignMat[merDesignMat$length > 20 & merDesignMat$length < 30 & merDesignMat$strand == '+',]
dim(merDesignMat)
yesMerDesignMat = merDesignMat[merDesignMat$label == 'YES', ]
dim(yesMerDesignMat)
noMerDesignMat = merDesignMat[merDesignMat$label == 'NO', ]
dim(noMerDesignMat)

noMerDesignMat = noMerDesignMat[sample(nrow(noMerDesignMat), size = 30, replace = FALSE),]
yesMerDesignMat = yesMerDesignMat[sample(nrow(yesMerDesignMat), size = 30, replace = FALSE),]
merDesignMat = rbind(yesMerDesignMat, noMerDesignMat)
# shuffle rows
merDesignMat = merDesignMat[sample(nrow(merDesignMat)),]
View(merDesignMat)

# compute distance matrix
distanceObject = stringDist(merDesignMat$seq, method = 'levenshtein', ignoreCase = TRUE, diag = FALSE, upper = FALSE)

# Compute K-means
View(merDesignMat)
distMat = as.matrix(distanceObject)
k <- 4
kmRes <- kmeans(distMat, k, nstart = 25)
attributes(kmRes)
View(kmRes$centers)

merDesignMat$kmclust = kmRes$cluster
j <- 0
repeat {
  j <- j+1
  kRows = merDesignMat[merDesignMat$kmclust == j,]
  yes = kRows[kRows$label == 'YES',]
  no = kRows[kRows$label == 'NO',]
  tot = nrow(kRows)
  yesn = nrow(yes)
  non = nrow(no)
  print('k:')
  print(k)
  print('tot:')
  print(tot)
  print('yes:')
  print(yesn/tot)
  print('no:')
  print(non/tot)
  print('_________________________________')
  if (j >= k) {
    break();
  }
}

# Visualize
fviz_cluster(kmRes, distMat)

# my own K-means!
test_k_means_model <- function(designMat, distMat, k, maxIterations) {
  # random choose k centers
  centersIndex <- sample(nrow(designMat), size = k, replace = FALSE)
  
  # assign initial center clusters
  j <- 0
  repeat {
    j <- j+1
    designMat[centersIndex[j],]$kmclust = j
    if (j >= k) {
      break();
    }
  }
  
  # assign initial clusters
  j <- 0
  repeat {
    j <- j+1
    distances <- distMat[j, centersIndex]
    minIndex <- which.min(distances)
    clust <- designMat[centersIndex[minIndex],]$kmclust
    if (clust != designMat[j,]$kmclust) {
      designMat[j,]$kmclust <- clust
    }
    if (j >= nrow(designMat)) {
      break()
    }
  }
  
  # iterate on new clusters and centroids
  iter <- 0
  repeat {
    iter <- iter+1
    # find new centers (centroid update)
    j <- 0
    repeat {
      j <- j+1
      points <- which(designMat[,]$kmclust == designMat[centersIndex[j],]$kmclust)
      distances <- distMat[points, points]
      sums <- rowSums(distances)
      centersIndex[j] = points[which.min(sums)]
      if (j >= length(centersIndex)) {
        break()
      }
    }
    # assign new clusters
    j <- 0
    isChanged <- FALSE
    repeat {
      j <- j+1
      distances <- distMat[j, centersIndex]
      minIndex <- which.min(distances)
      clust <- designMat[centersIndex[minIndex],]$kmclust
      if (clust != designMat[j,]$kmclust) {
        designMat[j,]$kmclust <- clust
        isChanged <- TRUE
      }
      if (j >= nrow(designMat)) {
        break()
      }
    }
    if (iter >= maxIterations || isChanged == FALSE) {
      break();
    }
  }
  # printing the results
  j <- 0
  repeat {
    j <- j+1
    kRows = designMat[designMat$kmclust == j,]
    yes = kRows[kRows$label == 'YES',]
    no = kRows[kRows$label == 'NO',]
    tot = nrow(kRows)
    yesn = nrow(yes)
    non = nrow(no)
    print('k:')
    print(j)
    print('tot:')
    print(tot)
    print('yes:')
    print(yesn/tot)
    print('no:')
    print(non/tot)
    print('_________________________________')
    if (j >= k) {
      break();
    }
  }
}

merDesignMat$kmclust <- 0
View(merDesignMat)
distMat = as.matrix(distanceObject)

test_k_means_model(merDesignMat, distMat, 5, 25)

