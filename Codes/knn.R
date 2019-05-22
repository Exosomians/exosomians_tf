library(Biostrings)

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

# my own KNN!
test_knn_model <- function(data, testIndex, k) {
  print('start...')
  print('total number of the rows:')
  print(nrow(data))
  
  trainIndex = 1:nrow(data)
  trainIndex = trainIndex[-testIndex]
  print('number of training rows:')
  print(length(trainIndex))
  
  distanceObject = stringDist(data$seq, method = 'levenshtein', ignoreCase = TRUE, diag = FALSE, upper = FALSE)
  distMat <- as.matrix(distanceObject)
  print('distance matrix computed.')
  testData <- data[testIndex,]
  print('number of test rows:')
  print(nrow(testData))
  
  testData$testLablel = lapply(testIndex, function(x) {
    distances <- distMat[x, -testIndex]
    indexes = c()
    i <- 0
    repeat {
      index = which.min(distances)
      indexes = c(indexes, index)
      distances <- distances[-index]
      i = i+1
      if (i >= k) {
        break
      }
    }
    minData = data[indexes,]
    yesLabels = nrow(minData[minData$label=='YES',])
    noLabels = nrow(minData[minData$label=='NO',])
    if (yesLabels >= noLabels) {
      return('YES')
    } else {
      return('NO')
    }
  })
  print('---------------')
  truePositive = nrow(testData[testData$label == 'YES' & testData$testLablel == 'YES',])
  print('truePositive:')
  print(truePositive)
  trueNegative = nrow(testData[testData$label == 'NO' & testData$testLablel == 'NO',])
  print('trueNegative:')
  print(trueNegative)
  falsePositive = nrow(testData[testData$label == 'NO' & testData$testLablel == 'YES',])
  print('falsePositive:')
  print(falsePositive)
  falseNegative = nrow(testData[testData$label == 'YES' & testData$testLablel == 'NO',])
  print('falseNegative:')
  print(falseNegative)
  sensitivity = truePositive / (truePositive + falseNegative)
  print('sensitivity:')
  print(sensitivity)
  specificity = trueNegative / (trueNegative + falsePositive)
  print('specificity:')
  print(specificity)
  print('---------------')
}
dim(merDesignMat)
View(merDesignMat)
testIndex = sample(1:nrow(merDesignMat), size=20, replace = FALSE)
test_knn_model(data=merDesignMat, testIndex = testIndex, k=4)

