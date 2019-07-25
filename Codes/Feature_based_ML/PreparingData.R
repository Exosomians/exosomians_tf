designMat = read.csv('Data/PrimaryDesignMat.csv')
labelsMat = read.csv('Data/SecondaryLabelsMat.csv')

finalMat = merge(labelsMat[, c('id', 'label')], designMat, all.x = T)

write.csv(finalMat, file = 'Data/PrimaryData.csv', quote = F, na = 'NA', row.names = F, col.names = T)

# install.packages('caret')
library(caret)

dataFolds = createFolds(y = finalMat$label, k = 10, list = T)

dataFolded = lapply(dataFolds, function(aFold) {
  testData = finalMat[aFold, ]
  trainData = finalMat[setdiff(seq(nrow(finalMat)), aFold), ]
  list(Train = trainData, Test = testData)
})

dir.create('Data/SecondaryDataFolded')

lapply(names(dataFolded), function(aData)
  {
  write.csv(dataFolded[[aData]]$Train, file = sprintf('Data/SecondaryDataFolded/%s_Train.csv', aData) , quote = F, na = 'NA', row.names = F, col.names = T)
  write.csv(dataFolded[[aData]]$Test, file = sprintf('Data/SecondaryDataFolded/%s_Test.csv', aData) , quote = F, na = 'NA', row.names = F, col.names = T)
})


#### Balancing class labels ####
finalMat = read.csv('Data/PrimaryData.csv')
BALANCING_WEIGHT = as.integer(floor(table(finalMat$label)[1]/table(finalMat$label)[2]))


finalMat$balancingWeight = 1
finalMat$balancingWeight[finalMat$label=='Yes'] = BALANCING_WEIGHT

finalMatBalanced = finalMat[rep(rownames(finalMat), finalMat$balancingWeight), ]
# rownames(finalMatBalanced) = NULL
finalMatBalanced = finalMatBalanced[sample(rownames(finalMatBalanced)), ]


write.csv(finalMatBalanced, file = 'Data/PrimaryDataBalanced.csv', quote = F, na = 'NA', row.names = F, col.names = T)


