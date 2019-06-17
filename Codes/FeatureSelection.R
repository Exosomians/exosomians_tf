source('Codes/Functions.R')
Initialize()


### loading data

TotalMatrixWithStruct <- read.csv('Data/MergedDesignMatLabel_SecondStruct_LenFilter.csv', stringsAsFactors = F)
ColumnsToDrop <- c('id','X','X.1','X.2','ic','ev','seq','DB','annotation','rnaType')
designMat <- TotalMatrixWithStruct[,!colnames(TotalMatrixWithStruct) %in% ColumnsToDrop ]



### making all the designMatrices 

seqKmer <- subset(designMat, select= (AAAA:UUUU))
secStr <- (subset(designMat, select=(A000:G111)))
colnames(secStr) <- designMat$label


### PCA on second structure features

secStr_feature_PCA_YES <- prcomp( t(secStr[designMat$label=='YES',]) , scale = TRUE)
fviz_eig(secStr_feature_PCA_YES)
fviz_pca_ind(secStr_feature_PCA_YES,
             col.ind = "cos2", # Color by the quality of representation
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE)



PCA_secStr <- autoplot(prcomp(secStr), data = designMat, colour = 'label')
PCA_seqKmer <- autoplot(prcomp(seqKmer), data = designMat, colour = 'label')


secStrDif <- ComputeDEbyLimma(designMat$label, secStr)
seqKmerDif <- ComputeDEbyLimma(designMat$label, seqKmer)


secStrDEfeatures <- rownames(subset(secStrDif, logFC > 2 | logFC < (-2)))
secStr_onlyDEincluded <- secStr[,as.character(colnames(secStr)) %in% secStrDEfeatures ]
autoplot(prcomp(secStr_onlyDEincluded, scale. = T), data = designMat, colour = 'label')







