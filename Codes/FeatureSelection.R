## in this script, we'll apply different feature selection 
## methods on the data to decrease the #features and eliminate  
## the redundant and uninformative ones

source('Codes/Functions.R')
set.seed(123)

Initialize()
FeatureSelectionLibs()


#### Functions 
### Visualizing data clusters for feature-set evaluation  

.checkFeatureByPCA <- function(InputSubset, designMat){
  autoplot(prcomp(InputSubset, scale. = T), data = designMat, colour = 'label')
}

.checkFeatureBytSNE <- function(InputSubset, designMat){
  tsne <- Rtsne(InputSubset, dims = 2, perplexity=30, verbose=TRUE, max_iter = 500, check_duplicates = FALSE)
  tsne_plot <- data.frame(x = tsne$Y[,1], y = tsne$Y[,2], col = designMat$label)
  ggplot(tsne_plot) + geom_point(aes(x=x, y=y, color=col))+theme_bw()+ggtitle('tSNE')
}

.checkFeatureByUMAP <- function(InputSubset, designMat){
  umap <- umap(InputSubset)
  umap_plot <- data.frame(x = umap$layout[,1], y = umap$layout[,2], col = designMat$label)
  ggplot(umap_plot) + geom_point(aes(x=x, y=y, color=col))+theme_bw()+ggtitle('UMAP')
}


### finding the features which are differentially expressed(classes: YES,NO) 
### then checking the pca plot to see if the classes are devided better with that subset of features
#### result -> wasn't much effective -> moving to traditional feature selection methods

.ComputeDEbyLimma <- function(labels, data){
  
  gr <- factor(labels)
  TransposedData <- t(data)
  
  #### limma
  data$description <- gr
  design <- model.matrix( ~description + 0, data)
  colnames(design) <- levels(gr)
  
  ## fit a linear model 
  fit <- lmFit(TransposedData, design)
  cont.matrix <- makeContrasts(YES-NO, levels = design)
  
  fit2 <- contrasts.fit(fit, cont.matrix)
  
  fit2 <- eBayes(fit2, 0.01) ## prior knowledge: fraction of genes to have significant difference in expression
  dif <- topTable(fit2, adjust='fdr', sort.by = 'B', number = Inf) 
  
  return( dif[order(as.numeric(dif$logFC), decreasing = T), ] )
}





#### loading data

TotalMatrixWithStruct <- read.csv('Data/MergedDesignMatLabel_SecondStruct_LenFilter_forgi_PlusAnnot.csv', stringsAsFactors = F)
deepbind <- read.delim('Data/deepbind_features.txt')

## data preprocessing
ColumnsToDrop <- c('id','X','X.1','X.2','ic','ev','seq','DB','annotation','rnaType','element_string','element_string_number')
designMat <- TotalMatrixWithStruct[,!colnames(TotalMatrixWithStruct) %in% ColumnsToDrop ]

## previous results(wrapper method) showed the trimer Second-struct features are not much informative
designMat <- designMat[,-c(which(colnames(designMat)=="A000"):which(colnames(designMat)=="G111"))]

designMat <- cbind(designMat, deepbind)

categorical_index = unlist(lapply(designMat, function(x) !class(x) %in% c('numeric', 'integer')))
colnames(designMat)[categorical_index]

numerical_designMat <- subset(designMat, select=-c(chr,strand))
numerical_designMat$label <- ifelse(numerical_designMat$label=='YES', 1, -1)

numerical_designMat_norm <- lapply(subset(numerical_designMat, select= -c(length, label)), 
                                   function(x) x/numerical_designMat$length)
numerical_designMat_norm <- data.frame(do.call(cbind, numerical_designMat_norm))

designMat_norm <- cbind(numerical_designMat_norm, 
                        length=designMat$length, 
                        chr=designMat$chr, 
                        strand=designMat$strand, 
                        label=designMat$label)

numerical_designMat_norm <- cbind(numerical_designMat_norm, 
                                  length=numerical_designMat$length, 
                                  label=numerical_designMat$label)


## contains all the features(numeric+categoriacl), categorical label
head(designMat)
## contains only the numeric features, numeric label
head(numerical_designMat_norm)




##########################################
### Correlation

correlationMatrix <- cor(subset(numerical_designMat_norm, select=-label))
correlation_melt <- melt(correlationMatrix)

### table of highly correlated features
x = subset(correlation_melt, Var1 != Var2 & (value>0.8 | value<(-0.8)) )


# if  (row_1= a b  row_2= b a ) then remove row_2
# we want to keep the one of the features(1st arguments)  
# such duplicates need to be removed before selection, 
# if not, both features will be removed

x$id <- paste0(x$Var1, '_', x$Var2)
x$id_rev <- paste0(x$Var2, '_', x$Var1)

uniq_correlated = c()

for(i in 1:nrow(x)){
  if( (!x[i,'id'] %in% uniq_correlated) & (!x[i,'id_rev'] %in% uniq_correlated)) 
    uniq_correlated = c(uniq_correlated, x[i,'id'])
}

x_filt <- x[x$id %in% uniq_correlated,1:3]
correlated_features_to_remove <- as.character(unique(x_filt$Var1) )


### removing redundant features
designMat <- designMat[,!colnames(designMat) %in% correlated_features_to_remove]
designMat_norm <- designMat_norm[,!colnames(designMat_norm) %in% correlated_features_to_remove]
numerical_designMat <- numerical_designMat[,!colnames(numerical_designMat) %in% correlated_features_to_remove]
numerical_designMat_norm <- numerical_designMat_norm[,!colnames(numerical_designMat_norm) %in% correlated_features_to_remove]




##########################################
#### DE approach
Dif_features <- .ComputeDEbyLimma(designMat$label, subset(numerical_designMat, select=-label))
DE_features <- rownames(subset(Dif_features, logFC > 2 | logFC < (-2)))

onlyDEincluded <- numerical_designMat[,as.character(colnames(numerical_designMat)) %in% DE_features ]
.checkFeatureByPCA(onlyDEincluded, designMat)
.checkFeatureBytSNE(onlyDEincluded, designMat)
.checkFeatureByUMAP(onlyDEincluded, designMat)

## all features together
.checkFeatureBytSNE(subset(numerical_designMat, select=-label), designMat)




##########################################
### Hypothesis testing  
# Null Hypothesis: no relationship

## t-test for numeric features
pvalues <- lapply(subset(numerical_designMat, select=-label), 
       function(x){
         t_test_res = t.test(x, numerical_designMat$label); t_test_res$p.value})

pvalues <- do.call(rbind, pvalues)
sum(pvalues>0.01)
colnames(numerical_designMat)[which((pvalues>0.01))]
### only three features where excluded -> not a discriminative approach


## chi-square for categorical features
.selectBychiSq <- function(feature, labels){
  tbl <- table(feature, labels)
  tbl <- tbl[rownames(tbl) != 'chrM', ]
  chi <- chisq.test(tbl)
  return(chi$p.value)
}

strand_pvalue <- .selectBychiSq(designMat$strand, designMat$label)
chr_pvalue <- .selectBychiSq(designMat$chr, designMat$label)

## 'chr' is more predictive than 'strand'



##########################################
### Information Value

factor_vars <- c ('chr', 'strand')  
all_iv <- data.frame(VARS=factor_vars, 
                     IV=as.numeric(length(factor_vars)), 
                     STRENGTH=as.character(length(factor_vars)), 
                     stringsAsFactors = F)  

for (factor_var in factor_vars){
  all_iv[all_iv$VARS == factor_var, "IV"] <- InformationValue::IV(X=designMat[, factor_var], Y=designMat$label)
  all_iv[all_iv$VARS == factor_var, "STRENGTH"] <- attr(InformationValue::IV(X=designMat[, factor_var], Y=designMat$label), "howgood")
}

all_iv <- all_iv[order(-all_iv$IV), ]  # sort

## 'chr', 'strand' -> IV=0 , NOT Predictive 


##########################################
### Information gain 

infor_gain <- information.gain(label~., designMat)
infor_gain$feature <- rownames(infor_gain)
infor_gain[order(infor_gain$attr_importance, decreasing = T),]
summary(infor_gain$attr_importance)
imp_features_raw <- infor_gain$feature[infor_gain$attr_importance > quantile(infor_gain$attr_importance, 0.85)]


infor_gain <- information.gain(label~., designMat_norm)
infor_gain$feature <- rownames(infor_gain)
infor_gain[order(infor_gain$attr_importance, decreasing = T),]
summary(infor_gain$attr_importance)
imp_features_norm <- infor_gain$feature[infor_gain$attr_importance > quantile(infor_gain$attr_importance, 0.85)]

intersect(imp_features_raw,imp_features_norm )




##########################################
### Generalized cross validation (GCV)

# variable importance based on GCV, 
# number of subset models the variable occurs (nsubsets) 
# and residual sum of squares (RSS).

marsModel <- earth(label ~ ., data=numerical_designMat) # build model
ev <- evimp (marsModel) # estimate variable importance
GCV_imp <- rownames(ev)

marsModel <- earth(label ~ ., data=numerical_designMat_norm) 
ev_norm <- evimp (marsModel) 
GCV_imp_norm <- rownames(ev)




########################################## time
### Forward-Selection 

# Creating sample instances
rdesc <- makeResampleDesc("Holdout")

# Setting the Sequential forward Search - "sfs" 
# for Sequential Backward Search - "sbs""
ctrl <- makeFeatSelControlSequential(method = "sfs", maxit = NA)

# Running the selection algorithm
task <- makeClassifTask(data = cbind(subset(numerical_designMat,select=-label),label=designMat$label),
                       target = 'label')

res <- selectFeatures("classif.rpart", task, rdesc, control = ctrl)
saveRDS(res, 'forwardSel.rds' )
analyzeFeatSelResult(res)  ### only selected feature:  D00120.001


data_norm <- cbind(subset(numerical_designMat_norm,select=-label),label=designMat$label)
task_norm <- makeClassifTask(data = data_norm, target = 'label')
res_norm <- selectFeatures("classif.rpart", task_norm, rdesc, control = ctrl)
saveRDS(res_norm, 'forwardSel_norm.rds' )



##########################################  time !!!!
### Random Forest

cf1 <- cforest(label ~ . , data= numerical_designMat_norm, control=cforest_unbiased(mtry=2,ntree=50)) 
varImpRF <- varimp(cf1) # var imp based on mean decrease in accuracy
saveRDS(varImpRF, 'varImpRF.rds')
varImpRF_adj <- varimp(cf1, conditional=TRUE)  # conditional=True, adjusts for correlations between predictors 
saveRDS(varImpRF_adj, 'varImpRF_adj.rds')



########################################## time !!
### Boruta 

boruta_output <- Boruta(label ~ ., data=na.omit(numerical_designMat_norm), doTrace=2)  # perform Boruta search
saveRDS(boruta_output, 'boruta_output.rds')
boruta_signif <- names(boruta_output$finalDecision[boruta_output$finalDecision %in% c("Confirmed", "Tentative")]) 
boruta_signif  # significant variables
rejected_features <- names(boruta_output$finalDecision)[boruta_output$finalDecision != 'Confirmed'] ## rejected variables


boruta_impHist  <- boruta_output$ImpHistory
boruta_impHist <- subset(boruta_impHist, select=-c(shadowMax,shadowMean, shadowMin))
boruta_mean <- sapply(1:ncol(boruta_impHist), function(i) mean(boruta_impHist[i]))

hist(boruta_means)
summary(boruta_means)

CUT_OFF = 20 
imp_index <- sapply(1:ncol(boruta_impHist), function(i) mean(boruta_imp[i])>CUT_OFF)
imp_features <- colnames(boruta_impHist)[unlist(imp_index)]
boxplot(boruta_impHist[,imp_features])




########################################## time !!!
### Relative Importance  

lmMod <- lm(label ~ . , data = numerical_designMat_norm)  # fit lm() model
relImportance <- calc.relimp(lmMod, type = "lmg", rela = TRUE)  # calculate relative importance scaled to 100 
saveRDS(relImportance, 'relativeImportance.rds')
sort(relImportance$lmg, decreasing=TRUE)  # relative importance




##########################################  time !!
### Recursive feature elimination 

inTrain <- createDataPartition(y = numerical_designMat[,'label'],
                               p = 0.75,
                               list = FALSE)
training <- numerical_designMat[ inTrain,]
testing <- numerical_designMat[-inTrain,]

training$label <- ifelse(training$label==1, 'YES', 'NO')
testing$label <- ifelse(testing$label==1, 'YES', 'NO')


# Setting the cross validation parameters
ctrl_param <- rfeControl(functions = rfFuncs,
                         method = "repeatedcv",
                         repeats = 5,
                         verbose = FALSE,
                         returnResamp = "all")

rfe_lm_profile <- rfe(training[, colnames(training)!='label'], training$label, #
                      sizes = c(1:20),
                      rfeControl = ctrl_param,
                      newdata = testing[, colnames(testing)!='label'])

saveRDS(rfe_lm_profile, 'recursiveFeatElim.rds')
predictors(rfe_lm_profile)
plot(rfe_lm_profile, type=c("g", "o"))


