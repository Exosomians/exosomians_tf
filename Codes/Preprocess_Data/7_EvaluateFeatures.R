source('Codes/Functions.R')
Initialize()

designMat <- read.csv('Data/oldDataRefined/DesignMatrices/7_DesignMat_SS_Kmer_DB_ForgiAnnot_Label.csv')
DeepBind_ids <- read.csv('Data/oldDataRefined/deepbind/deepBind_dictionary.csv',stringsAsFactors = F)
DeepBind_colnames <- colnames(designMat)[273:889]
DeepBind_df <- DeepBind_ids[DeepBind_ids$id %in% DeepBind_colnames,]
colnames(designMat)[273:889] <- DeepBind_df$protein_name


.Find_imp_feature <- function(featureToCheck, data ,minq, maxq){
  YES_data <- featureToCheck[data$label=='YES']
  NO_data <- featureToCheck[data$label=='NO']
  x = ifelse( quantile(YES_data, minq)>quantile(NO_data, maxq) | 
                quantile(NO_data, minq)>quantile(YES_data, maxq),T, F)
  return(x)
} 


isNumeric <- unlist(lapply(designMat, function(x) ifelse(class(x)%in% c('integer','numeric'), T, F)))


## raw data
feature_index <- sapply(1:ncol(designMat[,isNumeric]), 
              function(i) .Find_imp_feature(designMat[,isNumeric][,i], designMat ,0.35, 0.63))

raw_features <- colnames(designMat[,isNumeric])[which(feature_index)]


## normal data
designMat_normal <- designMat
designMat_normal[,isNumeric] <- designMat_normal[,isNumeric]/designMat_normal$length

feature_index <- sapply(1:ncol(designMat_normal[,isNumeric]), 
                        function(i) .Find_imp_feature(designMat_normal[,isNumeric][,i],designMat_normal ,0.35, 0.63))

normal_features <- colnames(designMat_normal[,isNumeric])[which(feature_index)]


## combining features
features <- unique(c(raw_features, normal_features))
features <- features[!features %in% c('ev','ic')]



.VisualFeature <- function(feature, data, dataType){
  
  tab <- data[,c(feature, 'label')]
  colnames(tab)[1] <- 'feature'
  p1=ggplot(tab, aes(x=label, y=feature))+geom_boxplot(aes(fill=label))+
    scale_fill_brewer(palette="Accent")+theme_bw()+ggtitle(paste0(feature,'_distribution ',dataType))
  p2=ggplot(tab, aes(x=feature))+geom_density(aes(fill=label),alpha=0.6,color='black')+
    scale_fill_brewer(palette="Accent")+theme_bw()
  
  grid.arrange(p1,p2,ncol=2) 
}


pdf('plots/oldDataRefined/check_features.pdf', width = 11,height = 6)
sapply(1:length(features), function(i){
  .VisualFeature(features[i], designMat, 'raw')
  .VisualFeature(features[i], designMat_normal, 'normal')
}, simplify = F) 
dev.off()


scoreTable <- read.csv('Data/oldDataRefined/featureSelection/featureScoreTable.csv')
selectedFeatures <- subset(scoreTable, sumScore>7)$X

Ds <- selectedFeatures[grep('D0*',selectedFeatures)] 
selectedFeatures[grep('D0*',selectedFeatures)] <- DeepBind_ids$protein_name[DeepBind_ids$id %in% Ds]

pdf('plots/oldDataRefined/selectedfeaturesDistrib.pdf', width = 11,height = 6)
sapply(1:length(selectedFeatures), function(i){
  .VisualFeature(selectedFeatures[i], designMat, 'raw')
  .VisualFeature(selectedFeatures[i], designMat_normal, 'normal')
}, simplify = F) 
dev.off()


### checking a good homer feature with deepbind data
pdf('YBX1.pdf')
.VisualFeature(YBX1_id[1], designMat, 'raw')
.VisualFeature(YBX1_id[1], designMat_normal, 'norm')
.VisualFeature(YBX1_id[2], designMat, 'raw')
.VisualFeature(YBX1_id[2], designMat_normal, 'norm')
dev.off()
