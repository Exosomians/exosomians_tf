# In this file, we will read the refined count files,
#   label the data, i.e. which sample/sequence is exported via exosome
#   Then, check the width distribution of the RNA regions
#   Then, filter the label and design matrix based on width distributions

#   Input: Raw data files (smRNA_counts) + primary design matrix(only a few features)
#   Output: Design matrix with added Labels


source('Codes/Functions.R')
Initialize()


#### Reading raw data from file ####
COUNTS_DATA_DIR = 'Data/oldDataRefined/Counts_and_Beds/'
filesPath = list.files(COUNTS_DATA_DIR,pattern = '*.txt' ,full.names = T)
countsFiles = lapply(filesPath, read.delim, header=F)
names(countsFiles) <- substr(list.files(COUNTS_DATA_DIR,pattern = '*.txt' ),1,2)

designMat = read.csv('Data/oldDataRefined/DesignMatrices/1_IC_PrimaryDesignMat.csv')

countsFiles <- lapply(countsFiles, function(aFile) {
  aFile <- subset(aFile, select=c(V1, V2, V3, V4, V5, V6, V7, V8, V9, V10))
  colnames(aFile) = c('chr', 'start', 'end', 'regionName', 'score', 'strand', 'coverage', 'bases', 'length', 'fraction')
  aFile$id <- paste0(aFile$chr, '_', aFile$start, '_', aFile$end, aFile$strand)
  aFile <- subset(aFile, id %in% designMat$id)
})

## checjing if the regions are in order
lapply(countsFiles, function(x) sum(x$id != designMat$id))
designMat$ic = countsFiles[['IC']]$coverage
designMat$ev = countsFiles[['EV']]$coverage


## defining labels
designMat$label <- ifelse(designMat$ev>quantile(designMat$ev, 0.75), 'YES', 'NO')
View(head(designMat,10))
write.csv(designMat, 'Data/oldDataRefined/DesignMatrices/2_IC_PrimaryDesignMat_label.csv', row.names = F)

## merging TCGA, serum and Hani designMat into one designMat
ICdesignMat = read.csv('Data/oldDataRefined/DesignMatrices/2_IC_PrimaryDesignMat_label.csv')
ICdesignMat$ic <- NULL
ICdesignMat$ev <- NULL
ICdesignMat$db <- 'H'
serumDesignMat = read.csv('Data/oldDataRefined/DesignMatrices/2_serum_PrimaryDesignMat_label.csv')
serumDesignMat$db <- 'S'
ICdesignMat <- rbind(ICdesignMat, serumDesignMat)
TCGAdesignMat = read.csv('Data/oldDataRefined/DesignMatrices/2_tcga_PrimaryDesignMat_label.csv')
TCGAdesignMat$db <- 'T'
ICdesignMat <- rbind(ICdesignMat, TCGAdesignMat)
write.csv(ICdesignMat, 'Data/oldDataRefined/DesignMatrices/2_all_PrimaryDesignMat_label_db.csv', row.names = F)


## Some plots
table(designMat$label)/sum(table(designMat$label))
label_p=ggplot(data.frame(table(designMat$label)), aes(x=Var1,y=Freq,color='black'))+
    geom_bar(stat = 'identity',color="dark blue", fill="cadetblue2",width=0.4)+xlab('label')+theme_bw()

pl1=ggplot(designMat, aes(x=label,y=length))+geom_boxplot(aes(fill=label))+theme_bw()+ggtitle('length distribution(refined data)')
pl2=ggplot(designMat, aes(x=label,y=length))+geom_violin(aes(fill=label))+theme_bw()+ggtitle('length distribution(refined data)')
pl3=ggplot(designMat, aes(x=length, color=label))+geom_density()+theme_bw()+ggtitle('length distribution(refined data)')
grid.arrange(pl1, pl2, pl3, ncol=3,nrow=1)

base_p1=DrawFeatureDistribution_LabelBased(subset(designMat, select=c(a,c,u,g,label)))
pheatmap(cor(subset(designMat,select=c(a,c,u,g))))
grid.arrange(label_p,base_p1, ncol=2,nrow=1)
