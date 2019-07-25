## checking if MBNL1 is a promising RBP in exRNA sorting machinery

source('Codes/Functions.R')
Initialize()

TotalMatrixWithStruct <- read.csv('Data/MergedDesignMatLabel_SecondStruct_LenFilter_forgi_PlusAnnot.csv', stringsAsFactors = F)
highScore_features <- read.csv('Data/featureSelection/designMat_selectedFeatures.csv')
dim(TotalMatrixWithStruct)
dim(highScore_features)


#### MBNL1 distribution
MBNL1_unmelt <- subset(highScore_features,select=c(D00120.001, label))
MBNL1 <- melt(MBNL1_unmelt)
head(MBNL1)

p1=ggplot(MBNL1, aes(y=value,x=label))+geom_boxplot(aes(fill=label))+
  theme_bw()+ggtitle('raw MBNL1 data')+ylab('deepbind score')
p2=ggplot(MBNL1, aes(x=value,color=label))+geom_density()+
  theme_bw()+ggtitle('raw MBNL1 data')+xlab('deepbind score')

MBNL1_norm <- MBNL1
MBNL1_norm$value <- MBNL1$value /highScore_features$length
p3=ggplot(MBNL1_norm, aes(y=value,x=label))+geom_boxplot(aes(fill=label))+
  theme_bw()+ggtitle('normalized MBNL1 data')+ylab('deepbind score/length')
p4=ggplot(MBNL1_norm, aes(x=value,color=label))+geom_density()+
  theme_bw()+ggtitle('normalized MBNL1 data')+xlab('deepbind score/length')

pdf('plots/MBNL1_distribution.pdf',width=11,height = 11)
gridExtra::grid.arrange(p1,p2,p3,p4,nrow=2,ncol=2)
dev.off()



### find the ids for each region
MBNL1_unmelt <- cbind(MBNL1_unmelt, TotalMatrixWithStruct[,c('id','ic','ev')]) 
head(MBNL1_unmelt)

### import the bed file to find the genomic positions for each id
bed <- read.table('Data/Annotation/combined.cluster.scoref10.bed')
colnames(bed) <- c('chr', 'start', 'end', 'id', 'score','strand')
bed <- subset(bed, select=-c(score))

MBNL1_regions <- merge(MBNL1_unmelt, bed, all.x=T, all.y=F, by.x='id',by.y='id')
MBNL1_regions <- subset(MBNL1_regions, select=c('chr','start','end','strand','id','D00120.001','ic','ev','label'))
head(MBNL1_regions,10)
write.table(MBNL1_regions,'Data/bedToCheck.bed',quote = F,row.names = F,col.names = F,sep = '\t')


MBNL1_regions <- read.delim('Data/MBNL1/bedToCheck.bed', stringsAsFactors = F,header = F)
colnames(MBNL1_regions) <- c('chr','start','end','strand','id','D00120.001','ic','ev','label')
head(MBNL1_regions)



## import the liftOvered(Hg19-to-Hg38) clip data and preprocess it
## cleaned bed files will be given to bedtools to find the intersection with our own regions

DATA_DIRECTORY = 'Data/MBNL1/liftOver/liftOvered_beds'
clip_data <- list.files(DATA_DIRECTORY,full.names = T)
clip <- lapply(clip_data, fread)
names(clip) <- c('B', 'C', 'peaks')
clip <- lapply(clip, function(x) subset(x, select=-V4))
clip <- lapply(clip, function(x) {colnames(x) = c('chr','start','end', 'count','strand');x})
clip <- lapply(clip, data.frame)

lapply(clip, head)
sapply(1:length(clip), function(i)
  write.table(clip[[i]],paste0('Data/MBNL1_',names(clip)[i],'_chip.bed'),
              quote = F,row.names = F,col.names = F,sep = '\t'))


## visualizing clip samples for comparison
samples <- rbind(data.frame(count=clip[['B']]$count,sample='B'),
      data.frame(count=clip[['C']]$count,sample='C'),
      data.frame(count=clip[['peaks']]$count,sample='peaks'))

head(samples)
ggplot(samples, aes(x=count))+geom_histogram(aes(fill=sample),color='black',bins=30,alpha=0.6)+
  scale_x_log10()+theme_bw()+ggtitle('Compare Clip Samples')+xlab('Read-count(log10)')+ylab('Frequency')





#### finding intersection and visualization

## How to find intersection with bedtools 
# bedtools intersect -wao -f 0.2 -r -a bedToCheck.bed  -b  MBNL1_C_chip.bed  > bedfile_intersect_0.2_C_clip.bed 


clip_column_names <- colnames(clip[['B']])
MBNL1_bed_colnames <- colnames(MBNL1_regions)

.preprocessIntersection <- function(intersect_file){
  colnames(intersect_file) <- c(MBNL1_bed_colnames,paste0(clip_column_names,'_clip'),'intersect')
  intersect_file$count_clip[intersect_file$count_clip=='.' ] = 0
  intersect_file$count_clip <- as.numeric(intersect_file$count_clip)
  return(intersect_file)
}


.visualizeIntersect <- function(intersect, Title){
  p1=ggplot(subset(intersect,select=c(count_clip,label)), aes(x=label,y=log10(count_clip+1)))+geom_boxplot(aes(fill=label))+theme_bw()+ggtitle(Title)
  p2=ggplot(subset(intersect,select=c(count_clip,label)), aes(color=label,x=log10(count_clip+1)))+geom_density()+theme_bw()+ggtitle(Title)
  p3=ggplot(subset(intersect,select=c(count_clip,label)), aes(x=label,y=log2(count_clip+1)))+geom_violin(aes(fill=label))+theme_bw()+ggtitle(Title)
  p4=ggplot(intersect,aes(x=D00120.001,y=log2(count_clip+1),color=label))+geom_point()+theme_bw()+ggtitle(Title)
  p5=ggplot(intersect,aes(x=log2(ev+1),y=count_clip,color=D00120.001))+geom_point()+theme_bw()+scale_color_gradient(low="blue", high="red")+ggtitle(Title)
  print(grid.arrange(p1,p2,p3, nrow=1,ncol=3))
  print(grid.arrange(p4,p5, nrow=1,ncol=2))
}



INTERSECTIONS_DIRECTORY <- 'Data/MBNL1/intersection'
intersect_files <- list.files(INTERSECTIONS_DIRECTORY, full.names = T)
intersect_files <- lapply(intersect_files, read.table)
intersect_files <- lapply(intersect_files, .preprocessIntersection)
names(intersect_files) <- c('B', 'C', 'peaks')
lapply(intersect_files, head)


intersect_0.2th_titles <- c('B sample(cutoff=0.2)', 'C sample(cutoff=0.2)')
pdf('plots/MBNLP1_clip_0.2.pdf', width = 12, height=6)
sapply(1:2,
  function(i){
    .visualizeIntersect(intersect_files[[i]], names(intersect_files)[i])
  }, simplify = F)
dev.off()




############# called peaks analysis
.cleanTable <- function(tab){
  tab = melt(tab)
  tab$count_clip[tab$count_clip == 1] <- 'Peak'
  tab$count_clip[tab$count_clip == 0] <- 'not-Peak'
  return(tab)
}

intersect_peaks <- intersect_files[['peaks']]
peaks_tab_init <- table(subset(intersect,select=c(count_clip,label)))

peaks_tab <- .cleanTable(peaks_tab_init)
p1=ggplot(peaks_tab, aes(x=count_clip,y=value,fill=label))+geom_bar(stat = 'identity',color='black')+
  theme_bw()+xlab('')+ylab('count')+scale_fill_manual(values=c('#999999','#E69F00'))+ggtitle('raw data')

## scaling by peak
pScaled_tab <- peaks_tab_init/rowSums(peaks_tab_init)
pScaled_tab <- .cleanTable(pScaled_tab)
p2=ggplot(pScaled_tab, aes(x=count_clip,y=value,fill=label))+geom_bar(stat = 'identity',color='black')+theme_bw()+
  scale_fill_manual(values=c('#999999','#E69F00'))+xlab('')+ylab('Freq')+ggtitle('scaled-by-peaks')


p3=ggplot(peaks_tab, aes(x=label,y=value,fill=count_clip))+geom_bar(stat = 'identity',color='black')+
  theme_bw()+scale_fill_brewer(palette="Blues")

lScaled_tab <- t(peaks_tab_init)/colSums(peaks_tab_init)
lScaled_tab <- .cleanTable(lScaled_tab)
p4=ggplot(lScaled_tab, aes(x=label,y=value,fill=count_clip))+geom_bar(stat = 'identity',color='black')+
  theme_bw()+scale_fill_brewer(palette="Blues")+xlab('')+ylab('Freq')+ggtitle('scaled-by-labels')

pdf('Data/MBNL1/called_peaks_analysis.pdf')
grid.arrange(p1, p2, nrow=1,ncol=2)
grid.arrange(p3, p4, nrow=1,ncol=2)
dev.off()




# bedtools intersect -wao -a bedToCheck.bed  -b  MBNL1_chip.bed  > bedfile_intersect_MBNL1chip_noFilter.bed
# bedtools intersect -wao -f 0.4 -r -a bedToCheck.bed  -b  MBNL1_chip.bed  > bedfile_intersect_MBNL1chip_0.4Filter_all.bed 


### compare diffrent filterings
x = rbind(data.frame(count=intersect_noFilter$count_chip,filter_type='NO Filter') , 
          data.frame(count=intersect_filter$count_chip,filter_type='Filter=0.4') , 
          data.frame(count=intersect$count_chip,filter_type='Filter=0.2') )
ggplot(x ,aes(x=log10(count+1)))+geom_density(aes(fill=filter_type),color='black',alpha=0.5)+
  theme_bw()+ggtitle('filtering effect')


