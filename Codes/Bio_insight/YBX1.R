## checking if YBX1 is a promising RBP in exRNA sorting machinery

source('Codes/Functions.R')
Initialize()


# bed <- read.table('Data/oldDataRefined/Counts_and_Beds/all_inrange_18l_500l_20covered_regions_annotated.srt.bed')

## refine the hg19 bed file for lifting Over 
YBX1_Counts = read.delim('YBX1_CLIP/YBX1_amir/YBX1_counts.txt', header=F)
YBX1_Counts <- subset(YBX1_Counts, select=c('V1', 'V2', 'V3', 'V6', 'V7'))
colnames(YBX1_Counts) <- c('chr', 'start', 'end', 'strand', 'count')
head(YBX1_Counts)
write.table(YBX1_Counts,'YBX1_CLIP/YBX1_amir/YBX1_counts_hg19.bed',quote = F,row.names = F,col.names = F,sep = '\t')


## liftOver in the terminal >> 352 regions are lost 
# liftOver YBX1_counts_hg19.bed  ~/delaram/liftOver/hg19ToHg38.over.chain YBX1_counts_hg38.bed unmapped_YBX1_counts.txt


## import the hg38 file
YBX1_Counts_hg38 <- read.delim('YBX1_CLIP/YBX1_amir/YBX1_counts_hg38.bed', header=F)
colnames(YBX1_Counts_hg38) <- c('chr', 'start', 'end', 'strand', 'count')
head(YBX1_Counts_hg38)
ggplot(YBX1_Counts_hg38, aes(x=count))+geom_density()+scale_x_log10()+theme_bw()+ggtitle('YBX1 clip count distribution')
summary(YBX1_Counts_hg38$count)


## find the related columns in the design matrix 
designMat <- read.csv('Data/oldDataRefined/DesignMatrices/7_DesignMat_SS_Kmer_DB_ForgiAnnot_Label.csv')

## find the deepbind ID for YBX1 protein 
deepBindDict <- read.csv('Data/oldDataRefined/deepbind/deepBind_dictionary.csv')

deepBindDict$id[deepBindDict$protein_name == 'YBX1'] #  "D00163.001" "D00163.002"
targetInfo <- subset(designMat, select=c('id', 'strand', 'label', 'ev', 'ic', "D00163.001", "D00163.002"))
head(targetInfo)
posInfo <- data.frame(str_split_fixed(targetInfo$id,'_' ,3))
colnames(posInfo) <- c('chr', 'start', 'end')
targetInfo <- cbind(posInfo, targetInfo)

write.table( targetInfo, 'YBX1_CLIP/YBX1_amir/initial_regions.bed',
            quote = F,row.names = F,col.names = F,sep = '\t') #subset(targetInfo,select=c('chr', 'start', 'end', 'strand'))



#### finding intersection and visualization
## How to find intersection with bedtools 
# bedtools intersect -wao -f 0.2 -r -a initial_regions.bed -b  YBX1_counts_hg38.bed  > initialReg_intersect_0.2_YBX1.bed 

.preprocessIntersection <- function(intersect_file, original_bed, clip_bed){
  colnames(intersect_file) <- c(original_bed,paste0(clip_bed,'_clip'),'intersect')
  intersect_file$count_clip[intersect_file$count_clip=='.' ] = 0
  intersect_file$count_clip[intersect_file$count_clip== -1 ] = 0
  intersect_file$count_clip <- as.numeric(intersect_file$count_clip)
  return(intersect_file)
}


.visualizeIntersect <- function(intersect, Title){
  p1=ggplot(subset(intersect,select=c(count_clip,label)), aes(x=label,y=log10(count_clip+1)))+geom_boxplot(aes(fill=label))+theme_bw()+ggtitle(Title)
  p2=ggplot(subset(intersect,select=c(count_clip,label)), aes(color=label,x=log10(count_clip+1)))+geom_density()+theme_bw()+ggtitle(Title)
  p3=ggplot(subset(intersect,select=c(count_clip,label)), aes(x=label,y=log2(count_clip+1)))+geom_violin(aes(fill=label))+theme_bw()+ggtitle(Title)
  
  p4=ggplot(intersect,aes(x=D00163.001,y=log2(count_clip+1),color=label))+geom_point()+theme_bw()+ggtitle(Title)
  p5=ggplot(intersect,aes(x=log2(ev+1),y=count_clip,color=D00163.001))+geom_point()+theme_bw()+scale_color_gradient(low="blue", high="red")+ggtitle(Title)
  
  p6=ggplot(intersect,aes(x=D00163.002,y=log2(count_clip+1),color=label))+geom_point()+theme_bw()+ggtitle(Title)
  p7=ggplot(intersect,aes(x=log2(ev+1),y=count_clip,color=D00163.002))+geom_point()+theme_bw()+scale_color_gradient(low="blue", high="red")+ggtitle(Title)
  
  print(grid.arrange(p1,p2,p3, nrow=1,ncol=3))
  print(grid.arrange(p4,p5, nrow=1,ncol=2))
  print(grid.arrange(p6,p7, nrow=1,ncol=2))
  
}

intersect <- read.delim('YBX1_CLIP/YBX1_amir/initialReg_intersect_0.2_YBX1.bed', header = F)
intersect <- .preprocessIntersection(intersect,colnames(targetInfo) , colnames(YBX1_Counts_hg38))
head(intersect)


pdf('plots/oldDataRefined/YBX1_clip.pdf', width = 12, height=6)
.visualizeIntersect(intersect, 'YBX1 clip')
dev.off()

## comparing 2 versions of YBX1 in deep bind
ggplot(intersect, aes(D00163.001,D00163.002))+geom_point()+theme_bw()
cor(intersect$D00163.001, intersect$D00163.002)





############# called peaks analysis
.cleanTable <- function(tab){
  tab = melt(tab)
  tab$count_clip[tab$count_clip == 1] <- 'Peak'
  tab$count_clip[tab$count_clip == 0] <- 'not-Peak'
  return(tab)
}

## comparing the called peaks -> does this file contain only the called regions?
peaks <-read.delim('YBX1_CLIP/shID-123.uniq.YBX1.CIMS.fdr10.c_hg38.bed', header= F)
colnames(peaks) <- c('chr', 'start', 'end', 'info', 'count', 'strand')

# bedtools intersect -wao -f 0.2 -r -a initial_regions.bed -b  shID-123.uniq.YBX1.CIMS.fdr10.c_hg38.bed  > initialReg_intersect_0.2_YBX1Peaks.bed 
peak_intersect <- read.delim('YBX1_CLIP/YBX1_amir/initialReg_intersect_0.2_YBX1Peaks.bed', header= F)
peak_intersect <- .preprocessIntersection(peak_intersect, colnames(targetInfo) , colnames(peaks))
head(peak_intersect)



intersect_peaks <- peak_intersect
peaks_tab_init <- table(subset(intersect_peaks,select=c(count_clip,label)))


peaks_tab <- .cleanTable(peaks_tab_init)
p1=ggplot(peaks_tab, aes(x=count_clip,y=value,fill=label))+geom_bar(stat = 'identity',color='black')+
  theme_bw()+xlab('')+ylab('count')+scale_fill_manual(values=c('#999999','#E69F00'))+ggtitle('raw data-0.2')

## scaling by peak
pScaled_tab <- peaks_tab_init/rowSums(peaks_tab_init)
pScaled_tab <- .cleanTable(pScaled_tab)
p2=ggplot(pScaled_tab, aes(x=count_clip,y=value,fill=label))+geom_bar(stat = 'identity',color='black')+theme_bw()+
  scale_fill_manual(values=c('#999999','#E69F00'))+xlab('')+ylab('Freq')+ggtitle('scaled-by-peaks')


p3=ggplot(peaks_tab, aes(x=label,y=value,fill=count_clip))+geom_bar(stat = 'identity',color='black')+
  theme_bw()+scale_fill_brewer(palette="Blues")+ggtitle('raw data-0.2')

lScaled_tab <- t(peaks_tab_init)/colSums(peaks_tab_init)
lScaled_tab <- .cleanTable(lScaled_tab)
p4=ggplot(lScaled_tab, aes(x=label,y=value,fill=count_clip))+geom_bar(stat = 'identity',color='black')+
  theme_bw()+scale_fill_brewer(palette="Blues")+xlab('')+ylab('Freq')+ggtitle('scaled-by-labels')

pdf('YBX1_CLIP/called_peaks_analysis.pdf')
grid.arrange(p1, p2, nrow=1,ncol=2)
grid.arrange(p3, p4, nrow=1,ncol=2)
dev.off()




