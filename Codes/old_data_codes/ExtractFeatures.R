# In this file, we will read the first annotation file (combined.cluster.scoref10.bed),
#   Try to extract some features from it
#   Then extract the sequences from genome based on start and end ranges,
#   Then read the second annotation file (combined.cluster.scoref10.annotated.bed),
#   Then extract some other features from it
#   Then label the data, i.e. which sample/sequence is exported via exosome

#   Input: Raw annotation files (combined.cluster.scoref10*) 
#   Output: Features Matrix (Design Matrix)


source('Codes/Functions.R')
Initialize()

library(BSgenome.Hsapiens.UCSC.hg38)

#### Extracting features using annotation files ####
ANNOTATION_DATA_DIR = 'Data/Annotation/'
annotsPath = list.files(ANNOTATION_DATA_DIR, full.names = T)
annotFile = read.delim(annotsPath[2], header = F)
head(annotFile)
colnames(annotFile) = c('seqnames', 'start', 'end', 'id', 'unknown', 'strand')

#### Extract sequencial features ####
hg38 = BSgenome.Hsapiens.UCSC.hg38

annotFile = subset(annotFile, seqnames!='chrEBV')
smRNAsRange = GRangesForBSGenome(genome = 'hg38',
                                 chrom = annotFile$seqnames,
                                 ranges = IRanges(start = annotFile$start,
                                                  end = annotFile$end),
                                 strand = annotFile$strand)

smRNAsSeq = getSeq(hg38, smRNAsRange)
names(smRNAsSeq) = annotFile$id

#### Extract GC content of sequences ####
acgtContent = floor(letterFrequency(smRNAsSeq, letters = 'ACGTN', OR = 0, as.prob = T)*100)

identical(annotFile$id, names(smRNAsSeq))
designMat = data.frame(id = annotFile$id,
                       chr = annotFile$seqnames,
                       seq = as.character(smRNAsSeq),
                       length = width(smRNAsSeq),
                       a = as.integer(acgtContent[,'A']),
                       c = as.integer(acgtContent[,'C']),
                       g = as.integer(acgtContent[,'G']),
                       t = as.integer(acgtContent[,'T']),
                       strand = annotFile$strand)

View(head(designMat))


#### Extract other features than sequences ####
annotsPath = list.files(ANNOTATION_DATA_DIR, full.names = T)
annotFileMoreDetailed = read.delim(annotsPath[1], header = F)
annotFileMoreDetailed = annotFileMoreDetailed[, c(4, 10)]
colnames(annotFileMoreDetailed) = c('id', 'annotation')
annotFileMoreDetailed$rnaType = gsub(annotFileMoreDetailed$annotation,
                                     pattern = '-.*', replacement = '')

annotFileMoreDetailed[
  annotFileMoreDetailed$rnaType=='hsa',
  'rnaType'] = 'miRNA'

# annotFileMoreDetailed[
#   grepl(annotFileMoreDetailed$rnaType, pattern = '14q|SNOR|snR|mg'),
#   'rnaType'] = 'snoRNA'

annotFileMoreDetailed[
  !grepl(annotFileMoreDetailed$rnaType, pattern = 'miRNA|tRNA'),
  'rnaType'] = 'snoRNA'

annotFileMoreDetailed = merge(annotFile, annotFileMoreDetailed, all.x = T)

dupIds = annotFileMoreDetailed$id[duplicated(annotFileMoreDetailed$id)]
annotFileMoreDetailed[annotFileMoreDetailed$id %in% dupIds, ]
annotFileMoreDetailed = annotFileMoreDetailed[!duplicated(annotFileMoreDetailed$id), ]
colnames(annotFileMoreDetailed)
annotFileMoreDetailed = annotFileMoreDetailed[, c('id', 'annotation', 'rnaType')]

### Cyto position of the sequences can be added!

#### Merge all annotations ####
designMat = merge(designMat, annotFileMoreDetailed, all.x = T)
write.csv(designMat, file = 'Data/PrimaryDesignMat.csv', quote = F, na = 'NA', row.names = F, col.names = T)


designMat2 = read.csv('Data/PrimaryDesignMat.csv')

