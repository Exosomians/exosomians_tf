# In this file, we will read the refined IC bed file,
#   Try to extract some features from it
#   Then extract the sequences from genome based on start and end ranges,
#   Then label the data, i.e. which sample/sequence is exported via exosome

#   Input: Raw bed files (IC_inrange_20l_500l_20covered_regions.bed) 
#   Output: Features Matrix (Design Matrix)


source('Codes/Functions.R')
Initialize()

library(BSgenome.Hsapiens.UCSC.hg38)

#### Extracting features using annotation files ####
ANNOTATION_DATA_DIR = 'Data/oldDataRefined'
annotsPath = list.files(ANNOTATION_DATA_DIR, pattern = '*.bed',full.names = T)
names(annotsPath) <- substr(list.files(ANNOTATION_DATA_DIR, pattern = '*.bed'),1,2)
annotFile = read.delim(annotsPath[['IC']], header = F)
colnames(annotFile) = c('seqnames', 'start', 'end')
head(annotFile)
dim(annotFile)

#### Extract sequencial features ####
hg38 = BSgenome.Hsapiens.UCSC.hg38

annotFile = subset(annotFile, seqnames %in% paste0('chr', c(1:22,'X', 'Y')))
dim(annotFile)
smRNAsRange = GRangesForBSGenome(genome = 'hg38',
                                 chrom = annotFile$seqnames,
                                 ranges = IRanges(start = annotFile$start,
                                                  end = annotFile$end),
                                 strand = annotFile$strand)

smRNAsSeq = getSeq(hg38, smRNAsRange)


#### Extract GC content of sequences ####
acgtContent = floor(letterFrequency(smRNAsSeq, letters = 'ACGTN', OR = 0, as.prob = T)*100)

designMat = data.frame(id = paste0(annotFile$seqnames, '_', annotFile$start, '_', annotFile$end),
                       chr = annotFile$seqnames,
                       seq = as.character(smRNAsSeq),
                       length = width(smRNAsSeq),
                       a = as.integer(acgtContent[,'A']),
                       c = as.integer(acgtContent[,'C']),
                       g = as.integer(acgtContent[,'G']),
                       t = as.integer(acgtContent[,'T']))

head(designMat)
write.csv(designMat, file = 'Data/oldDataRefined/1_PrimaryDesignMat.csv', quote = F, na = 'NA', row.names = F)


