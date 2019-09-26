# In this file, we will read the refined IC bed file,
#   Try to extract some features from it
#   Then extract the sequences from genome based on start and end ranges,
#   Then label the data, i.e. which sample/sequence is exported via exosome

#   Input: Raw bed files (IC_inrange_18l_500l_20covered_regions.bed) 
#   Output: Features Matrix (Design Matrix)


source('Codes/Functions.R')
Initialize()

library(BSgenome.Hsapiens.UCSC.hg38)

#### Extracting features using annotation files ####
ANNOTATION_DATA_DIR = 'Data/oldDataRefined/Counts_and_Beds'
annotsPath <- list.files(ANNOTATION_DATA_DIR, pattern = '*.bed',full.names = T)
annotFile <- read.delim(annotsPath, header = F)
colnames(annotFile) <- c('seqnames', 'start', 'name', 'name', 'score', 'strand')


## add strand in the id notation and re-generate the designMats 
annotFile$name <- paste0(annotFile$seqnames, '_', annotFile$start, '_', annotFile$end)

head(annotFile)
dim(annotFile)

#### Extract sequencial features ####
hg38 <- BSgenome.Hsapiens.UCSC.hg38

annotFile <- subset(annotFile, seqnames %in% paste0('chr', c(1:22,'X', 'Y')))
dim(annotFile)
smRNAsRange = GRangesForBSGenome(genome = 'hg38',
                                 chrom = annotFile$seqnames,
                                 ranges = IRanges(start = annotFile$start,
                                                  end = annotFile$end),
                                 strand = annotFile$strand)

smRNAsSeq <- getSeq(hg38, smRNAsRange)
smRNAsSeq <- RNAStringSet(complement(smRNAsSeq))

#### Extract GC content of sequences ####
acgtContent = floor(letterFrequency(smRNAsSeq, letters = 'ACGUN', OR = 0, as.prob = F))

designMat = data.frame(id = annotFile$name,
                       chr = annotFile$seqnames,
                       strand = annotFile$strand,
                       seq = as.character(smRNAsSeq),
                       length = width(smRNAsSeq),
                       a = as.integer(acgtContent[,'A']),
                       c = as.integer(acgtContent[,'C']),
                       g = as.integer(acgtContent[,'G']),
                       u = as.integer(acgtContent[,'U']))

head(designMat)
designMat <- designMat[!duplicated(designMat$id),]
write.csv(designMat, file = 'Data/oldDataRefined/DesignMatrices/1_IC_PrimaryDesignMat.csv', quote = F, na = 'NA', row.names = F)


