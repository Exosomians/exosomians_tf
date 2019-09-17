#  In this file, we will try to find the ngrams(K-mers) 
#     of RNA sequences to use them as features. 
#     We'll make all the possible 4-grams and 2-grams
#     then, we'll check the count of the ngrams in each RNA sequence



source('Codes/Functions.R')
Initialize()
designMat = read.csv('Data/oldDataRefined/DesignMatrices/3_DesignMat_SS_Label.csv', stringsAsFactors = F)
dim(designMat)

K_MER_SIZE = 4
All_Possible_Kmers <- Make_All_Possible_Kmers(K_MER_SIZE, 4, c('A','U','C','G'))


RNA_NucleotideKmers <- sapply(as.character(designMat$seq), makeKmerForNucleotide , K_MER_SIZE)
summary(sapply(RNA_NucleotideKmers, length))


Nucleotide_KmerMatrix <- MakeFeatureSpecificMatrix(All_Possible_Kmers, RNA_NucleotideKmers, designMat$id)
DrawFeatureDistribution(Nucleotide_KmerMatrix)

summary(rowSums(Nucleotide_KmerMatrix))
write.csv(Nucleotide_KmerMatrix,'Data/oldDataRefined/Ngrams/ngramMatrix.csv',quote = F)

designMat <- cbind(designMat, Nucleotide_KmerMatrix)
# write.csv(designMat,'Data/oldDataRefined/DesignMatrices/4_DesignMat_SS_Kmer_Label.csv',quote = F,row.names = F)
designMat <- read.csv('Data/oldDataRefined/DesignMatrices/4_DesignMat_SS_Kmer_Label.csv', stringsAsFactors = F)



######### adding deepbind scores
## writing sequences
write.table(designMat$seq, file = '~/Sequences.seq',quote = F, row.names = F, col.names = F)
## run deepbind script and import the results
deepbind <- read.table('Data/oldDataRefined/deepbind/deepBind_refined.txt', header = T)

designMat <- cbind(designMat, deepbind)
write.csv(designMat,'Data/oldDataRefined/DesignMatrices/5_DesignMat_SS_Kmer_DB_Label.csv',quote = F,row.names = F)

### mapping deepbind ids to their protein-family names to be more readable 
DeepBind_ids <- read.delim('Data/oldDataRefined/deepbind/models.ids', header=F)
DeepBind_ids <- subset(data.frame(str_split_fixed(DeepBind_ids$V1,' ' ,3)), select=-X2)
names(DeepBind_ids) <- c('id', 'protein_name')
write.csv(DeepBind_ids, 'Data/oldDataRefined/deepbind/deepBind_dictionary.csv',quote = F,row.names = F)

