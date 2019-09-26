
# Required Libraries: parallel, stringr
library(parallel)
library(stringr)

options(stringsAsFactors = F)

SEED = 1
THREADS = detectCores()-2

DESIGN_MATRIX_PREP4 = 'Data/oldDataRefined/DesignMatrices/4_DesignMat_SS_Kmer_Label.csv'
DESIGN_MATRIX_PREP5 = 'Data/oldDataRefined/DesignMatrices/5_DesignMat_SS_Kmer_Deepbind_Label.csv'

seqsLen = nchar(designMatrix$seq)
dotBrsLen = nchar(designMatrix$DB)
identical(seqsLen, dotBrsLen)


source('Codes/Functions.R')
# Initialize()

designMatrix = read.csv(DESIGN_MATRIX_PREP4, stringsAsFactors = F)

## Adding DeepBind scores
# Bash
# cd ~/R/Projects/exosomians; mkdir Apps; cd Apps
# wget http://tools.genes.toronto.edu/deepbind/deepbind-v0.11-linux.tgz
# tar xf deepbind-v0.11-linux.tgz

## Extracting all DeepBind models name related to Homo Sapiens.
# ./deepbind/db/db.tsv

deepbindModels = read.delim('Apps/deepbind/db/db.tsv', comment.char = '#')
# unique(deepbindModels$Species) # "Homo sapiens"
deepbindModels = deepbindModels[deepbindModels$Species=='Homo sapiens', ]
deepbindModelsId = deepbindModels$ID

deepbindParams = list.files('Apps/deepbind/db/params/')
deepbindParams = gsub('.txt', '', deepbindParams)

deepbindAvailableModels = deepbindModelsId[deepbindModelsId%in%deepbindParams]

write.table(deepbindAvailableModels, file = 'Apps/deepbind/models.ids', 
            quote = F, row.names = F, col.names = F)


# Writing sequences
write.table(designMatrix$seq, file = 'Apps/deepbind/sequences.seq', 
            quote = F, row.names = F, col.names = F)

## Not feasible for many sequences
# system('cd Apps/deepbind; ./deepbind models.ids < sequences.seq > deepbind.output')

## Parallelize it by splitting models id into several parts,
# and run deepbind on multicores

# totalNumberOfModels = length(deepbindAvailableModels)
# numberOfModelsForEachThread = floor(totalNumberOfModels/THREADS)

# Bash
# cd ~/R/Projects/exosomians/Apps/deepbind
# mkdir models; cp models.ids models/
# THREADS=`nproc`; THREADS=`expr $THREADS - 2`
# split --additional-suffix=".ids" --number=l/$THREADS models.ids
# rm models.ids
# (CAUTION: Run the code below in a TMUX environment! It might take a while...)
# for model in *.ids; do ../deepbind $model < ../sequences.seq > ${model/.ids/.output} & done
# paste *.output | column -s $'\t' -tn > all.output

## run deepbind script and import the results
# deepbind <- read.table('Data/oldDataRefined/deepbind/deepBind_refined.txt', header = T)

deepbindScores = read.table('Apps/deepbind/models/all.output', header = T)
designMatrix <- cbind(designMatrix, deepbindScores)

write.csv(designMatrix, DESIGN_MATRIX_PREP5, quote = F,row.names = F)

### Mapping deepbind ids to their protein-family names to be more readable 
deepbindDictionary = deepbindModels[deepbindModels$ID %in% deepbindAvailableModels, ]
write.csv(deepbindDictionary, 'Apps/deepbind/deepBindDictionary.csv', quote = F, row.names = F)
