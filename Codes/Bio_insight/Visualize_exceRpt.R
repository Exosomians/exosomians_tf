library(stringr)
library(reshape)
library(reshape2)
library(ggplot2)
library(gridExtra)

load('Data/exceRpt Pipeline - ALL-EV-IC_exceRpt_smallRNAQuants_ReadCounts.RData')
biotypes <- read.table('Data/exceRpt Pipeline - ALL-EV-IC_exceRpt_biotypeCounts.txt')

sumRow <- nrow(exprs.circRNA) + nrow(exprs.gencode) + nrow(exprs.miRNA) + nrow(exprs.piRNA) + nrow(exprs.tRNA)

sampleNames <- str_split_fixed(colnames(exprs.gencode),'_',3)[,2]
sampleNames <- substr(sampleNames, 1,nchar(sampleNames)-2)

colnames(biotypes) <- substr(sampleNames,1,2)
biotypes_melted <- melt(t(biotypes))
colnames(biotypes_melted) <- c('sample', 'rna_type', 'value')

pdf('biotypes_ic_ev.pdf', height=14, width=10)
ggplot(biotypes_melted, aes(x=rna_type,y=value, fill=sample))+
  geom_boxplot(position=position_dodge(1))+coord_flip()+theme_bw()+
  scale_fill_manual(values=c("#56B4E9", "#E69F00"))
dev.off()


### tRNA
colnames(exprs.tRNA) <- substr(sampleNames,1,2)
exprs.tRNA_melted <- melt(t(exprs.tRNA))
colnames(exprs.tRNA_melted) <- c('sample', 'RNA', 'value')
pdf('tRNA_ic_ev.pdf', height=8, width=7)
ggplot(exprs.tRNA_melted, aes(x=RNA,y=value, fill=sample))+
  geom_boxplot(position=position_dodge(1))+coord_flip()+theme_bw()+
  scale_fill_manual(values=c("slateblue3", "pink2"))+ggtitle('tRNA distribution EV vs IC')
dev.off()

### circRNA
colnames(exprs.circRNA) <- substr(sampleNames,1,2)
exprs.circRNA_melted <- melt(t(exprs.circRNA))
colnames(exprs.circRNA_melted) <- c('sample', 'RNA', 'value')
ggplot(exprs.circRNA_melted, aes(x=RNA,y=value, fill=sample))+
  geom_boxplot(position=position_dodge(1))+coord_flip()+theme_bw()+
  scale_fill_manual(values=c("pink2", "#56B4E9"))+ggtitle('circRNA distribution EV vs IC')

### genCode
pdf('geneCode_dstrib_moreEVthanIC.pdf', height = 10, width=7)
colnames(exprs.gencode) <- substr(sampleNames,1,2)
.plotDifRNA(exprs.gencode, 'genCode distrib (EV>IC)')
dev.off()

## piRNA
colnames(exprs.piRNA) <- substr(sampleNames,1,2)
.plotDifRNA(exprs.piRNA, 'piRNA distrib (EV>IC)')

### exprs.miRNA
colnames(exprs.miRNA) <- substr(sampleNames,1,2)
.plotDifRNA(exprs.miRNA, 'miRNA distrib (EV>IC)')




.plotDifRNA <- function(data, title){
  ev = data[,colnames(data)=='EV']
  ic = data[,colnames(data)=='IC']
  diff <- sapply(1:nrow(ev), function(i) quantile(ev[i,],0.75)-quantile(ic[i,],0.75))
  data_filter <- data[diff>5,]
  
  data_melted <- melt(t(data_filter))
  
  colnames(data_melted) <- c('sample', 'RNA', 'value')
  ggplot(data_melted, aes(x=RNA,y=value, fill=sample))+
    geom_boxplot(position=position_dodge(1))+coord_flip()+theme_bw()+
    scale_fill_manual(values=c("pink3", "#999999")) + ggtitle(title)
}

