library(TCGAbiolinks)
library(SummarizedExperiment)
library(dplyr)
library(DT)
library(tibble)
rm(list=ls())

setwd('D:\\test\\Five_key_lncRNAs\\chapter1\\download\\expression')

exp <- GDCquery(project = "TCGA-PAAD", 
                data.category = "Transcriptome Profiling", 
                data.type = "Gene Expression Quantification", 
                workflow.type = "HTSeq - FPKM")

GDCdownload(exp)

expdat <- GDCprepare(query =exp)

count_matrix= as.data.frame(assay(expdat))




require("rtracklayer")
require("SummarizedExperiment")
setwd("D:\\test\\Five_key_lncRNAs\\chapter1")



gtf1 <- rtracklayer::import('gencode.v30.long_noncoding_RNAs.gtf')

head(gtf1$gene_id ) 

head(gtf1$gene_type)




library(tibble)

a<-cbind(gtf1$gene_id,gtf1$gene_name,gtf1$gene_type)

colnames(a)<-c("gene_id","gene_name","gene_type")

gtf_df <- as.data.frame(a)

gtf_df <- gtf_df %>% 
  tidyr::separate(gene_id,into = c("gene_id","drop"),sep="\\.") %>% 
  dplyr::select(-drop) 


gtf_df <- gtf_df %>% 
  tidyr::unite(gene_id,gene_name,gene_id,gene_type,sep = " | ")%>% 
   dplyr::distinct() %>% 
   tidyr::separate(gene_id, c("gene_name","gene_id","gene_biotype"), sep = " \\| ")






fpkmToTpm <- function(fpkm)
{
  exp(log(fpkm) - log(sum(fpkm)) + log(1e6))
}


expr <- as.data.frame (apply(count_matrix  , 2, fpkmToTpm))

library(tibble)

expr <-  expr %>% rownames_to_column("gene_id")


# =======================================================


ncRNA <- c("sense_overlapping","lincRNA","3prime_overlapping_ncRNA",
           "processed_transcript","sense_intronic","bidirectional_promoter_lncRNA",
           "non_coding")

ncRNA %in% unique(gtf_df$gene_biotype)



LncRNA_exprSet <- gtf_df %>% 
  dplyr::inner_join(expr,by ="gene_id") %>%
  distinct(gene_name, .keep_all = TRUE)


# 
# write.csv(LncRNA_exprSet,file='LncRNA_exprSet.csv')
# 

save(LncRNA_exprSet,file = "LncRNA_exprSet.Rda")


# =======================================================



load("LncRNA_exprSet.Rda")
library(tidyverse)
LncRNA_exprSet <- LncRNA_exprSet %>% 
  remove_rownames %>%
  column_to_rownames(var="gene_name")


LncRNA_exprSet$gene_id <- NULL
LncRNA_exprSet$gene_biotype <- NULL


#=======================================================


metadata <- data.frame(colnames(LncRNA_exprSet))


for (i in 1:length(metadata[,1])) {
  num <- as.numeric(as.character(substring(metadata[i,1],14,15)))
  if (num == 1 ) {metadata[i,2] <- "T"}
  if (num != 1) {metadata[i,2] <- "N"}
}
names(metadata) <- c("id","group")
metadata$group <- as.factor(metadata$group)



#=======================================================


group = factor(metadata$group,levels=c('N','T'))

design <- model.matrix(~group)

library(edgeR)

y <- DGEList(counts=LncRNA_exprSet,group=group)
y <- calcNormFactors(y)
y <- estimateCommonDisp(y)
y <- estimateTagwiseDisp(y)
et <- exactTest(y,pair = c('N','T'))
topTags(et)
ordered_tags <- topTags(et, n=100000)

allDiff=ordered_tags$table
allDiff=allDiff[is.na(allDiff$FDR)==FALSE,]
diff=allDiff
write.csv(diff,file="edgerOut.csv")
newData=y$pseudo.counts


foldChange=1
padj=0.05
getwd()

diffSig <- diff[(diff$PValue  < padj & (diff$logFC > foldChange | diff$logFC < (-foldChange))),]
write.csv(diffSig,file="diffSig.csv")
diffUp <- diff[(diff$FDR < padj & (diff$logFC > foldChange)),]
write.csv(diffUp,file="up.csv")
diffDown <- diff[(diff$FDR < padj & (diff$logFC < (-foldChange))),]
write.csv(diffDown,file="down.csv")

normalizeExp <- rbind(id=colnames(newData),newData)
diffExp <- rbind(id=colnames(newData),newData[rownames(diffSig),])
diffExp <- rbind(id=colnames(newData),newData[rownames(diffSig),])
write.table(diffExp,file="diffmRNAExp.txt",sep="\t",quote=F,col.names=F)



