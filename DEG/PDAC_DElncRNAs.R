library(GDCRNATools)
library(tidyverse)
library(caret)
library(DMwR)

options(stringsAsFactors = F)

####### download data
setwd('D:\\PAAD\\111data')
project <- 'TCGA-PAAD'
rnadir <-paste(project, 'RNAseq', sep='/')#dir
####### Download RNAseq data #######
gdcRNADownload(project.id ='TCGA-PAAD',data.type='RNAseq',write.manifest =FALSE,directory= rnadir)


#######clinical data
metaMatrix.RNA <-gdcParseMetadata(project.id ='TCGA-PAAD',data.type  ='RNAseq',write.meta =FALSE)
####### Filter duplicated samples in RNAseq metadata #######
metaMatrix.RNA <-gdcFilterDuplicate(metaMatrix.RNA)
####### Filter non-Primary Tumor and non-Solid Tissue Normal samples in RNAseq metadata #######
metaMatrix.RNA <-gdcFilterSampleType(metaMatrix.RNA)


####### Merge RNAseq data #######
rnaCounts <-gdcRNAMerge(metadata  = metaMatrix.RNA, path = rnadir, data.type ='RNAseq')
####### Normalization of RNAseq data #######
rnaExpr <-gdcVoomNormalization(counts = rnaCounts, filter =FALSE)


####### clinical data
metaMatrix.RNA$sample_type <- factor(metaMatrix.RNA$sample_type)
count <- as.data.frame(t(rnaCounts))
metadata <- metaMatrix.RNA[,-c(1,2,3,5,6,8,9,10,11,12,13,14,15)]


####### merge
count <- count %>% 
  rownames_to_column(var = 'sample') %>% 
  inner_join(metadata,by = 'sample') %>% 
  column_to_rownames(var = 'sample')


####### SMOTE
set.seed(20)
smote_train <- SMOTE(sample_type ~ ., data  = count)
smote_train <- na.omit(smote_train)
table(smote_train$sample_type)

metadata_smote <- smote_train[,ncol(smote_train)]
input_smote <- as.data.frame(t(smote_train[,-ncol(smote_train)]))

####### upsample
set.seed(20)
up_train <- upSample(x = count[, -ncol(count)],
                     y = count$sample_type,yname = T)   
up_train <- na.omit(up_train)#
table(up_train[,ncol(up_train)]) 

metadata_up <- up_train[,ncol(up_train)]
input_up <- as.data.frame(t(up_train[,-ncol(up_train)]))


####### differential analysis
DEGAll_smote <-gdcDEAnalysis(counts= input_smote, 
                       group=metadata_smote, 
                       comparison ='PrimaryTumor-SolidTissueNormal',
                       method='edgeR',filter = F)

DEGAll_up <-gdcDEAnalysis(counts= input_up, 
                       group=metadata_up, 
                       comparison ='PrimaryTumor-SolidTissueNormal',
                       method='edgeR',filter = F)



#######  DEGs
diffSig_smote <- DEGAll_smote[(DEGAll_smote$group == 'long_non_coding' & (DEGAll_smote$PValue  < 0.05 & (DEGAll_smote$logFC > 1 | DEGAll_smote$logFC < (-1)))),]
diffSig_up <- DEGAll_up[(DEGAll_up$group == 'long_non_coding' & (DEGAll_up$PValue  < 0.05 & (DEGAll_up$logFC > 1 | DEGAll_up$logFC < (-1)))),]


rnaExpr <- as.data.frame(rnaExpr)
rnaExpr <- rnaExpr %>% 
  rownames_to_column(var = 'id')
diffSig_smote <- diffSig_smote %>% 
  rownames_to_column(var = 'id')
diffSig_up <- diffSig_up %>% 
  rownames_to_column(var = 'id')

diffmRNAExp_smote <- diffSig_smote %>% 
  inner_join(rnaExpr,by = 'id')
diffmRNAExp_up <- diffSig_up %>% 
  inner_join(rnaExpr,by = 'id')


####### 
diffmRNAExp <- diffmRNAExp_smote[diffmRNAExp_smote$id %in% diffmRNAExp_up$id,]
write.csv(diffmRNAExp,file = 'difflncRNAExp.csv',quote = F)#
