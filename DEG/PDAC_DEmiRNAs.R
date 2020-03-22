library(GDCRNATools)
library(tidyverse)
library(caret)
library(DMwR)

options(stringsAsFactors = F)

####### download data
setwd('D:\\PAAD\\111data')
project <- 'TCGA-PAAD'
mirdir <-paste(project, 'miRNAs', sep='/')#dir
####### Download RNAseq data #######
gdcRNADownload(project.id ='TCGA-PAAD',data.type='miRNAs',write.manifest =FALSE,directory= mirdir)


####### Parse miRNAs metadata #######
metaMatrix.MIR <-gdcParseMetadata(project.id ='TCGA-PAAD',data.type  ='miRNAs',write.meta =FALSE)
####### Filter duplicated samples in miRNAs metadata #######
metaMatrix.MIR <-gdcFilterDuplicate(metaMatrix.MIR)
####### Filter non-Primary Tumor and non-Solid Tissue Normal samples in miRNAs metadata #######
metaMatrix.MIR <-gdcFilterSampleType(metaMatrix.MIR)


####### Merge miRNAs data #######
mirCounts <-gdcRNAMerge(metadata  = metaMatrix.MIR,path= mirdir,data.type ='miRNAs')
####### Normalization of miRNAs data #######
mirExpr <-gdcVoomNormalization(counts = mirCounts, filter =FALSE)

####### clinical data
metaMatrix.MIR$sample_type <- factor(metaMatrix.MIR$sample_type)
count <- as.data.frame(t(mirCounts))
metadata <- metaMatrix.MIR[,-c(1,2,3,5,6,8,9,10,11,12,13,14,15)]#保留分组信息


####### merge
count <- count %>% 
  rownames_to_column(var = 'sample') %>% 
  inner_join(metadata,by = 'sample') %>% 
  column_to_rownames(var = 'sample')


####### SMOTE
set.seed(14)
smote_train <- SMOTE(sample_type ~ ., data  = count)
smote_train <- na.omit(smote_train)#删除NA行
table(smote_train$sample_type)

metadata_smote <- smote_train[,ncol(smote_train)]
input_smote <- as.data.frame(t(smote_train[,-ncol(smote_train)]))
DEGAll_smote <-gdcDEAnalysis(counts= input_smote, #Counts数据
                             group=metadata_smote, #样本分组信息
                             comparison ='PrimaryTumor-SolidTissueNormal',#比较组信息 
                             method='edgeR',filter = F)
DEGAll_smote <- DEGAll_smote %>% 
  rownames_to_column(var = 'id')

####### upsample
set.seed(14)
up_train <- upSample(x = count[, -ncol(count)],
                     y = count$sample_type,yname = T)   
up_train <- na.omit(up_train)#删除NA行
table(up_train[,ncol(up_train)]) 

metadata_up <- up_train[,ncol(up_train)]
input_up <- as.data.frame(t(up_train[,-ncol(up_train)]))


####### differential analysis
DEGAll_smote <-gdcDEAnalysis(counts= input_smote, #Counts数据
                       group=metadata_smote, #样本分组信息
                       comparison ='PrimaryTumor-SolidTissueNormal',#比较组信息 
                       method='edgeR',filter = F)
DEGAll_smote <- DEGAll_smote %>% 
  rownames_to_column(var = 'id')

DEGAll_up <-gdcDEAnalysis(counts= input_up, #Counts数据
                       group=metadata_up, #样本分组信息
                       comparison ='PrimaryTumor-SolidTissueNormal',#比较组信息 
                       method='edgeR',filter = F)


####### DEG
diffSig_smote <- DEGAll_smote[(DEGAll_smote$PValue  < 0.05 & (DEGAll_smote$logFC > 1 | DEGAll_smote$logFC < (-1))),]
diffSig_up <- DEGAll_up[(DEGAll_up$PValue  < 0.05 & (DEGAll_up$logFC > 1 | DEGAll_up$logFC < (-1))),]


rnaExpr <- as.data.frame(mirExpr)
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


####### interaction
diffmRNAExp <- diffmRNAExp_smote[rownames(diffmRNAExp_smote) %in% rownames(diffmRNAExp_up),]
write.csv(diffmRNAExp,file = 'diffmiRNAExp.csv',quote = F)#
