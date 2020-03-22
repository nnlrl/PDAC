# =======================================================
#########################################################
#########################################################
#作者：工程师2号 ########################################
#简书笔记博客（柳叶刀与小鼠标）##########################
# https://www.jianshu.com/u/619b87e54936 ###############
#########################################################
#########################################################
# =======================================================




# =======================================================
####step1  下载临床数据包################################
# =======================================================


library(TCGAbiolinks)
library(SummarizedExperiment)
library(dplyr)
library(DT)
library(tibble)
rm(list=ls())


setwd('D:\\PAAD\\111data')
rm(list=ls())
query <- GDCquery(project = "TCGA-PAAD", 
                  data.category = "Clinical", 
                  file.type = "xml")
GDCdownload(query)
clinical <- GDCprepare_clinic(query, clinical.info = "patient")



clinical_trait <- clinical  %>%
  dplyr::select(bcr_patient_barcode,gender,vital_status,                            
                days_to_death,days_to_last_followup,race_list,
                person_neoplasm_cancer_status,age_at_initial_pathologic_diagnosis,
             neoplasm_histologic_grade,stage_event_pathologic_stage,             
                stage_event_tnm_categories  ) %>%
  distinct( bcr_patient_barcode, .keep_all = TRUE)  


#整理死亡患者的临床信息
dead_patient <- clinical_trait  %>%
  dplyr::filter(vital_status == 'Dead') %>%
  dplyr::select(-days_to_last_followup) %>%
  reshape::rename(c(bcr_patient_barcode = 'Barcode',
                    gender = 'Gender',
                    vital_status = 'OS',
                    days_to_death='OS.Time',
                    race_list = 'Race',
                    person_neoplasm_cancer_status='cancer_status',
                    age_at_initial_pathologic_diagnosis = 'Age',
                    neoplasm_histologic_grade = 'Grade',
                    stage_event_pathologic_stage = 'Stage',
                    stage_event_tnm_categories = 'TNM' )) %>%
  mutate(OS=ifelse(OS=='Dead',1,0))%>%
  mutate(OS.Time=OS.Time/365)



#整理生存患者的临床信息
alive_patient <- clinical_trait %>%
  dplyr::filter(vital_status == 'Alive') %>%
  dplyr::select(-days_to_death) %>%
  reshape::rename(c(bcr_patient_barcode = 'Barcode',
                    gender = 'Gender',
                    vital_status = 'OS',
                    days_to_last_followup='OS.Time',
                    race_list = 'Race',
                    person_neoplasm_cancer_status='cancer_status',
                    age_at_initial_pathologic_diagnosis = 'Age',
                    neoplasm_histologic_grade = 'Grade',
                    stage_event_pathologic_stage = 'Stage',
                    stage_event_tnm_categories = 'TNM' )) %>%
  mutate(OS=ifelse(OS=='Dead',1,0))%>%
  mutate(OS.Time=OS.Time/365)

#合并两类患者，得到肾透明细胞癌的临床信息
survival_data <- rbind(dead_patient,alive_patient)

write.csv(survival_data , file = 'survival.csv')




# =======================================================
####step2 单因素#########################################
# =======================================================


library(survival)
library(future.apply)

plan(multiprocess)
setwd('D:\\PAAD\\cox')

rm(list=ls())

#diffmRNAExp.txt文件是上一节生成的差异基因文件

univariate_data <- read.table('difflncRNAExp.csv',header = T,row.names = 1,
                              sep = ',')

#univariate_data  <- log2(univariate_data +0.001)
univariate_data <- univariate_data %>%
  remove_rownames() %>% 
  column_to_rownames(var = 'symbol')
univariate_data <- univariate_data[,-c(1:7)]

metadata <- data.frame(colnames(univariate_data))





for (i in 1:length(metadata[,1])) {
  num <- as.numeric(as.character(substring(metadata[i,1],14,15)))
  if (num == 1 ) {metadata[i,2] <- "T"}
  if (num != 1) {metadata[i,2] <- "N"}
}
names(metadata) <- c("id","group")
metadata$group <- as.factor(metadata$group)


metadata <- subset(metadata,metadata$group == "T")
metadata

exprSet <- univariate_data[,which(colnames(univariate_data) %in% metadata$id)]

colnames(exprSet)  <- substr(x=colnames(exprSet),start = 1,stop = 12)
colnames(exprSet)  <- chartr(old='.',new = '-',x=colnames(exprSet) )

survival <- read.csv('survival.csv',header = T,row.names = 1)




select_colname <- function(expr,name){
  expr <- expr[,which(colnames(expr) %in% name[,1])]
  expr  <- expr[,order(names(expr))]   
  name <- name[which(name[,1] %in% colnames(expr)),]
  name <- name[order(name[,1]),]
  expr_name <- list( expr,name)
  return( expr_name)
}

dat <- select_colname(exprSet,survival )


data <- as.data.frame(t(dat[[1]]))


data$Barcode <- row.names(data)
survival <- survival %>%
  dplyr::select(Barcode,OS,OS.Time)


survival_dat  <- merge(survival,data,by='Barcode')

#这一步是为了改变基因名，如果不改基因名的话会报错
#原始的基因名，包含‘-’，‘.‘，’：‘这样的特殊符号，
#这样的特殊符号会导致下面的运算出错
#所以我们统一将这些符号改成’_‘

colnames(survival_dat) <- gsub("\\-", "", colnames(survival_dat))

colnames(survival_dat) <- gsub("\\.", "_", colnames(survival_dat))

colnames(survival_dat) <- gsub("\\:", "_", colnames(survival_dat))

covariates <- as.character(colnames(survival_dat))[-(1:3)]

univ_formulas <- sapply(covariates,
                        function(x){
                          as.formula(paste('Surv(OS_Time, OS)~', x))})

univ_models <- future_lapply( univ_formulas,
                              function(x){coxph(x, data = survival_dat,x = T,y = T)})


univ_results <-  lapply(univ_models,
                        function(x){ 
                          x <- summary(x)
                          p.value <- signif(x$wald["pvalue"], digits = 2)
                          beta <- signif(x$coef[1], digits = 2)
                          HR <- signif(x$coef[2], digits = 2)
                          HR.confint.lower <- signif(x$conf.int[,"lower .95"], 2)
                          HR.confint.upper <- signif(x$conf.int[,"upper .95"], 2)
                          HR <- paste0(HR, " (", 
                                       HR.confint.lower, "-", HR.confint.upper, ")")
                          res <- c(beta, HR, p.value)
                          names(res) <- c("coef", "HR (95% CI for HR)", "p.value")
                          return(res)
                        })


res <- as.data.frame(t(do.call(cbind, univ_results)))

res[1:6,]

res$p.value <- as.numeric(as.character(res$p.value))

res$coef <- as.numeric(as.character(res$coef ))

res[1:6,]

table(res$p.value < 0.001)

res  <- res[res$p.value <= 0.001, ]

res  <- res[order(res$p.value), ]

single_gene <- rownames(res)
#single_gene即为得到的单因素cox分析中P小于0.05的基因


# =======================================================
####step3 多因素cox分析##################################
# =======================================================



survival_data <- survival_dat %>%
  select(OS,OS_Time, !!single_gene,)

cox <- coxph(Surv(OS_Time, OS) ~ ., data = survival_data,x = T,y = T)

#我这里使用的是‘both’，原作者使用的是后向

cox=step(cox,direction = "both")


tmp <- summary(cox)

coxgene <- rownames(tmp$coefficients)

# =======================================================
####使用LOO-BootStrap CV#################################
# =======================================================
library(survival)
library(riskRegression)

x1 = Score(list('TCGA-PDAC'  = cox),formula=Surv(OS_Time, OS) ~ 1,
           data=survival_data,times=5,split.method = 'loob',
           plots=c("roc","calibration"),metrics="auc",seed = 2)

plotCalibration(x1)

cindex <- cindex(list('TCGA-PDAC'  = cox),formula=Surv(OS_Time, OS) ~ .,
                 data=survival_data,eval.times=seq(1,5,0.1),split.method = 'bootcv',
                 plots=c("roc","calibration"),metrics="auc",seed = 2)

plot(cindex,xlim = c(0,5))

#   times      model  AUC se.conservative lower upper
#4:     1  TCGA-PDAC 78.6            0.08  69.3  87.9
#5:     3  TCGA-PDAC 79.9            0.20  66.6  93.2
#6:     5  TCGA-PDAC 83.8            0.31  65.6 100.0

# Efron's leave-one-out-bootstrap based on 100 bootstrap samples (drawn with replacement) each of size 177.
# The 'confidence intervals' and 'p-values' are obtained with the delta method after bootstrap.

# =======================================================
####step4 lasso分析######################################
# =======================================================




library(glmnet)


x <- data.matrix(survival_dat[,4:dim(survival_dat)[2]])

time=as.double(survival_dat$OS_Time)
status=as.double(survival_dat$OS)


y <- data.matrix(Surv(time,status))
cv.fit <- cv.glmnet(x, y, family="cox",maxit = 1000,nfolds = 5)
plot(cv.fit)
cv.fit$lambda.min
cv.fit$lambda.1se


coef.min = coef(cv.fit, s = "lambda.min")


# requires tibble.
tidy_coef <- function(x){
  coef(x,s = "lambda.min") %>%
    as.matrix  %>%   # Coerce from sparse matrix to regular matrix.
    data.frame %>%  # Then dataframes.
   rownames_to_column %>%  # Add rownames as explicit variables.
    setNames(c("term","estimate")) %>%
    filter(estimate != 0)
}


coef_data <- tidy_coef(cv.fit)
coef_data 


multi_gene <- coef_data$term


######################################################################################

#其实不管是基于AIC信息的前向或者后向cox，或者基于lasso分析，两个分析的目的
#均是为了筛选合适作为预测模型的基因，两者并无优劣。

######################################################################################




# =======================================================
###绘图##################################################
# =======================================================




survival_data <- survival_dat %>%
  remove_rownames() %>%
  column_to_rownames(var = 'Barcode')%>%
  select(OS,OS_Time, !!coxgene)



prognostic_index = predict(cox,type="risk",newdata=survival_data)


index =as.vector(ifelse(prognostic_index>median(prognostic_index),"high","low"))

index_data <- cbind(survival_data,prognostic_index,
                                      index)


index_data



# =======================================================
###FIGURE3A##################################################
# =======================================================



library(ggplot2)
a <- index_data[order(index_data$prognostic_index),]
a$id <- 1:length(rownames(a))

ap <- ggplot(a,aes(x=id,y=prognostic_index))+
  geom_point(aes(colour = factor(index))) + ylim(0,3)+

 theme(legend.title=element_blank(),legend.justification=c(1,.045),
           legend.position=c(1,.045),legend.key.size=unit(1.5,'cm')
)+geom_point(aes(y=prognostic_index, x=id,
                 colour = factor(index)),
             size=3)+theme(panel.grid =element_blank())+
  theme(axis.text.x = element_blank()) +
  theme(axis.ticks.x = element_blank())+
  theme(panel.border = element_blank())+
  theme(axis.line = element_line(size=1))


pdf('FIGURE3A.pdf',height = 4,width = 4)
ap
dev.off()



# =======================================================
###FIGURE3B##################################################
# =======================================================

a <- index_data[order(index_data$prognostic_index),]
a$id <- 1:length(rownames(a))
a$OS  <- ifelse(a$OS == 0, 'Alive','Dead')
bp <-qplot(id,OS_Time,data=a,col = as.factor(OS))
theme_set(theme_bw())

bp <- bp + theme(legend.title=element_blank(),legend.justification=c(1,1),
           legend.position=c(1,1),legend.key.size=unit(1.5,'cm')
)+geom_point(aes(y=OS_Time, x=id),size=3)+theme(panel.grid =element_blank())+
  theme(axis.text.x = element_blank()) +theme(axis.ticks.x = element_blank())+
  theme(panel.border = element_blank())+ylab('SurvivalTime(Year)')+xlab('')+
  theme(axis.line = element_line(size=1, colour = "black"))



pdf('FIGURE3B.pdf',height = 4,width = 4)
bp
dev.off()






# =======================================================
###FIGURE3C##################################################
# =======================================================

heatmap <- a  %>%
  dplyr::select(!!multi_gene,index)



str(heatmap)


annotation_col  <- data.frame(row.names(heatmap),heatmap$index)
row.names(annotation_col) <- annotation_col$row.names.heatmap.
annotation_col$row.names.heatmap. <- NULL
colnames(annotation_col) <- "Subtype"
heatmap$index <- NULL
heatmap <- as.data.frame(t(heatmap))



annotation_col$Subtype <- as.factor(annotation_col$Subtype)
Subtype = c("#BC3C28", "#0072B5")
names(Subtype) = c("high", "low")
ann_colors = list(Subtype = Subtype )
# bk = unique(c(seq(-4,4, length=100)))


pdf('FIGURE3C.pdf',height = 4,width = 4)
library(pheatmap)
par(oma=c(2,2,1,2),mar=c(5,4,4,2))
p1 = pheatmap(heatmap, 
              scale = 'row',
              cluster_cols = FALSE,
              # breaks = bk,
              cluster_row = T,
              annotation_colors = ann_colors,
              show_colnames     = FALSE,
              show_rownames     = FALSE,
              color = colorRampPalette(c("#DC143C",
                                         "#F8F8FF",
                                         "#000080"))(100), 
              legend = T,
              annotation_col = annotation_col,
              annotation_legend = T,
              treeheight_row=0,
              treeheight_col=0,
              fontsize = 16,
              fontsize_row=6, 
              fontsize_col = 6)
dev.off()





# =======================================================
###FIGURE5###############################################
# =======================================================



library(survivalROC)
pdf(file="FIGURE5.pdf")
par(oma=c(0.5,1,0,1),font.lab=1.5,font.axis=1.5)
roc=survivalROC.C(Stime=index_data$OS_Time, 
                  status=index_data$OS, 
                  marker = index_data$prognostic_index, 
                  predict.time =5)



plot(roc$FP, roc$TP, type="l", xlim=c(0,1), ylim=c(0,1),col='red', 
     xlab="False positive rate", ylab="True positive rate",
     main=paste("ROC curve (", "AUC = ",round(roc$AUC,3),")"),
     lwd = 2, cex.main=1.3, cex.lab=1.2, cex.axis=1.2, font=1.2)
abline(0,1)


dev.off()


# =======================================================
###FIGURE6###############################################
# =======================================================

library(survminer)



fit  <- survfit(Surv(OS_Time, OS) ~ index, data = index_data,type = 'fh2')


pdf(file="FIGURE6.pdf",height = 4, width = 4)
par(oma=c(2,2,1,2),mar=c(5,4,4,2))

ggsurv <- ggsurvplot(fit, data = index_data,
                     pval = T,
                     xlim = c(0,5),  
                     break.time.by = 1,  
                     xlab = "Time in years",
                     palette = c( "#377EB8",
                                  "#42B540",
                                  "#E41A1C"),
                     legend.labs =c("Low Expression",
                                    "High Expression"))



ggsurv <- ggpar( ggsurv,
                 font.y  = c(16, "bold"), 
                 font.x  = c(16, "bold"),
                 legend = "top",
                 font.legend = c(6, "bold"))

ggsurv
dev.off()


# =======================================================
#########################################################
#########################################################
#作者：工程师2号 ########################################
#简书笔记博客（柳叶刀与小鼠标）##########################
# https://www.jianshu.com/u/619b87e54936 ###############
#########################################################
#########################################################
# =======================================================



