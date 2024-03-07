### This script is used to validate the early screening model 
### based on an independent validation cohort

#################################################################
######## import library && data set #############################
#################################################################
### import library
library(readxl)
library(dplyr)
library(factoextra) 
library(ggplot2)
library(pROC)
library(stats)
library(caret)
library(pheatmap)
library(clusterProfiler)
library(ggstatsplot)
### import the function
source('./used_function.R')
### import data set
lung_vali_77 <- read_excel("/Users/pingyi/Desktop/BDF/data/lung_vali_77.xlsx")
matrix_tss <- read.csv("~/Desktop/XJ_succ/output/data/matrix_tss.csv", row.names=1)
model_clinic2 <- read_excel("/Users/pingyi/Desktop/BDF/data/clinic_vali_H.xlsx")

### make the vali-dataset & save
FT_vali_tss <- cbind(lung_vali_77,matrix_tss[,model_clinic2$sample_name])
dim(FT_vali_tss) # 153
write.csv(FT_vali_tss ,file = "~/Desktop/BDF/data/FT_vali_tss.csv")

#############################################################
######## validate predict model #############################
#############################################################

### import dataset
FT_vali_tss <- read.csv("~/Desktop/BDF/data/FT_vali_tss.csv", row.names=1)
### import the model
load(file = '/Users/pingyi/Desktop/BDF/script/out_data/CAnonCA200f')
#training_range_site,CAnonCAp,result_matrix_mean4, modeltss_CA200,CAtrain,CAtest
FT_vali_tss$site <- paste0('X',rownames(FT_vali_tss))
FT_vali_tss2 <- FT_vali_tss[which(FT_vali_tss$site %in% rownames(training_range_site)),]
FT_vali_tss2 <- data.frame(t(FT_vali_tss2[,c(3:153)]))
### pre_treat ### z-score for sites as training set did
FT_vali_tss2_zscore <- newmatrix_z_norm(new_matrix=FT_vali_tss2,training_range_site)
FT_vali_zscore_range <- range_cal(matrix = FT_vali_tss2_zscore,roworcol = 2)
range(FT_vali_zscore_range$median) # -7.189484 43.905928
### choose the 200 features of model use
rank_sites  = rownames(data.frame(modeltss_CA200$finalModel$importance))
vali_topbins<- FT_vali_tss2_zscore[,rank_sites]
dim(vali_topbins) # 151 200
### model predict based on the external dataset
pred <- predict(modeltss_CA200,newdata = vali_topbins,type = "prob")
### make the prediction table
pred$sample = rownames(pred)
pred$res = 1
pred[which(pred$CA<0.4),'res'] = 0 

pred$true = 1
pred$true[78:151] = 0

pred$equ = 1
pred[which(pred$res != pred$true),'equ'] = 0 

roc_data <- roc(pred$true, pred[,1],type="prob")
print(roc_data)
write.csv(pred,file = '/Users/pingyi/Desktop/BDF/script/out_data/valid_pred3.csv')









