### try to normalize with CA non-CA model mean/sd
### import the library 
library(readxl)
library(dplyr)
library(factoextra) 
library(ggplot2)
library(pROC)
library(stats)
library(caret)
library(pheatmap)
library(clusterProfiler)
library(org.Hs.eg.db)
library(msigdbr)
library(ggrepel)
### import the function
source('~/Desktop/BDF/script/BD_lung/used_function.R')
source('/Users/pingyi/Desktop/BDF/script/BD_lung/function_path.R')

####################### try #######################
load(file = '~/Desktop/BDF/script/out_data/CAnonCA200f') 
training_range_siteCA = training_range_site
### try training_range_site to normalize all data (train && validate)
####################### try end####################
#load(file = '~/Desktop/BDF/script/out_data/lung_tss_pred')
lung_clinic_withvar <- read.csv("~/Desktop/BDF/script/out_data/lung_clinic_withvar.csv", row.names=1)
bridge_matrix <- read.csv("~/Desktop/BDF/script/out_data/bridge_matrix.csv")
KEGG_path <-  read.gmt("~/Desktop/BDF/data/public/path/NV/KEGG/c2.cp.kegg_medicus.v2023.2.Hs.symbols.gmt")
length(unique(KEGG_path$term)) # 619
lung_clinic <- lung_clinic_withvar
#dim(lung_tss_pred) #  196 90808
dim(lung_clinic) # 65 17
lung_clinic$group = lung_clinic$EGFR
CA_norm <- rbind(CAtrain,CAtest)
rownames(CA_norm)
lung_tss_pred <- CA_norm[lung_clinic$sample,]
dim(lung_tss_pred) # 65 90788

### 70% for training 30% for testing
clinic = lung_clinic
matrix = lung_tss_pred
set.seed(123)
inTrain=createDataPartition(y=clinic$group,p=0.70,list=FALSE)
#set.seed(1349)
#inTrain=createDataPartition(y=clinic$group,p=0.60,list=FALSE)

training=matrix[inTrain,]
testing=matrix[-inTrain,]
### save the group of training && testing 
trainLabels <- clinic[inTrain,'group']
testLabels <- clinic[-inTrain,'group']
trainsample <- clinic[inTrain,]
testsample <- clinic[-inTrain,]

table(trainLabels)
table(testLabels)

test_clinic_matrix<- clinic[-inTrain,]
dim(training) #    47 90808
dim(testing)  #  18 90808
### 
training_zscore = training
testing_zscore = testing
alltss_zscore <- rbind(training_zscore,testing_zscore)
nEtrain = training_zscore
nEtest = testing_zscore

#########################################################
################# train-path-p  ###########
#########################################################
CAall <- nEtrain
dim(CAall) #  65 90770
### 
CA_100f <- data.frame(sites = colnames(CAall))
CA_100f <- left_join(CA_100f,bridge_matrix,by = 'sites')
CA_100f_gene <- unique(CA_100f$gene)
length(CA_100f_gene) # 51988
genes_used_path <- data.frame(gene = unique(KEGG_path$gene)) # 2727
inter_genes = intersect(genes_used_path$gene,CA_100f_gene) # 2673

### the site number of 'EGFR'
site_N = CA_100f[which(CA_100f$gene %in% inter_genes),'sites']#7653
used_gene = data.frame(sites = site_N)
used_gene= left_join(used_gene,bridge_matrix,by = 'sites')
dim(used_gene) #  7654    4

gene_table <- data.frame(CAall[,site_N])
rownames(gene_table) = rownames(CAall)
gene_table$sample = rownames(gene_table)
gene_table <- left_join(gene_table,lung_clinic[,c('sample','group','EGFR')],by='sample')
gene_table[is.na(gene_table)] = 0
Nsites = ncol(gene_table) - 3
gene_table2 <- gene_table[,c(1:Nsites,Nsites+3)]

### calculate the p-values for gene level (EGFR -vs nonEGFR)
alltss_zscore <- gene_table2[,1:Nsites]
group_train = data.frame(disease = gene_table2$EGFR)
test_gene <- cbind(group_train$disease,alltss_zscore)
dim(group_train)
tss_p_genes_EGFR <- p_cal(trainTransformed = alltss_zscore,
                          trainLabels = group_train$disease,
                          trainLabels_sig = 'EGFR',paired = FALSE) #"0.000687248266510784"     "0.994189112386891"            
dim(tss_p_genes_EGFR) #  7651    2
### add gene symbol to the tss_p_genes_EGFR
colnames(tss_p_genes_EGFR)[2] = 'sites'
tss_p_genes_EGFR <- left_join(tss_p_genes_EGFR,bridge_matrix,by='sites')
### sort by the p ### gene-p-value-train
tss_p_genes_EGFR <- tss_p_genes_EGFR[order(tss_p_genes_EGFR$tout),]
#write.csv(tss_p_genes_EGFR ,file = '～/Desktop/BDF/script/output/var_more/final/tss_p_genes_EGFR.csv')
##########################################
### unique genes for pathway calculate
tss_uniq_genes_EGFR <- tss_p_genes_EGFR[!duplicated(tss_p_genes_EGFR$gene),]
### choose genes in EGFR pathways
genes_used_path <- data.frame(gene = unique(KEGG_path$gene)) # 2727
inter_genes = data.frame(gene = intersect(genes_used_path$gene,CA_100f_gene)) # 2673
inter_genes <- left_join(inter_genes,tss_uniq_genes_EGFR,'gene')
inter_genes_p <- inter_genes[which(inter_genes$tout<0.05),] # 325
inter_genes_p <- inter_genes[which(inter_genes$tout<0.01),] # 65
write.csv(inter_genes,
          file = '/Users/pingyi/Desktop/BDF/data/ATAC/inter_genes_for_EGFR_model.csv')
