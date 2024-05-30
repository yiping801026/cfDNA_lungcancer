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
### import the in-house function
source('~/Desktop/BDF/script/BD_lung/used_function.R')
source('/Users/pingyi/Desktop/BDF/script/BD_lung/function_path.R')

################## normalize with CA non-CA model mean/sd#################
load(file = '~/Desktop/BDF/script/out_data/CAnonCA200f') 
### use training_range_site to normalize all data (train && validate)
training_range_siteCA = training_range_site
lung_clinic_withvar <- read.csv("~/Desktop/BDF/script/out_data/lung_clinic_withvar.csv", row.names=1)
bridge_matrix <- read.csv("~/Desktop/BDF/script/out_data/bridge_matrix.csv")
### Import the pathway information
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

######################################################################
################# calculate the gene p-values in training  ###########
######################################################################
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
### unique genes for pathway calculate
tss_uniq_genes_EGFR <- tss_p_genes_EGFR[!duplicated(tss_p_genes_EGFR$gene),]
### choose genes in EGFR pathways
genes_used_path <- data.frame(gene = unique(KEGG_path$gene)) # 2727
inter_genes = data.frame(gene = intersect(genes_used_path$gene,CA_100f_gene)) # 2673
inter_genes <- left_join(inter_genes,tss_uniq_genes_EGFR,'gene')
#inter_genes_p <- inter_genes[which(inter_genes$tout<0.05),] # 325
inter_genes_p <- inter_genes[which(inter_genes$tout<0.01),] # 65
write.csv(inter_genes,
          file = '/Users/pingyi/Desktop/BDF/data/ATAC/inter_genes_for_EGFR_model.csv')

#########################################################################
######### build model based on the p<0.01 genes - 65 totally  ###########
#########################################################################
train_s = clinic[inTrain,'sample']
test_s = clinic[-inTrain,'sample']
### build the model based on features:rank_sites
train_topbins<-as.data.frame(nEtrain[,inter_genes_p$sites])
dim(train_topbins) 
test_topbins<-as.data.frame(nEtest[,inter_genes_p$sites])
dim(test_topbins) 
train_topbins$group <- factor(trainLabels,levels = c('EGFR','zero'))
test_topbins$group <- factor(testLabels,levels = c('EGFR','zero'))

set.seed(123) 
ctrl = trainControl(method = "repeatedcv", number = 10, repeats = 3,
                    allowParallel = TRUE,classProbs = TRUE,savePredictions = TRUE,search='grid',summaryFunction = twoClassSummary)
set.seed(123) 
ctrl$sampling <- "up"
tunegrid <- expand.grid(mtry = c(5))
rf_default1 <- train(group~., 
                     data=train_topbins, 
                     method='rf', 
                     #metric='Accuracy',  #Metric compare model is Accuracy
                     metric = "ROC",
                     tuneGrid=tunegrid,
                     ntree=c(2500),
                     trControl=ctrl)
rf_default1

pred <- predict(rf_default1,newdata = test_topbins,type = "prob")
pred$sample = rownames(pred)
lung_clinic_withvar$group = lung_clinic_withvar$EGFR
pred <- left_join(pred,lung_clinic_withvar[,c('sample','group')],by='sample')
pred$res = 1
pred[which(pred$EGFR<cutoff),'res'] = 0 

pred$true = 1
#pred$true[1:25] = 0
pred[which(pred$group != 'EGFR'),'true'] = 0 

pred$equ = 1
pred[which(pred$res != pred$true),'equ'] = 0 

roc_data <- roc(pred$true, pred[,1],type="prob",ci = T)
roc_data


#########################################################################
######### plot the heatmap based on p<0.01 genes - 65 totally  ##########
#########################################################################

used_100f <- rbind(train_topbins, test_topbins)
### PCA for CA vs nonCA
samples <- data.frame(sample = rownames(used_100f))
samples <- left_join(samples,lung_clinic_withvar[,c('sample','EGFR')],by = 'sample')
library(pheatmap)
heat_plot <- used_100f[,1:65]

range(heat_plot)

#heat_plot <- rbind(heat_plot,raw_plot_data[model_clinic[which(model_clinic$stage_less == '0'),'sample'],]) #(263)
heat_plot <- rbind(used_100f[samples[which(samples$EGFR == 'zero'),'sample'],],
                   used_100f[samples[which(samples$EGFR == 'EGFR'),'sample'],]) # 100,32s

dim(heat_plot) 
#annotation_row <- data.frame( group = c(rep('Healthy',79),rep('Cancer',94)))
annotation_row <- data.frame(
  group = c(rep('EGFR-',36),rep('EGFR+',29)))
row.names(annotation_row) <- rownames(heat_plot)
annotation_row$sample = rownames(heat_plot)

annotation_row <- left_join(annotation_row,samples,'sample')
colnames(annotation_row)
annotation_row <- data.frame(group = (annotation_row[,c("group")]))
row.names(annotation_row) <- rownames(heat_plot)
colnames(annotation_row) <- c('group')

ann_colors = list(
  
  group = c(`EGFR-`="#ADDB88" ,`EGFR+` = '#B4B4D5'))


p1 <- pheatmap(heat_plot[,1:65],
               #legend_breaks = 0:,
               cluster_rows = F,
               cluster_cols = T,
               #scale = 'column',
               #scale = 'row',
               scale = 'none',
               annotation_row = annotation_row,
               annotation_colors = ann_colors,
               show_colnames = F,
               show_rownames = F,
               cell_width = 0.1 ,
               cell_height = 0.01,
               breaks = -4:5,
               color = colorRampPalette(c("darkblue", "white", "tomato2"))(8),
               
               #color = colorRampPalette(c("darkblue", "wheat3", "tomato2"))(8),
               #color = colorRampPalette(c("darkblue", "white", "red4"))(100),
               #annotation_col=annotation_col
               
               #gaps_col=col_gap,
               gaps_row = c(36),
               border_color = 'white'
               #border = TRUE
)

p1

#########################################################################
################# calculate the pathway p-values in training  ###########
#########################################################################

EGFR_path_used <- unique(KEGG_path$term) # 619
path_p_EGFR <- data.frame(EGFR = gene_table2$EGFR)
#save(gene_table2,file = '～/Desktop/BDF/script/output/var_more/final/gene_table2')
site_used_path <- colnames(gene_table2)[1:7653]
calculate_matrix_path <- rbind(nEtrain,nEtest)
dim(calculate_matrix_path)
### calculate path score for train
path_p_EGFR <- path_score_calculate(new_matrix = calculate_matrix_path,
                                    site_used_path,
                                    KEGG_path=KEGG_path,
                                    lung_clinic_withvar)

path_names <- colnames(path_p_EGFR)[2:ncol(path_p_EGFR)]
colnames(path_p_EGFR)[2:ncol(path_p_EGFR)] <- path_names

### calculate the p-values for pathway level (EGFR -vs nonEGFR)
alltss_zscore <- path_p_EGFR[,2:ncol(path_p_EGFR)]
dim(alltss_zscore)
group_train = data.frame(disease = path_p_EGFR$EGFR)
dim(group_train)
tss_p_path_EGFR <- p_cal(trainTransformed = alltss_zscore,
                         trainLabels = group_train$disease,
                         trainLabels_sig = 'EGFR',paired = FALSE)           
dim(tss_p_path_EGFR) #   619   2
tss_p_path_EGFR$path_name <- as.character(EGFR_path_used)
tss_p_path_EGFR <- tss_p_path_EGFR[order(tss_p_path_EGFR$tout),]

tss_p_path_EGFR <- tss_p_path_EGFR[order(tss_p_path_EGFR$tout),]
write.csv(tss_p_path_EGFR,'/Users/pingyi/Desktop/BDF/script/out_data/tss_p_path_EGFR.csv')


#########################################################################
####################### plot boxplot for pathways  ######################
#########################################################################

path_names <- colnames(path_p_EGFR)[2:ncol(path_p_EGFR)]
for (i in 1:length(path_names)){
  #i = 2
  new_name <-  sub("KEGG_MEDICUS_",'',path_names[i])
  new_name <- sub("_SIGNALING_PATHWAY",'',new_name)
  new_name <- gsub("_", " ", new_name)
  path_names[i] <- new_name
}


colnames(path_p_EGFR)[2:ncol(path_p_EGFR)] <- path_names
long_df <- pivot_longer(path_p_EGFR, cols = -EGFR, names_to = "variable", values_to = "value")
#all_path_genes_tss2 <- all_path_genes_tss[,which(grepl("EGFR",colnames(all_path_genes_tss)))] # 26paths+group(EGFR)

long_df2 <- long_df[which(grepl("EGFR",long_df$variable)),]

p <- ggboxplot(long_df2, x = "variable", y = "value",
               color = "EGFR", palette = "jco",
               add = "jitter")

p + stat_compare_means(aes(group = EGFR),label =  "p.signif")+ 
  ggplot2::theme(axis.text = element_text( color = "black", size = 10),
                 plot.title = element_text(size = 5), 
                 plot.margin = unit(c(4,12, 4, 12), "cm"))+coord_flip()













