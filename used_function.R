### import library
library(readxl)
library(dplyr)
library(factoextra) 
library(ggplot2)
library(pROC)
library(stats)
library(caret)

### tssMatrix_rmNA2 should be the result of function 'pre_treat_rm0pNA' 
### && model_clinic$platform includes the inform of group to be coloured in the plot
### path be the output dir '/Users/pingyi/Desktop/XJ/XJ_paper/output/figures/batch/tss_raw_all/'
#platform = model_clinic$platform
#path = '/Users/pingyi/Desktop/XJ/XJ_paper/output/figures/batch/tss_raw_all/'
PCA_plot_10PCs <- function(tssMatrix_rmNA2 ,platform,path,xlimn,ylimn){
  ### PCA plots of top10 PCs to check the batch effect
  numbers <- 1:10
  all_combinations <- expand.grid(numbers, numbers)
  filtered_combinations <- all_combinations[all_combinations$Var1 > all_combinations$Var2, ]
  result <- as.matrix(filtered_combinations)
  rownames(result) <- 1:nrow(result)
  ### calculate the PCA results
  combat_tss_rmt <- data.frame(t(tssMatrix_rmNA2))
  pca_result  <- prcomp(combat_tss_rmt)
  
  for (i in 1:nrow(result)) {
    #i =1
    p <- fviz_pca_ind(pca_result ,axes = c(result[i,1], result[i,2]),
                      label=" ", habillage=platform,
                      addEllipses=T, ellipse.level=0.5,
                      palette = c('brown3','#FFB90F','red'))+
      ggtitle(" ") +
      ggplot2::theme_bw() +
      ggplot2::ylim(-ylimn,ylimn)+
      ggplot2::xlim(-xlimn,xlimn)+
      coord_fixed(ratio = 1)+ 
      ggplot2::theme(axis.text = element_text( color = "black", size = 15),
                     text = element_text(size = 25),  # 调整文字大小
                     plot.title = element_text(size = 20),  # 调整标题文字大小
                     plot.margin = unit(c(2, 2, 2, 2), "cm"),  # 调整边距
                     panel.grid.major=ggplot2::element_line(colour=NA),panel.grid.minor = ggplot2::element_blank())
    p
    #p2 <- p + ggplot2::ggtitle(" ") + ggplot2::theme_bw() + ggplot2::theme(panel.grid=ggplot2::element_blank())
    #p2
    pdf(file=paste0(path,result[i,1],'vs',result[i,2],'.pdf'),bg = "white")
    print(p)
    dev.off()
    
  }
  
}

### calculate the median_z-score for the matrix(ignore outliers)
### eg. try <- data.frame(t(apply(data.frame(tssMatrix_rmNA2), 1, custom_z_score)))
custom_z_score <- function(x) {
  median_value <- median(x)
  mad_value <- mad(x)
  z_scores <- (x - median_value) / mad_value
  return(z_scores)
}


### This function is used to treat the new_matrix the same as the training set
### eg. newmatrix_z_norm(new_matrix=testing,training_range_site)
newmatrix_z_norm <- function(new_matrix,training_range_site){
  
  #new_matrix = testing
  Nnew_samples <- dim(new_matrix)[1]
  Nmedian = Nnew_samples + 1 
  Nmad = Nnew_samples + 2
  #range(new_matrix)
  #new_matrix <- new_matrix[,-1]
  new_matrix <- new_matrix[,rownames(training_range_site)]
  new_matrix <- rbind(new_matrix,data.frame(t(training_range_site[,c('median','mad')])))
  new_matrix_norm <- data.frame(apply(data.frame(new_matrix), 2, function(x)return((x[1:Nnew_samples]-x[Nmedian])/x[Nmad]) ))
  return(new_matrix_norm)
}

### calculate the row or col range, median, mean data in the matrix
### eg. try_range <- range_cal(matrix = try,roworcol = 1)
range_cal <- function(matrix,roworcol){
  range_matrix <- data.frame(t(apply(data.frame(matrix), roworcol, function(x)return(range(x)))))
  colnames(range_matrix) <- c('min','max')
  range_matrix$median <- apply(data.frame(matrix), roworcol, function(x)return(median(x)))
  range_matrix$mean <- apply(data.frame(matrix), roworcol, function(x)return(mean(x)))
  range_matrix$mad <- apply(data.frame(matrix), roworcol, function(x)return(mad(x)))
  
  return(range_matrix)
}


### calculate the p-values between groups
### eg. p_cal(trainTransformed = trainTransformed,trainLabels,trainLabels_sig = 'CA')
p_cal <- function(trainTransformed,trainLabels,trainLabels_sig,paired = FALSE){
  ### use p-value to make features smaller
  binNum<-dim(trainTransformed)[2]
  tout<-rep(1,binNum)
  CA<-trainLabels==trainLabels_sig
  if(paired == FALSE ){
    print('Use the t.test')
    for(i in 1:binNum){tout[i]<-t.test(trainTransformed[CA,i], trainTransformed[!CA,i])$p.value}
    } else {
      print('Use the paired t.test')
      for(i in 1:binNum){tout[i]<-t.test(trainTransformed[CA,i], trainTransformed[!CA,i],paired = TRUE)$p.value}}
  
  tout_withID<-data.frame(cbind(tout,colnames(trainTransformed)))
  tout_withID[,1] <- as.numeric(tout_withID[,1])
  #tout_sort<-tout_withID[order(tout_withID[,1]),]
  tout_sort<-tout_withID
  print(c('The range of the p-value is: ',range(tout_withID[,1]))) # 3.443190e-24 9.999481e-01
  return(tout_sort)
}

### calculate the fold_change between groups
### eg. FD_cal(trainTransformed = BRCA_norm_70counts,trainLabels= Basaltype,trainLabels_sig = 'Basal',paired = FALSE)
FD_cal <- function(trainTransformed,trainLabels,trainLabels_sig,paired = FALSE){
  #trainTransformed = BRCA_norm_70counts
  #trainLabels= Basaltype
  #trainLabels_sig = 'Basal'
  #paired = FALSE
  ### use p-value to make features smaller
  binNum<-dim(trainTransformed)[2]
  tout<-rep(1,binNum)
  CA<-trainLabels==trainLabels_sig
  
  for(i in 1:binNum){
    
    tout[i]<- mean(trainTransformed[CA,i])/mean(trainTransformed[!CA,i])
      
    }

  tout_withID<-data.frame(cbind(tout,colnames(trainTransformed)))
  tout_withID[,1] <- as.numeric(tout_withID[,1])
  #tout_sort<-tout_withID[order(tout_withID[,1]),]
  tout_sort<-tout_withID
  print(c('The range of the fold-change is: ',range(tout_withID[,1]))) # 3.443190e-24 9.999481e-01
  return(tout_sort)
}







### columns of this input matrix are the samples with rows be the features
pre_treat_rm0pNA <- function(matrix){
  #### change -1 into 0 
  matrix[matrix == -1] = 0
  #### remove all 0s
  print('columns of this input matrix are the samples with rows be the features')
  print(c('The structure of the raw matrix is:',dim(matrix)))
  tss_table_raw <- data.frame(t(matrix))
  tssMatrix_rowSum<-rowSums(tss_table_raw,na.rm=T)
  tssMatrix_colSum<- colSums(tss_table_raw,na.rm=T)
  tssMatrix_rm0<-tss_table_raw[tssMatrix_rowSum!=0,]
  tssMatrix_rm0<-tssMatrix_rm0[,tssMatrix_colSum!=0]
  
  print(c('The structure of the matrix after remove all NA & all zero is:',dim(tssMatrix_rm0) ))
  #### remove rows or columns with many NA and 0s 
  zero_matrix<-tssMatrix_rm0==0
  zero_colSum<-colSums(zero_matrix,na.rm=T)
  zero_rowSum<-rowSums(zero_matrix,na.rm=T)
  na_matrix<-is.na(tssMatrix_rm0)
  tss_naSum<-rowSums(na_matrix)
  bin_naSum<-colSums(na_matrix)
  
  tss_sampleNum<-dim(tssMatrix_rm0)[2]
  tss_binNum<-dim(tssMatrix_rm0)[1]
  #tssMatrix_rmNA<-t(tssMatrix_rm0[(tss_naSum+zero_rowSum)<(tss_sampleNum/2),])#
  tssMatrix_rmNA<-data.frame(tssMatrix_rm0[,(bin_naSum+zero_colSum)<(tss_binNum/2)])#
  
  print(c('The structure of the matrix after remove many NA & all zero is:',dim(tssMatrix_rmNA) ))
  print(c('The range of the matrix after remove many NA & all zero is:',range(tssMatrix_rmNA) ))
  
  return(tssMatrix_rmNA)
}

nearZeroVar_for_matrix<- function(tssMatrix_rmNA){
  
  nzv <- nearZeroVar(tssMatrix_rmNA, freqCut=80/20, uniqueCut=30)
  if(length(nzv)!=0)
  {
    tssMatrix_rmLowVar <- tssMatrix_rmNA[, -nzv]
  } else {
    tssMatrix_rmLowVar <- tssMatrix_rmNA
  }
  print(dim(tssMatrix_rmLowVar)) #263 90819
  
  matrix <-  data.frame(tssMatrix_rmNA)
  return(matrix)
}


                            






