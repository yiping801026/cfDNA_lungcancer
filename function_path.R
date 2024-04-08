
#################################################################
######## NO.1 function for KEGG path matrix calculate############
#################################################################

# eg. path_score_calculate(new_matrix = nEtest,site_used_path,KEGG_path=KEGG_path,lung_clinic_withvar)
path_score_calculate <- function(new_matrix,site_used_path,KEGG_path,lung_clinic_withvar){
  
  ### choose sites
  new_matrix <- new_matrix[,site_used_path]
  dim(new_matrix)
  ### add group 
  new_matrix$sample = rownames(new_matrix)
  new_matrix <- left_join(new_matrix,lung_clinic_withvar[,c('sample','EGFR')],by='sample')
  Nsites = ncol(new_matrix) - 2
  ### for all path to calculate the path_feature_table
  path_p_EGFR_new <- data.frame(EGFR = new_matrix$EGFR)
  KEGG_path_EGFR <- KEGG_path
  EGFR_path_used <- unique(KEGG_path$term) # 619
  ### for all path in EGFR_path_used
  for (i in 1:length(EGFR_path_used)){
    #i = 1
    genes_one_path <- KEGG_path_EGFR[which(KEGG_path_EGFR$term == EGFR_path_used[i]),]
    genes_one_path <- left_join(genes_one_path,tss_uniq_genes_EGFR,by='gene')
    genes_one_path <- genes_one_path[!is.na(genes_one_path$sites),]
    gene_table3 <- data.frame(t(new_matrix[,1:Nsites]))
    gene_table3$sites <- colnames(new_matrix[1:Nsites])
    genes_one_path <- left_join(genes_one_path,gene_table3[1:Nsites,],by='sites')
    #genes_one_path[nrow(genes_one_path)+1,6:ncol(genes_one_path)] = gene_table2$EGFR
    genes_one_path <- genes_one_path[,7:ncol(genes_one_path)]
    genes_one_path_t<- data.frame(t(genes_one_path))
    genes_one_path_t$mean <- apply(genes_one_path_t,1,mean)
    genes_one_path_t$EGFR <- new_matrix$EGFR
    path_p_EGFR_new[,i+1] <-  genes_one_path_t$mean
  }
  ### change path names
  colnames(path_p_EGFR_new)[2:ncol(path_p_EGFR_new)] <- as.character(EGFR_path_used)
  rownames(path_p_EGFR_new) = new_matrix$sample
  return(path_p_EGFR_new)

}


#########################################################
######## NO.2 function for  path name change ############
#########################################################

#path_names <- colnames(path_p_EGFR_new)[2:ncol(path_p_EGFR_new)]
### eg. path_names <- KEGG_less_name(path_names)
KEGG_less_name <- function(path_names){
  for (i in 1:length(path_names)){
    new_name <-  sub("KEGG_MEDICUS_",'',path_names[i])
    new_name <- sub("_SIGNALING_PATHWAY",'',new_name)
    new_name <- gsub("_", " ", new_name)
    path_names[i] <- new_name}
  return(path_names)
  }

#colnames(path_p_EGFR_new)[2:ncol(path_p_EGFR_new)] <- path_names


#########################################################
######## NO.3 function for rank gap calculate ###########
#########################################################
#eg. path_gap <- rank_gap_groups(clinic_rank,path_p_EGFR)
rank_gap_groups <- function(clinic_rank,path_p_EGFR){
  ### rank matrix for train_topbins
  train_topbins_rank <- data.frame(apply(path_p_EGFR[,-1],2,function(x) return(rank(x))))
  train_topbins_rank$sample = rownames(train_topbins_rank)
  ### sum of several groups
  train_topbins_rank <- left_join(train_topbins_rank,clinic_rank,by='sample')
  colnames(train_topbins_rank)[ncol(train_topbins_rank)]
  Nsites = ncol(train_topbins_rank)-2
  Nsites
  group_mean <- function(value_col,stage_less){
    matrix_col <- data.frame(stage_less = stage_less,
                             value = value_col)
    result_matrix <- data.frame(aggregate(value ~ stage_less,data = matrix_col, FUN = mean))
    return(result_matrix$value)
  }
  
  result_matrix_mean <- data.frame(site=colnames(train_topbins_rank)[1:Nsites])
  result_matrix_mean$zero = 1
  result_matrix_mean$nE = 1
  
  for (i in 1:Nsites) {
    result_matrix_mean[i,2:3] <- group_mean(train_topbins_rank[,i],stage_less = train_topbins_rank$EGFR)
  }
  
  ### 
  rownames(result_matrix_mean) <- result_matrix_mean$site
  result_matrix_mean$gap = result_matrix_mean[,3]-result_matrix_mean[,2]
  
  ### sort by absolute gap && choose top 200
  NV_feature_gap <- result_matrix_mean[order(-abs(result_matrix_mean$gap)),]
  colnames(NV_feature_gap)[1] <- 'path_name'
  return(NV_feature_gap)
}




