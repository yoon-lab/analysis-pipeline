#Function for matching clusters to cell type 
#cellTypes is a vector of factor that have cell type information
#clusters is a vector of number that have cluster information

cluster_assign <- function(cellTypes,clusters){
  #Measure F1 Score matrix
  tbl_celltype <- table(cellTypes)
  tbl_cluster <- table(clusters)
  F1_mat <- matrix(NA, nrow=length(tbl_cluster),ncol=length(tbl_celltype))
  
  for (i in 1:length(tbl_cluster)){
    for (j in 1:length(tbl_celltype)){
      
      true_positives <- sum(clusters == i &
                              cellTypes == levels(cellTypes)[j])
      
      detected <- sum(clusters == i)
      truth <- sum(cellTypes == levels(cellTypes)[j])
      #calculate precision, recall and F1 Score
      
      precision_ij <- true_positives / detected
      recall_ij <- true_positives / truth
      F1_ij <- 2 * (precision_ij * recall_ij) / (precision_ij + recall_ij)
      
      if (F1_ij == "NaN") F1_ij <- 0
      

      F1_mat[i, j] <- F1_ij
    }
  }
  
rownames(F1_mat) <- names(tbl_cluster)
colnames(F1_mat) <- names(tbl_celltype)
  
  # Get cluster-cell type table 
  labels_matched <- rep(NA, nrow(F1_mat))
  for(i in 1:nrow(F1_mat)){
    labels_matched[i] <- names(which.max(F1_mat[i,]))
  }
  
  
  return(list("F1_mat" = F1_mat,
              "labels_matched" = labels_matched))
}