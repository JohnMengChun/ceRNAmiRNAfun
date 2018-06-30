#' calculate correlation coefficients with ceRNAs of a miRNA by sliding windows. 
#'
#' This function will calculate correlation coefficients within each window, a sliding mover that contains putative ceRNA triplets composed of a miRNA and several genes, after determing a window size and the number of windows. And the corresponding miRNA expression of each window was the average expression of samples within different windows. 
#' @return dataframe format with correlation coefficients within each window for each corresponding miRNA average expression. 
#'
#' @param w window size of the each calculation of triplets,the default is 105.
#' @param miRNA_total miRNA dataset in a dataframe format


#' @examples
#' ## Use the internal dataset
#' data("mirna_sam", package = "anamiR", envir = environment())
#' data("gene_sam", package = "anamiR", envir = environment())
#'
#'
#' ## use sliding window to calculate correlations
#' Datalist_sliding_w(w=105,miRNA_total=mirna_total)
#'   
#'
#' @import doParallel
#' @export


Datalist_sliding_w <- function(w,miRNA_total){
  w=105
  parallel_d <- foreach(mir=miRNA_total, .export = c('dictionary','mirna_sam', 'gene_sam'))  %dopar%  {
    print(which(mir==miRNA_total))
    Gene <- as.character(data.frame(dictionary[dictionary[,1]==mir,][[2]])[,1])
    Gene <- gsub("-",".",Gene)
    Gene <- intersect(Gene,rownames(gene_sam))
    gene_mir <- data.frame("miRNA"=t(mirna_sam[rownames(mirna_sam)==mir,]), t(gene_sam[intersect(Gene,rownames(gene_sam)),]))
    colnames(gene_mir)=c("miRNA",colnames(gene_mir)[-1])
    
    N=dim(gene_mir)[1]
    gene_pair <- combn(Gene,2)
    total_pairs <- choose(length(Gene),2)
    cor_all <- matrix(0,N-w+1,total_pairs)
    for(pair in 1:total_pairs){
      print(pair)
      
      r=gene_pair[1,pair]
      s=gene_pair[2,pair]
      
      y <- gene_mir[,c("miRNA",r,s)]
      y <- y[order(y$miRNA),]
      corr <- NULL
      miRNA <- NULL
      for(i in 1:(N-w+1)){  
        
        rand_dist <- data.frame(y[i:(i+w-1),])
        corr[i] <- cor(rand_dist[,2],rand_dist[,3])
        miRNA[i] <- mean(rand_dist$miRNA)
        
      }
      data=data.frame(cbind(miRNA,corr))
      data <- data[order(data$miRNA),]
      cor_all[,pair] <- data$corr
    }
    
    triplet = data.frame("miRNA"=data$miRNA,"corr"=cor_all)
    triplet
  }
  
  return(parallel_d)
  
}
