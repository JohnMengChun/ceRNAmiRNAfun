#' calculate correlation coefficients with ceRNAs of a miRNA by sliding windows.
#'
#' This function will calculate correlation coefficients within each window, a sliding mover that contains putative ceRNA triplets composed of a miRNA and several genes, after determining a window size and the number of windows. And the corresponding miRNA expression of each window was the average expression of samples within different windows.
#' @return Realdata a list of dataframe formats with correlation coefficients of genes within each window for each corresponding miRNA average expression.
#'
#' @param w window size of the each calculation of triplets,the default is 10.
#' @param dictionary miRNA and the corresponding genes' combinatino in list format, with miRNA name in the first column and the corresponding genes in the second column
#' @param miRNA_total miRNA list in a vector format, the order of miRNA_total should be the same as the first column of dictionary
#' @param mirna_sam expression data of miRNA in dataframe format, with miRNA's name in rows and sample name in columns
#' @param gene_sam expression data of gene in dataframe format, with gene's name in rows and sample name in columns
#'
#' @examples
#' ## Use the internal dataset
#'
#' data("dictionary", package = "ceRNAmiRNAfun", envir = environment())
#' data("miRNA_total", package = "ceRNAmiRNAfun", envir = environment())
#' data("mirna_sam", package = "ceRNAmiRNAfun", envir = environment())
#' data("gene_sam", package = "ceRNAmiRNAfun", envir = environment())
#'
#' ## use sliding window to calculate correlations
#' Datalist_sliding_w(w=10,dictionary,miRNA_total,mirna_sam,gene_sam)
#'
#'
#' @import doParallel
#' @import foreach
#' @import data.table
#' @export

Datalist_sliding_w <- function(w,dictionary,miRNA_total,mirna_sam,gene_sam){
  datap <- foreach(mir=miRNA_total)  %dopar%  {
    print(which(mir==miRNA_total))
    Gene <- as.character(data.frame(dictionary[dictionary[,1]==mir,][[2]])[,1])
    Gene <- gsub("-",".",Gene)
    Gene <- intersect(Gene,rownames(gene_sam))
    Gene <- unlist(lapply(Gene, function(k) if(identical(k, character(0)))NULL else k))
    gene_mir <- data.frame("miRNA"=t(mirna_sam[rownames(mirna_sam)==mir,]), t(gene_sam[Gene,]))
  }


  dataps=which(lengths(datap)>=3)
  namep=gsub("\\.",replacement="-",lapply(datap[dataps],colnames))

  miRSPLIT=strsplit(miRNA_total[dataps],",")

  nameuse=mapply(grepl,miRSPLIT,namep)
  miRSPLIT1=unlist(miRSPLIT[nameuse])
  geneuse=lapply(datap[dataps][nameuse],colnames)
  lapply(geneuse,function(t){geneuse[1]=NULL;t})
  datause=(datap[dataps][nameuse])




  parallel_d=foreach(i=c(1:sum(nameuse))) %dopar% {
    datapdf=data.frame(datause[i])
    colnames(datapdf)=c("miRNA",colnames(datapdf)[-1])
    N=dim(datapdf)[1]
    genef=colnames(datapdf)[-1]
    gene_pair <- combn(genef,2)
    total_pairs <- choose(length(genef),2)
    cor_all <- matrix(0,N-w+1,total_pairs)
    for(pair in 1:total_pairs){
      print(pair)

      r=gene_pair[1,pair]
      s=gene_pair[2,pair]

      y <- datapdf[,c("miRNA",r,s)]
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

