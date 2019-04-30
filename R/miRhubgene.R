#' sort the information to find out the triplets.
#'
#' This function will sort the information of miRNA and ceRNA to find out triplets.
#' @return outputm a dataframe formats with miRNA names with the number of bridging ceRNA triplets and the corresponding genes with the number of ceRNA triplets.
#'
#' @param dictionary miRNA and the corresponding genes' combination in list format, with miRNA name in the first column and the corresponding genes in the second column 
#' @param Result_TCGA_LUSC a list of dataframe formats with miRNA,candidate ceRNAs,and their locations and number of segmentation based on segcluster_peakmerge function. 
#'
#' @examples
#' ## Use the internal dataset
#' 
#' data("dictionary", package = "ceRNAmiRNAfun", envir = environment())
#' data("Result_TCGA_LUSC", package = "ceRNAmiRNAfun", envir = environment())
#' 
#' ## use miRhubgene to merge out the bridging miRNA and hub genes
#' miRhubgene(dictionary,Result_TCGA_LUSC)
#'   
#'
#' @export

miRhubgene <- function(dictionary,Result_TCGA_LUSC){
all_putative_p <- unlist(lapply(dictionary[,2], function(x) unlist(c(x))))
length(unique(all_putative_p))   

sum(unlist(lapply(dictionary[,2], function(x) choose(length(unlist(c(x))),2) ))) 

ceRNA_pair_count_2 <- unlist(lapply(Result_TCGA_LUSC,function(x) dim(x)[1]))
rank_count_2 <- sort(ceRNA_pair_count_2,decreasing = TRUE) 

rank_miRNA_2=list()
rank_count_2 <- unique(rank_count_2)#oringinal:rank_count_2 <- unique(rank_count)
for(i in 1:length(rank_count_2)){
  rank_miRNA_2[[i]] <- which(lapply(Result_TCGA_LUSC,function(x) dim(x)[1]==rank_count_2[i])==TRUE) 
}
rank_count_2 <- sort(ceRNA_pair_count_2,decreasing = TRUE)
rank_count_2 <- unique(rank_count_2)
summary_miR_2=c()
for(i in 1:length(rank_count_2)){
  rank_miRNA <- which(lapply(Result_TCGA_LUSC,function(x) dim(x)[1]==rank_count_2[i])==TRUE) 
  miR_name <- unlist(lapply(Result_TCGA_LUSC[rank_miRNA], function(x) unique(x[,1])))
  summary_miR_2 <- rbind(summary_miR_2,c(list(miR_name),rank_count_2[i]))
}

sum(unlist(summary_miR_2[,2]))/dim(summary_miR_2)[1] 



rank_count_2 <- sort(ceRNA_pair_count_2,decreasing = TRUE) 
rank_count_2 <- unique(rank_count_2)
hubgene_name=list()
for(i in 1:length(rank_count_2)){
  rank_miRNA <- which(lapply(Result_TCGA_LUSC,function(x) dim(x)[1]==rank_count_2[i])==TRUE)  
  tmp=c()
  for(j in rank_miRNA){
    a <- names(sort(table(unlist(strsplit(unlist(Result_TCGA_LUSC[[j]][,2]), " "))),decreasing=TRUE)[1])
    count <- sort(table(unlist(strsplit(unlist(Result_TCGA_LUSC[[j]][,2]), " "))),decreasing=TRUE)[1]
    
    tmp=rbind(tmp,c(a,count))
  }
  hubgene_name[[i]] <- tmp
}
matrix(unlist(hubgene_name),ncol=2,byrow=TRUE)
hubg_onlyname_2 <- unlist(sapply(hubgene_name, function(x) x[,1]))
hubmiR=matrix(0,length(unlist(hubgene_name))/2,2)
hubmiR[,1]=unlist(hubgene_name)[which(nchar(unlist(hubgene_name))>1)]
hubmiR[,2]=unlist(hubgene_name)[which(nchar(unlist(hubgene_name))==1)]




Result_TCGA_LUSC_dataframe=c()     
for(i in 1:length(Result_TCGA_LUSC)){
  Result_TCGA_LUSC_dataframe=rbind(Result_TCGA_LUSC_dataframe, Result_TCGA_LUSC[[i]])
}

hubgene_all_TCGA_LUSC <- sort(table(unlist(strsplit(as.character(Result_TCGA_LUSC_dataframe[,2])," "))), decreasing = TRUE)
sum(hubgene_all_TCGA_LUSC)/length(hubgene_all_TCGA_LUSC)  
length(hubgene_all_TCGA_LUSC)  


hg_2 <- names(hubgene_all_TCGA_LUSC)
ceRNA_acrossmiR <- list()
for(i in 1:length(hg_2)){
  ceRNA_acrossmiR[[i]] <- unlist(lapply(Result_TCGA_LUSC, function(x){
    if(class(x)=="matrix"){
      tmp=sum(unlist(strsplit(unlist(x[,2]),' ')) %in% hg_2[i])
      return(tmp)
    }else return(NA)
  }))
}

unlist(lapply(ceRNA_acrossmiR,function(x) sum(x!=0, na.rm = TRUE)))  


miroutputm=matrix(0,length(unlist(summary_miR_2[,1])),2)
j=1
for (i in 1:4)
{
miroutputm[c(j:(j-1+length(summary_miR_2[,1][i][[1]]))),2]=summary_miR_2[i,2][[1]]
for (p in (j:(j-1+length(summary_miR_2[,1][i][[1]])))){
miroutputm[p,1]=unlist(summary_miR_2[,1])[p]
}
j=length(summary_miR_2[,1][i][[1]])+j
}

cn=min(length(miroutputm[,1]),length(hubmiR[,1]))
outputm=cbind(miroutputm[(1:cn),],hubmiR[(1:cn),])
colnames(outputm)=c("miRNA","bridging ceRNA triplets","Gene","ceRNA triplets")
outputm=as.data.frame(outputm)
return(outputm)
}



