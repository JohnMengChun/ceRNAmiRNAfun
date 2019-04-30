#' sort the information to find out the triplets and plot the top miRNA expression level of modulated ceRNA interaction happened
#'
#' This function will sort the information to find out the triplets and plot the top miRNA expression level of modulated ceRNA interaction happened
#' @return plotp a plot format with miRNA expression in x-axis and the interaction ceRNA in y-axis.
#'
#' @param dictionary miRNA and the corresponding genes' combination in list format, with miRNA name in the first column and the corresponding genes in the second column 
#' @param mirna_sam expression data of miRNA in dataframe format, with miRNA's name in rows and sample name in columns
#' @param Result_TCGA_LUSC a list of dataframe formats with miRNA,candidate ceRNAs,and their locations and number of segmentation based on segcluster_peakmerge function. 
#' @param miRhubgeneoutput a dataframe formats with miRNA names with the number of bridging ceRNA triplets and the corresponding genes with the number of ceRNA triplets.
#' @param w window size of the each calculation of triplets,the default is 10.
#' @param N number of the samples,the default is 475.
#'
#' @examples
#' ## Use the internal dataset
#' 
#' data("dictionary", package = "ceRNAmiRNAfun", envir = environment())
#' data("mirna_sam", package = "ceRNAmiRNAfun", envir = environment())
#' data("Result_TCGA_LUSC", package = "ceRNAmiRNAfun", envir = environment())
#' data("miRhubgeneoutput", package = "ceRNAmiRNAfun", envir = environment())
#' 
#' ## use miRhubgeneplot to merge out the bridging miRNA and hub genes and then make the plot of the top miRNA expression level of modulated ceRNA interaction happened
#' miRhubgeneplot(dictionary,mirna_sam,Result_TCGA_LUSC,miRhubgeneoutput,w=10,N=475)
#'   
#' @import data.table
#' @import ggplot2
#' @export

miRhubgeneplot <- function(dictionary,mirna_sam,Result_TCGA_LUSC,miRhubgeneoutput,w,N){

all_putative_p <- unlist(lapply(dictionary[,2], function(x) unlist(c(x))))
length(unique(all_putative_p))   

sum(unlist(lapply(dictionary[,2], function(x) choose(length(unlist(c(x))),2) ))) 

ceRNA_pair_count_2 <- unlist(lapply(Result_TCGA_LUSC,function(x) dim(x)[1]))
rank_count_2 <- sort(ceRNA_pair_count_2,decreasing = TRUE) # get all

rank_miRNA_2=list()
rank_count_2 <- unique(rank_count_2)
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

mir=(as.character(miRhubgeneoutput$miRNA))[1]
w=10

N=475
a=Result_TCGA_LUSC_dataframe[Result_TCGA_LUSC_dataframe[,1]==mir ,2:4]


hubgene_2=matrix(unlist(sapply(Result_TCGA_LUSC_dataframe[,2], function(x) strsplit(as.character(x)," "))),ncol=2,byrow = TRUE)
hubgene_2=as.data.frame(hubgene_2)
colnames(hubgene_2)=c('ceRNA1','ceRNA2')
dim(hubgene_2) 
df_LUSC=data.frame("miRNA"=unlist(Result_TCGA_LUSC_dataframe[,1]), 'ceRNA1'=hubgene_2[,1],'ceRNA2'=hubgene_2[,2])
df_LUSC=sapply(df_LUSC, function(x) as.character(x))
gp2=apply(df_LUSC, 1, function(x) paste(sort(c(x[2],x[3]))[1],sort(c(x[2],x[3]))[2])  )
df_LUSC=data.frame(df_LUSC[,1],gp2,stringsAsFactors = FALSE)

overlap_ceRNA_pair=c()
for(i in unlist(summary_miR_2[,1])){
  tmp <- unlist(df_LUSC[df_LUSC[,1]==i, 2]) 
  overlap_ceRNA_pair=rbind(overlap_ceRNA_pair, c("miRNA"=i, 'ceRNA'=list(tmp)))
  
}



hubg_table_overlap <- sort(table(unlist(strsplit(unlist(overlap_ceRNA_pair[,2])," "))),decreasing = TRUE)
hubg_top <- names(hubg_table_overlap)
Num_of_including_hubg <- apply(sapply(hubg_top, grepl,unlist(a[,1])),1,sum)
Index_including_bubg <- which(Num_of_including_hubg >1)
a=a[Index_including_bubg,]


a_order <- a[order(unlist(lapply(a[,2], function(x) last(x[,1]))),decreasing = TRUE),]
start_ls <- lapply(a_order[,2], function(x) na.omit(x[,1]))
end_ls <- lapply(a_order[,2], function(x) na.omit(x[,2]))
start <- unlist(start_ls)
end <- unlist(end_ls)
location <- end_ls  
ystart=rep(0:(dim(a_order)[1]-1),unlist(lapply(location, length)))
yend=ystart+1

d=data.frame(x1=start, x2=end, y1=ystart, y2=yend, 
             ylab=unlist(rep(a_order[,1],unlist(lapply(location, length)))), 
             yloc=rep(seq(0.5,dim(a_order)[1],by=1),unlist(lapply(location, length))),
             Count=as.matrix(na.omit(unlist(a_order[,3]))))

x_range <- c(mean(as.numeric(mirna_sam[mir,order(mirna_sam[rownames(mirna_sam)==mir,])][1:w])),
             mean(as.numeric(mirna_sam[mir,order(mirna_sam[rownames(mirna_sam)==mir,])][N:(N-w+1)])))

p <- ggplot(data=d,aes(fill=Count^2))+ xlim(x_range)+
  geom_rect(mapping=aes(xmin=x1, xmax=x2, ymin=y1, ymax=y2), 
            colour="steelblue", alpha=0.9)+ 
  theme(axis.text.y = element_blank(),axis.ticks.y = element_blank())+ 
  scale_fill_gradient2(low="black", high="steelblue")


plotp=p+ scale_y_continuous(breaks = d$yloc, labels = d$ylab)+
  theme_bw()+
  xlab("miRNA")+
  theme(axis.text.x=element_text(size=15),axis.text.y = element_text(size=8), axis.title=element_text(size=15))

return(plotp)
}


