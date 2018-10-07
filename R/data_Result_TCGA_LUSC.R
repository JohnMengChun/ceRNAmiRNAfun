#' a list of dataframe formats with miRNA,
#' candidate ceRNAs,and their locations and  
#' number of segmentation based on 
#' segcluster_peakmerge function. 
#'
#' This dataframe is generated from \link{segcluster_peakmerge},
#' 
#'
#' @format A large dataframe with miRNA's name in 
#' the first column, candidate ceRNA names, the ceRNA's
#' locations, and the number of segmentation as order in
#' in the following columns.
#' 
#' 
#' @return list
#' @source TCGA data portal:
#'  \url{https://tcga-data.nci.nih.gov/docs/publications/tcga/}
"Result_TCGA_LUSC"
