#' evaluate the correlation coefficients between two genes to find out the pair ceRNAs of a miRNA.
#'
#' This function will evaluate the fraction of samples which their gene correlations are under the condition weighted by their correlation and the fractions that are not presnt up to a given position in all samples. And then, we sort the genes with increasing expression level of a miRNA, by permuting the miRNA expression for 1000 times, and computed S, the maximum deviation from zero of the fractions of samples which their gene correlations are under the condition and those are not, of the permuted data. The empirical P value is calculated relative to this null distribution. Based on the P value, we can find out the ceRNAs of a miRNA. 
#' @return numeric format with P value of each pair of genes with each miRNA in overall triplet.  
#'
#' @param triplet two genes and a miRNA combination.
#' @param cor_shreshold_rw the shreshold of correlation of coefficient between two genes.


#' @examples
#' ## Use the internal dataset
#' data("triplet", package = "anamiR", envir = environment())
#'
#' ## evaluate correlation coefficents between two genes.
#' Random_walk(triplet,cor_shreshold_rw=0.5)
#'
#' @export

Random_walk <- function(triplet,cor_shreshold_rw){
    triplet <- triplet[order(triplet$miRNA),]
    colnames(triplet) <- c("miRNA","corr")
    n=dim(triplet)[1]
    score <- c(0,rep(NA,n))
    Nr <- sum(abs(triplet[abs(triplet$corr) > cor_shreshold_rw,"corr"]))
    for(i in 1:n){
      
      if(abs(triplet$corr[i]) > cor_shreshold_rw){ 
        score[i+1] <- score[i]+abs(triplet$corr[i])/Nr
      }else{
        score[i+1] <- score[i]-1/sum(abs(triplet$corr) < cor_shreshold_rw)  
      }
    }
    EnScore <- max(abs(score))  
    
    
    # permutation
    score_per <- c(0,rep(NA,n))
    EnScore_per <- rep(0,1000) 
    for(sim in 1:1000){       
      
      triplet_ran <- unique.data.frame(triplet)
      triplet_ran[,1] <- sample(triplet[,1],n)
      triplet_ran <- triplet_ran[order(triplet_ran[,1]),]
      
      for(i in 1:n){
        
        if(abs(triplet_ran$corr[i]) > cor_shreshold_rw){    
          score_per[i+1] <- score_per[i]+abs(triplet_ran$corr[i])/Nr
          
        }else{
          score_per[i+1] <- score_per[i]-1/sum(abs(triplet_ran$corr) < cor_shreshold_rw)   
        } 
      }
      EnScore_per[sim] <- max(abs(score_per))
    }
    p.value <- min(sum(EnScore_per<=EnScore)/1000,sum(EnScore_per>=EnScore)/1000) 
    return(p.value)
  }

