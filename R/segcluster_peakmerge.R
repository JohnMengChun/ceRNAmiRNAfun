#' do the segment clustering and peak merging to get specific miRNA expression levels of ceRNA-ceRNA interactions.
#'
#' The function will do segmentation with the method called "circular binary segmentation" to cluster noisy correlation into neighboring regions of distinct correlation levels. And then do the peak merging by considering low probability of multi-peaks occurring. Triplets with more than double peaks of segments if the distance of two peaks were smaller than a fixed value and there were significantly different between the average correlation of each peak would be merged until single or double peaks. Also, the segments with few samples detected by circular binary segmentation which are less than three and make no sense to the interactions would merge to a single peak. The sampling correlation coefficients were compared to that from the whole samples after being normalized by using Fisher transformation. mRNA pairs with peaks significantly different from the threshold and the segment would be reported as the candidate ceRNAs.       
#' @return list format with miRNAs,candidate ceRNA-ceRNA pirs,peak locations, and the number of samples occuring ceRNA-ceRNA interactions.
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


segcluster_peakmerge=function(cor_shreshold_peak,gene_sam,miRNA_total,d){
cor_shreshold_peak=0.85

start.time=proc.time()

count=1
Realdata2_tmp_ls=list()
for(mir in miRNA_total){
  print(which(mir==miRNA_total))
  Gene <- as.character(data.frame(dictionary[dictionary[,1]==mir,][[2]])[,1])
  Gene <- intersect(Gene,rownames(gene_sam))
  
  gene_pair <- combn(Gene,2)
  total_pairs <- choose(length(Gene),2)
  d <- readRDS(paste("/tempwork173/linwang/data_TCGA-LUSC_split/data_TCGA-LUSC_",which(mir==miRNA_total),".rds",sep=""))
  
  tmp <- NULL
  tmp <- tryCatch({
    foreach(p=1:total_pairs, .combine = "rbind")  %dopar%  {
      
      print(p)
      cand.ceRNA=c()
      location=list()
      r=gene_pair[1,p]
      s=gene_pair[2,p]
      triplet <- d[,c(1,p+1)]
      colnames(triplet) <- c("miRNA","corr")
      
      if(sum(is.na(triplet$corr)) ==0){

    
        CNA.object <- CNA(triplet$corr,rep(1,dim(triplet)[1]),triplet$miRNA)
        colnames(CNA.object) <- c("chrom","maploc",paste("gene",r,"and",s)) 
        result <- segment(CNA.object)
        

        if(sum(result$output$num.mark<=3)>=1){ ### 2
          tooshort <- which(result$output$num.mark<=3)
          num.mark <- c(0,cumsum(result$output$num.mark),last(cumsum(result$output$num.mark)))

          
          if(1 %in% diff(tooshort)){
            cc=1
            lag=0
            for(q in 1:(length(tooshort)-1)){
              if(tooshort[q+1]-tooshort[q]==1){
                result$output[tooshort[q],"loc.end"] <- result$output[tooshort[q+1],"loc.end"]
                result$output[tooshort[q],"seg.mean"] <- t(matrix(result$output[tooshort[q]:tooshort[q+1],"num.mark"]))%*%matrix(result$output[tooshort[q]:tooshort[q+1],"seg.mean"])/sum(result$output[tooshort[q]:tooshort[q+1],"num.mark"])
                result$output[tooshort[q],"num.mark"] <- sum(result$output[tooshort[q]:tooshort[q+1],"num.mark"])
                result$output[tooshort[q+1],] <- result$output[tooshort[q],]
                
                lag[cc] <- tooshort[q]
                cc <- cc+1
              }
              
            }
            result$output <- result$output[-lag,]
            row.names(result$output) <- 1:dim(result$output)[1]
            tooshort <- which(result$output$num.mark<=3)
          }
          if(length(tooshort)>=1){
            cc=1
            lag=c()
            for(t in 1:length(tooshort)){
              long_seg <- which(result$output$num.mark>3)
              diff=abs(tooshort[t]-long_seg)
              closest_seg <- long_seg[which(diff==min(diff))]
              if(length(closest_seg)>=2){
                b <- abs(result$output$seg.mean[closest_seg]-result$output$seg.mean[tooshort[t]])
                closest_seg <- closest_seg[b==min(b)]
              }
              result$output[tooshort[t],"loc.start"] <- min(result$output[tooshort[t],"loc.start"],result$output[closest_seg,"loc.start"])
              result$output[tooshort[t],"loc.end"] <- max(result$output[tooshort[t],"loc.end"],result$output[closest_seg,"loc.end"])
              result$output[tooshort[t],"seg.mean"] <- t(matrix(result$output[tooshort[t]:closest_seg,"num.mark"]))%*%matrix(result$output[tooshort[t]:closest_seg,"seg.mean"])/sum(result$output[tooshort[t]:closest_seg,"num.mark"])
              result$output[tooshort[t],"num.mark"] <- sum(result$output[tooshort[t]:closest_seg,"num.mark"])
              result$output[closest_seg,] <- result$output[tooshort[t],]
              
              lag[cc] <- tooshort[t]
              cc <- cc+1
   
            }
            result$output <- result$output[-lag,]
            row.names(result$output) <- 1:dim(result$output)[1]
            
            
          }
          
        }
        
        
        cand.corr <- c(-1,result$output$seg.mean,-1)  #add 0 to detect peaks happened at head and tail
        peak.loc <- findPeaks(cand.corr)-2

        no_merg_loc <- c()
        no_merg_count <- 1
        if(sum(cand.corr[peak.loc+1] > cor_shreshold_peak) >=2){ ### para 0.5
          
          for(i in 1:(length(peak.loc)-1)){
            if(sum(result$output[(peak.loc[i]+1):(peak.loc[i+1]-1),"num.mark"]) > w){
              no_merg_loc[no_merg_count] <- peak.loc[i]
            }
            
          }
          tryCatch({
            peak.loc <- peak.loc[-which(peak.loc==no_merg_loc)]
          },error=function(e){})
          
          
          while(sum(cand.corr[peak.loc+1] > cor_shreshold_peak) >=2){  ### para 0.5
            num.mark <- c(0,cumsum(result$output$num.mark),last(cumsum(result$output$num.mark)))
            TestPeak.pval <- c()
            for(i in 1:(length(peak.loc)-1)){
              z1 <- fisherz(mean(triplet$corr[(num.mark[peak.loc[i]]+1):num.mark[peak.loc[i]+1]],na.rm=T))
              z2 <- fisherz(mean(triplet$corr[(num.mark[peak.loc[i]]+1):num.mark[peak.loc[i+1]+1]],na.rm=T))
              N1 <- length(triplet$corr[(num.mark[peak.loc[i]]+1):num.mark[peak.loc[i]+1]])
              N2 <- length(triplet$corr[(num.mark[peak.loc[i]]+1):num.mark[peak.loc[i+1]+1]])
              TestPeak.pval[i] <- 2*pnorm(abs(z1-z2)/sqrt(1/(N1-3)+1/(N2-3)),lower.tail = FALSE)
            }
            
            if(sum(TestPeak.pval>0.05)!=0){ 
              TestPeak.p <- TestPeak.pval[TestPeak.pval>0.05]
              mergp.loc <- which(TestPeak.pval%in%TestPeak.p)
              distance <- c()
              for(i in 1:(length(peak.loc)-1)){
                distance[i] <- sum(result$output[(peak.loc[i]+1):(peak.loc[i+1]-1),"num.mark"])
              }
              peak_min <- mergp.loc[distance[mergp.loc]==min(distance[mergp.loc])]
              p_merg <- intersect(mergp.loc,peak_min)
              if(length(peak_min)>=2){
                peak_min <- mergp.loc[TestPeak.pval[p_merg]==min(TestPeak.pval[p_merg])]
              }
              peak_min <- peak_min[1]
              
              result$output[peak.loc[peak_min],"loc.end"] <- result$output[peak.loc[peak_min+1],"loc.end"]
              result$output[peak.loc[peak_min],"seg.mean"] <- t(matrix(result$output[peak.loc[peak_min]:peak.loc[peak_min+1],"num.mark"]))%*%matrix(result$output[peak.loc[peak_min]:peak.loc[peak_min+1],"seg.mean"])/sum(result$output[peak.loc[peak_min]:peak.loc[peak_min+1],"num.mark"])
              result$output[peak.loc[peak_min],"num.mark"] <- sum(result$output[peak.loc[peak_min]:peak.loc[peak_min+1],"num.mark"])
              result$output <- result$output[-c((peak.loc[peak_min]+1):peak.loc[peak_min+1]),]
              row.names(result$output) <- 1:dim(result$output)[1]
              
              cand.corr.new <- c(-1,result$output$seg.mean,-1)  #add 0 to detect peaks happened at head and tail
              peak.loc.new <- findPeaks(cand.corr.new)-2
              
              tryCatch({
                no_merg_loc <- c()
                no_merg_count <- 1
                for(i in 1:(length(peak.loc.new)-1)){
                  if(sum(result$output[(peak.loc.new[i]+1):(peak.loc.new[i+1]-1),"num.mark"])> w){
                    no_merg_loc[no_merg_count] <- peak.loc.new[i]
                  }  
                }
                peak.loc.new <- peak.loc.new[-no_merg_loc]
              },error=function(e){})
              
              if(length(peak.loc.new)==length(peak.loc)) break
              peak.loc <- peak.loc.new
              
              
            }else break
            
          }
          
        }
        num.mark <- c(0,cumsum(result$output$num.mark),last(cumsum(result$output$num.mark)))
        max_seg <- which(result$output$seg.mean==max(result$output$seg.mean))
        min_seg <- which(result$output$seg.mean==min(result$output$seg.mean))
        
        z1 <- fisherz(result$output$seg.mean[max_seg])
        z2 <- fisherz(result$output$seg.mean[min_seg])
        N1 <- result$output[max_seg,"num.mark"]
        N2 <- result$output[min_seg,"num.mark"]
        Test <- 2*pnorm(abs(z1-z2)/sqrt(1/(N1-3)+1/(N2-3)),lower.tail = FALSE)
  
        if(Test < 0.05){
          if(sum(cand.corr[peak.loc+1] > cor_shreshold_peak) >0 && sum(cand.corr[peak.loc+1] > cor_shreshold_peak) <=2){  ### para 0.5
            cand.ceRNA=paste(r,s) 
            
            tryCatch({
              peak.loc=sort(c(peak.loc,no_merg_loc)) #put back the peak that can't be merged 
            },error=function(e){})
            
            True_peak <- peak.loc[cand.corr[peak.loc+1] > cor_shreshold_peak]
            location=result$output[True_peak,c("loc.start","loc.end")]
            
            if(!is.null(cand.ceRNA)){
              lst <- list(miRNA=mir,cand.ceRNA=cand.ceRNA,location=location,numOfseg=result$output$num.mark[True_peak])
              lst
              
              
            }
          }
          
          
        }
        
      }
      
            
    }
    
  },error=function(e){
    e

  })

  Realdata2_tmp_ls[[count]] <- tmp
  count=count+1

}
return(Realdata2_tmp_ls)
}