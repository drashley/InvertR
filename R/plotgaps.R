#' Plots gaps in output png  
#'
#' Plots gapfile above histogram in plotROIFrequencies
#'
#' @param gapfile File with gap locations
#' @param chr Genomic coordinate chromosome
#' @param binSize The number of reads in each bin used to calculate wcRatio
#' @param startLoc Genomic coordinate start position
#' @param endLoc Genomic coordinate end position
#' @param padding Number of bases to extend beyond genomic coordinates listed in regionTable
#' @author Ashley D. Sanders, Mark Hills
#' @export
plotGaps <- function(gapfile, chr, startLoc, endLoc, padding)
{
  gapByChr <- gapfile[which(gapfile[,2] == chr & gapfile[,3] >= startLoc-padding & gapfile[,4] <= endLoc+padding),]
  # gapByChr <- gapfile[which(gapfile[,2] == chr),]
  if(nrow(gapByChr) != 0)  
  {   
    #message('  (plotting gaps in the region)')
    for(gap in seq(1, nrow(gapByChr)))
    {
      if(as.character(gapByChr[gap,8]) != 'centromere' | as.character(gapByChr[gap,8]) != 'heterochromatin')
      {
        #rect(gapByChr[gap,3], 1.01, gapByChr[gap,4], 1.1, border = "gold2", col= "gold2", lwd=0.1)
        rect(gapByChr[gap,3], 1.03, gapByChr[gap,4], 1.1, border = "grey48", col= "grey48", lwd=0.1)
      }else{
        #rect(gapByChr[gap,3], 1.02, gapByChr[gap,4], 1.1, border = "gold4", col= "gold4", lwd=0.1)
        rect(gapByChr[gap,3], 1.02, gapByChr[gap,4], 1.1, border = "grey20", col= "grey20", lwd=0.1)
      }  
    }  
  }else{
    message('  (no gaps here!)')
  }
}