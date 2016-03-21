#' Generates plot of wcRatios 
#'
#' Generates single plot of wcRatios for fileFrequencies
#'
#' @param fileLocation Specifies directory to save file
#' @param fileName Unique identifier of input file
#' @param fileFrequencies Dataframe of wcRatio used to generate histrogram
#' @param chr Genomic coordinate chromosome
#' @param binSize The number of reads in each bin used to calculate wcRatio
#' @param startLoc Genomic coordinate start position
#' @param endLoc Genomic coordinate end position
#' @param ROIname Genomic coordinate number 
#' @param padding Number of bases to extend beyond genomic coordinates listed in regionTable  e.g. [10000]
#' @param gapfile Input file of gaps to plot
#' @param callingThreshold Threshold line used to find ROIs
#' @param ROI ROI location table 
#' @param strand if\code{TRUE} plot reads above histrogram
#' @param png if\code{TRUE} Save png
#' @author Ashley D. Sanders, Mark Hills
#' @export
plotROIFrequencies <- function(fileLocation, fileName, fileFrequencies, chr, binSize, startLoc, endLoc, ROIname='wholeChr', gapfile=0, callingThreshold=0.66, ROI=0, padding=0, strand=TRUE, png=TRUE) 
{
  #create an empty png
    if (png==TRUE)
  	{ png(paste(fileLocation, ROIname,'_', fileName, '_', chr,  '_bin', binSize, '_t', callingThreshold, '.png', sep=""), width=6, height=2.5, units='in', res=800) }
  
  if(length(gapfile) != 1 | strand==TRUE)
  {
    ylim=c(-0.1,1.1)
  }else{
    ylim=c(-0.1,1)
  }
  xlim=c(signif(startLoc-padding, digits=4), signif(endLoc+padding, digits=4))

  plot(NA,  xlab="chromosome location", ylab="WC Ratio", ylim=ylim, xlim=xlim, yaxt="n", xaxt="n", cex.lab=0.8)
  abline(h=c(0,1), col="black", lty=2, lwd=0.3)
  abline(v=c(startLoc, endLoc), col="gray36", lty=3, lwd=0.3)
  axis(2, at=c(0, 0.5, 1), labels=NULL, lty=1, las=2, cex.axis=0.5)
  axis(1, cex.axis=0.5)
  title(main= paste(ROIname,'_', fileName, '_', chr,  '_bin', binSize, sep=""),cex.main=0.8, col.main="gray36")

  # colour Scheme
  ## het inversion:
  col_a <- "cyan4"
  ## hom inversion:
  col_b <- "darkorchid4"
  ## pb text and lines
  col_c <- "darkred"
  
  
  if(callingThreshold != 'na')
  {
    abline(h=callingThreshold, col='gray60', lty=4, lwd=0.2)
  }
  startLoc<- as.numeric(startLoc)
  endLoc<- as.numeric(endLoc)
  if(ROIname != 'wholeChr'){
    size<- round(((endLoc+padding) - (startLoc-padding))/1000, digits=2)
    sizeText <- paste(size, ' kb', sep="")
  }else{
    size<- round(((endLoc+padding) - (startLoc-padding))/1000000, digits=2)
    sizeText <- paste(size, ' Mb', sep="")
  }  
  
  # plots ROIs below the graph
  if(length(ROI) != 1)
  {    
    HetroiNo <- 0
    HomroiNo <- 0
    roiNo <-0
    for( row in seq(1, nrow(ROI)))
    {
      ## this will plot and number each ROI, colour is red (col_c) here
      mtext(row, cex=0.2, side=1, col=col_c, at=(as.numeric(ROI[row,4])+((as.numeric(ROI[row,5])-as.numeric(ROI[row,4]))/2)), line=0.4)
      #het inv (col_a)
      #ROIlocationTable <- data.frame(index=vector(), chr=vector(), callingThreshold=vector(), ROIstart=vector(), ROIend=vector(), ROIsize= vector(), deltaWC=vector(), roiReads=vector())                
      if(ROI[row,7] >= 0.33 & ROI[row,7] <= 0.66)
      {
        rect(ROI[row,4], -0.12, ROI[row,5], -0.03, border = col_c, lwd=0.1, col= col_c)				
       # abline(v=c(ROI[row,4], ROI[row,5]), col=col_c, lty=3, lwd=0.25) # vertical line at breakpoint
        # this will place bp location at each ROI
        if(ROIname != 'wholeChr'){
        #  mtext(ROI[row,4], cex=0.25, side=1, col=col_c, at=ROI[row,4])
        # mtext(ROI[row,5], cex=0.25, side=1, col=col_c, line= 0.25, at=ROI[row,5])
        #  mtext(paste('ΔWC =', ROI[row,7], sep=""), cex=0.35, side=3, col=col_a, at=ROI[row,4]+((ROI[row,5]-ROI[row,4])/2), line=0.4)
        }
        
        HetroiNo <- HetroiNo+1
        roiNo <- roiNo+1
      #hom inv (col_b)
      } else if( ROI[row,7] > 0.66) {
        rect(ROI[row,4], -0.12, ROI[row,5], -0.03, border = col_c, lwd=0.1, col= col_c)					
      #  abline(v=c(ROI[row,4], ROI[row,5]), col=col_c, lty=3, lwd=0.25) # vertical line at breakpoint
        # this will place bp location of each ROI
        if(ROIname != 'wholeChr'){
         # mtext(ROI[row,4], cex=0.25, side=1, col=col_c, line= 0, at=ROI[row,4])
         # mtext(ROI[row,5], cex=0.25, side=1, col=col_c, line= 0.25, at=ROI[row,5])
         # mtext(paste('ΔWC =', ROI[row,7], sep=""), cex=0.35, side=3, col=col_b, at=ROI[row,4]+((ROI[row,5]-ROI[row,4])/2), line=0.4)
         } 
        HomroiNo <- HomroiNo+1
        roiNo <- roiNo+1
      }
    }
    title(sub = paste(sizeText, ' region | Av reads/Mb = ', round(nrow(fileFrequencies)/(((endLoc+padding)-(startLoc-padding))/1000000), digits=1), ' | ', roiNo, ' ROIs (', HetroiNo, ' het & ', HomroiNo, ' hom)', sep=""), cex.sub=0.4, col.sub="maroon", adj=1)

  } else {
    roiNo <- 0
    title(sub = paste(sizeText, ' region | Av reads/Mb = ', round(nrow(fileFrequencies)/(((endLoc+padding)-(startLoc-padding))/1000000), digits=1), ' | ', roiNo, ' ROIs', sep=""), cex.sub=0.4, col.sub="maroon", adj=1)
    
  }

  # removes bp lines from gap file
  rect(startLoc,1.05,endLoc,1.1, col="white", border=0.01)
  mtext('Het Inv', cex=0.35, side=3, col=col_a, adj=1, line=0.4)
  mtext('Hom Inv', cex=0.35, side=3, col=col_b, adj=1, line=0.1)
  

  # plots gaps above the graph
  if(length(gapfile) != 1)
  {
    plotGaps(gapfile, chr, startLoc, endLoc, padding)
    abline(h=c(0,1), col="black", lty=1, lwd=0.3)
  }
    
  # plots Strand-seq strands above the graph
  if (strand==TRUE)
  {
    for (line in seq(1:nrow(fileFrequencies)))
    {
      if (fileFrequencies[line,4] == '-')
      {
        rect(fileFrequencies[line,3], 1.015, fileFrequencies[line,3]+100, 1.1, border = "sandybrown", col= "sandybrown", lwd=0.5)
      }else{
        rect(fileFrequencies[line,3], 1.015, fileFrequencies[line,3]+100, 1.1, border = "paleturquoise4", col= "paleturquoise4", lwd=0.5)
      }
    }
  }

  # plot the WC ratios for each library
  libNo <- 0
  for(lib in seq(1:length(unique(fileFrequencies[,1])) ) )
  {
    l <- unique(fileFrequencies[,1])[lib]
    lineIndex <- fileFrequencies[which(fileFrequencies[,1] == l),]
    lines(lineIndex[,3], round(lineIndex[,6], digits=1), lwd=0.4, col="black")
    libNo <- libNo + 1
  
  }

   if (png==TRUE){ graphics.off() }
}

