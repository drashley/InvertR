#' finds ROI locations in collapsed wcRatio file
#'
#' locates putative inversions based on distribution of wcRatios and generates an ROI file
#' 
#' @param collapseFile Dataframe of collapsed wcRatios 
#' @param fileFrequencies Bedfile of reads used to calculate wcRatios
#' @param chrState Identify ROIs that deviate from chrState
#' @param baselineThreshold Minimum threshold used to call ROIs from wcRatio file 
#' @param elongationPenalty Extends roiRowList
#' @param minReads Minimimum reads required to call the ROI 
#' @param verbose Verbose messages if\code{TRUE} 
#' @author Ashley D. Sanders, Mark Hills
#' @export
findROILocation <- function(collapseFile, fileFrequencies, chrState=0, verbose=TRUE, baselineThreshold=0.8, elongationPenalty=3, minReads=4)
{
  
    #If collapseFile is in bedgraph format...
    if(collapseFile[1,1] == "track type=bedGraph")
    {
      collapseFile <- collapseFile[2:nrow(collapseFile),]
    }
  
    #set rownames from 1 to nrow of file
    rownames(collapseFile) <- c(1:nrow(collapseFile))
    #first calculate WC ratio of all 
    WCaverage <- mean(as.numeric(fileFrequencies[which(fileFrequencies[,6] > baselineThreshold),6]))
    Th <- round(WCaverage - 0.2, digits =2)
    callingThreshold <- Th*100
    
    rowCalls <- as.numeric(rownames(collapseFile[which(as.numeric(collapseFile[,4]) < callingThreshold),]))
    rowCalls <- c(rowCalls[1]-1, rowCalls, rowCalls[length(rowCalls)]+1)    
    
    roiRowList <- sort(c(rowCalls[1], rowCalls[which(diff(rowCalls) >= elongationPenalty)], rowCalls[which(diff(rowCalls) >= elongationPenalty)+1], rowCalls[length(rowCalls)]))
      
    ROIframe <- data.frame(chr=vector(), start=vector(), end=vector())
    ROIlocationTable <- data.frame(index=vector(), chr=vector(), callingTh=vector(), ROIstart=vector(), ROIend=vector(), ROIsize= vector(), deltaWC=vector(), roiReads=vector())
    
    if (length(roiRowList) > 0)
    {
      if(roiRowList[1] == 0){roiRowList[1] <- 1}
      #if(roiRowList[1] == 0){cycleStart<-2}else{cycleStart<-1}
        #for(row in seq(cycleStart,length(roiRowList)))
        for(row in seq(1,length(roiRowList)))
        {
          #Take modulo to pull out odd and even rows
          if(row %% 2 == 1)
          {
            ROIlocation <- data.frame(chr=collapseFile[roiRowList[row],1], start=as.numeric(collapseFile[roiRowList[row],2]), end=as.numeric(collapseFile[roiRowList[row+1], 3]), stringsAsFactors=FALSE)
            ROIframe <- rbind(ROIframe, ROIlocation)
          }
        }
  
      if(nrow(ROIframe) != 0)
      {
        if(verbose==T){message('List of ROIs generated')}
  
  ### 
        for(row in seq(1, nrow(ROIframe)))
        {
          if(row == 1)
          {
            #For the first row, analyze the region with the highest WC ratio from the start of the region (0), to the start of the ROI for upstream coordinate

            bestCallBefore <- unique(max(fileFrequencies[which(fileFrequencies[,3] <= ROIframe[row,2] & fileFrequencies[,3] >= 0),6]))
            preROIread <- tail(fileFrequencies[which(fileFrequencies[,3] <= ROIframe[row,2] & fileFrequencies[,3] >= 0 & fileFrequencies[,6] == bestCallBefore),],1)
            #message('first ROIstart found')  
            if(nrow(ROIframe) != 1)
            {
              #For first row, analyze the region with the highest WC ratio from the end of ROI from row 1, to the start of the ROI for row2
              bestCallAfter <- unique(max(fileFrequencies[which(fileFrequencies[,3] >= ROIframe[row,3] & fileFrequencies[,3] <= ROIframe[row+1, 2]),6]))			
              postROIread <- head(fileFrequencies[which(fileFrequencies[,3] >= ROIframe[row,3] & fileFrequencies[,3] <= ROIframe[row+1,2] & fileFrequencies[,6] == bestCallAfter),],1)
              #message('first ROI found')
            }else{
              #if only one row, analyze the region with the highest WC ratio from the end of ROI from row 1, to the end of chr
              bestCallAfter <- unique(max(fileFrequencies[which(fileFrequencies[,3] >= ROIframe[row,3] & fileFrequencies[,3] <= fileFrequencies[nrow(fileFrequencies), 3]),6]))
              postROIread <- head(fileFrequencies[which(fileFrequencies[,3] >= ROIframe[row,3] & fileFrequencies[,3] <= fileFrequencies[nrow(fileFrequencies), 3] & fileFrequencies[,6] == bestCallAfter),],1)
              #message('first ROI found')            
            }
          }
          else if(row == nrow(ROIframe))
          {
            #For the last row do the same as in 'else' for the before calls			
            bestCallBefore <- unique(max(fileFrequencies[which(fileFrequencies[,3] <= ROIframe[row,2] & fileFrequencies[,3] >= ROIframe[(row-1),3]),6]))
            preROIread <- tail(fileFrequencies[which(fileFrequencies[,3] <= ROIframe[row,2] & fileFrequencies[,3] >= ROIframe[row-1,3] & fileFrequencies[,6] == bestCallBefore),],1)
            #For last row, analyze the region with the highest WC ratio from the end of ROI from last row, to the end of the reads in file Frequencies
            bestCallAfter <- unique(max(fileFrequencies[which(fileFrequencies[,3] >= ROIframe[row,3] & fileFrequencies[,3] <= fileFrequencies[nrow(fileFrequencies), 3]),6]))
            postROIread <- head(fileFrequencies[which(fileFrequencies[,3] >= ROIframe[row,3] & fileFrequencies[,3] <= fileFrequencies[nrow(fileFrequencies), 3] & fileFrequencies[,6] == bestCallAfter),],1)
            #message('last ROI found')  
          }else{
            #For the all other rows, analyze the region with the highest WC ratio from the end of the previous ROI to the start of the  next ROI
            bestCallBefore <- unique(max(fileFrequencies[which(fileFrequencies[,3] <= ROIframe[row,2] & fileFrequencies[,3] >= ROIframe[(row-1),3]),6]))
            preROIread <- tail(fileFrequencies[which(fileFrequencies[,3] <= ROIframe[row,2] & fileFrequencies[,3] >= ROIframe[row-1,3] & fileFrequencies[,6] == bestCallBefore),],1)
            #For all other rows, analyze the region with the highest WC ratio from the end of ROI from row to the start of ROI from next row
            bestCallAfter <- unique(max(fileFrequencies[which(fileFrequencies[,3] >= ROIframe[row,3] & fileFrequencies[,3] <= ROIframe[row+1, 2]),6]))
            postROIread <- head(fileFrequencies[which(fileFrequencies[,3] >= ROIframe[row,3] & fileFrequencies[,3] <= ROIframe[row+1,2] & fileFrequencies[,6] == bestCallAfter),],1)
            #message('next ROI found')
          }
         
    ###      
          #This is the last read before inversion that has the highest ratio.
          #preROIState <- preROIread[1,4]
          if(chrState == 'cc'){preROIState<- '+'}else if(chrState == 'ww'){preROIState <- '-'}else{preROIState <- preROIread[1,4]}

          #pull out locations of all reads mapping to the opposite strand to the preROI state
          oppositeStrands <- fileFrequencies[which(fileFrequencies[,3] < postROIread[,3] & fileFrequencies[,3] > preROIread[,3] & fileFrequencies[,4] != preROIState & fileFrequencies[,3] <= ROIframe[row,3]), ]
          #First, if there's less than 4 reads, ignore this ROI; it will never have >20% WC
          if(nrow(oppositeStrands) > 4)
          {
            
            ######~~~~~~ ANALYSIS FOR RETRIEVING ROI START LOCATION ~~~~~~~######
            #####################################################################
            # First, work out 20% WCratio downstream from start
            upWCratio <- head(fileFrequencies[which(fileFrequencies[,3] >= oppositeStrands[1,3]),],20)
            upWCNum <- nrow(upWCratio[which(upWCratio[,4] != preROIState),])*5
            
            #while these criteria are NOT met, move forward to the next read that doesn't match preROI state...
            while(upWCNum <= 20 & nrow(upWCratio) == 20)
            {
              upWCratio <- head(fileFrequencies[which(fileFrequencies[,3] > upWCratio[1,3] & fileFrequencies[,4] != preROIState),],20)
              upWCNum <- nrow(upWCratio[which(upWCratio[,4] != preROIState),])*5				
            }
            #Now calculate upstream reads. Make sure all 10 reads upstream match the preROIState
            upPureRatio <- tail(fileFrequencies[which(fileFrequencies[,3] < upWCratio[1,3]),],10)
            upPureNum <- nrow(upPureRatio[which(upPureRatio[,4] == preROIState),])*10
            
            #while these criteria are NOT met, move backward to the next read that does match preROI state...
            while(upPureNum != 100 & nrow(upPureRatio) == 10 )
            {
              upPureRatio <- tail(fileFrequencies[which(fileFrequencies[,3] < upPureRatio[10,3] & fileFrequencies[,4] == preROIState),],10)
              upPureNum <- nrow(upPureRatio[which(upPureRatio[,4] == preROIState),])*10
            }
            
            #Add caveat that if <10 reads available, make upPureNum =100
            if(head(upPureRatio$pos,1) == head(collapseFile$start,1)){upPureNum <- 100}

            ######~~~~~~ ANALYSIS FOR RETRIEVING ROI END LOCATION ~~~~~~~######
            ###################################################################
            # And 20% WC upstream from end
            downWCratio <- tail(fileFrequencies[which(fileFrequencies[,3] <= oppositeStrands[nrow(oppositeStrands),3]),],20)
            downWCNum <- nrow(downWCratio[which(downWCratio[,4] != preROIState),])*5
            
            #while these criteria are NOT met, move forward to the next read that doesn't match preROI state...
            while(downWCNum <= 20 & nrow(downWCratio) == 20)
            {
              downWCratio <- tail(fileFrequencies[which(fileFrequencies[,3] < downWCratio[1,3] & fileFrequencies[,4] != preROIState),],20)
              downWCNum <- nrow(downWCratio[which(downWCratio[,4] != preROIState),])*5				
            }
            #Now calculate upstream reads. Make sure all 10 reads upstream match the preROIState
            downPureRatio <- head(fileFrequencies[which(fileFrequencies[,3] > downWCratio[nrow(downWCratio),3]),],10)
            downPureNum <- nrow(downPureRatio[which(downPureRatio[,4] == preROIState),])*10
            
            #while these criteria are NOT met, move backward to the next read that does match preROI state...
            while(downPureNum != 100 & nrow(downPureRatio) == 10)
            {
              downPureRatio <- head(fileFrequencies[which(fileFrequencies[,3] > downPureRatio[1,3]),],10)
              downPureNum <- nrow(downPureRatio[which(downPureRatio[,4] == preROIState),])*10
            }
            
      ###########
            if(upWCNum >= 20 & upPureNum == 100 & downWCNum >= 20 & downPureNum == 100)
            {
              ROIstartLoc <- tail(upPureRatio[,3],1)
              ROIendLoc <- head(downPureRatio[,3],1)
            
              if(ROIstartLoc < ROIendLoc)
              {
                #Calculate a new WC ratio Average for the refined ROI
                #ROIWCaverage <- mean(fileFrequencies[which(fileFrequencies[,3] > ROIstartLoc & fileFrequencies[,3] < ROIendLoc), 6])
                ## CORRECTION WORKING!
                ROIoppositereads <- nrow(fileFrequencies[which(fileFrequencies[,3] > ROIstartLoc & fileFrequencies[,3] < ROIendLoc & fileFrequencies[,4] != preROIState),])
                ROIallreads <- nrow(fileFrequencies[which(fileFrequencies[,3] > ROIstartLoc & fileFrequencies[,3] < ROIendLoc),])
                ROIWCaverage <- 1-(ROIoppositereads/ ROIallreads)
               
                #To account for background, subtract the 'background' WCaverage from the ROIWCaverage.
                deltaWC <- round(WCaverage - ROIWCaverage, digits=3)
                #Calculate read depth of this region
                roiReads <- round(nrow(fileFrequencies[which(fileFrequencies[,3] >= ROIstartLoc & fileFrequencies[,3] <= ROIendLoc ),]), digits=2)
                
                #now make it into a spectacular table.  Hurrah!
                ROIlocationLine <- data.frame(index=upPureRatio[1,1], chr=upPureRatio[1,2], callingTh=Th, ROIstart=ROIstartLoc, ROIend=ROIendLoc, ROIsize=(ROIendLoc-ROIstartLoc), deltaWC=deltaWC, roiReads=roiReads)
                #ROIlocationTable <- rbind(ROIlocationTable, ROIlocationLine)
                
                if (ROIlocationLine$deltaWC > 0.3 && ROIlocationLine$roiReads > minReads){ ROIlocationTable <- rbind(ROIlocationTable, ROIlocationLine) }
                if(verbose==T){message(paste('ROI mapped: (No ', row, '/', nrow(ROIframe), ' possible ROIs)', sep=""))}
              
              }else{ if(verbose==T){message("no inversion here!  move along...")} }

            }else{ if(verbose==T){message("inverion breakpoints not clear!  move along...")} }  
    ####         
          }else{ if(verbose==T){message("region too small! moving along...")} }
        }
  
  
      }else{ if(verbose==T){message("no putative inversions seen in this region")} }
  
  
  
      if (nrow(ROIlocationTable) == 0)
      {
        if(verbose==T){message("CHROMOSOME DEVOID OF INVERSIONS") }
        ROIlocationTable <- 0
      }    
      return(list(ROIlocationTable, Th))
    }
  if(verbose==T){message("NO ROIs identified, this region is DEVOID OF INVERSIONS") }
  ROIlocationTable <- 0

  return(list(ROIlocationTable, Th))
}
