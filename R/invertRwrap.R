
#' Wrapper function for InvertR
#'
#' This script will move through .bam or .bed files in a folder and perform several steps (see Details).
#'
#' 1. calculate the WCratio chromosome-by-chromosome
#' 2. Locate the ROIs in chromosomes passing the WCcutoff
#' 3. write a bedgraph file of wcRatios -> can upload on to UCSC Genome browser
#' 4. write an ROI file for each index with all chromosomes included
#' 
#' @param regionTable Genomic coordinates to be analyzed (ROI list or Chr Table)
#' @param dataDirectory Output directory. If non-existent it will be created
#' @param binSize The number of reads in each bin used to calculate wcRatio
#' @param WCcutoff The number of watson or crick reads used to define chrStates
#' @param gapfile Input txt file of gaps in the genome
#' @param type File input type, either 'bed' or 'bam'
#' @param dup If \code{TRUE}, removes duplicate reads
#' @param qual Filter reads based on specified qulaity score
#' @param padding Number of bases to extend beyond genomic coordinates listed in regionTable
#' @param verbose If \code{TRUE} Verbose messages
#' @param strand If \code{TRUE} Plot the crick and watson reads above the histrogram
#' @param png If\code{TRUE} Generates a png figure of each file
#' @param findROIs If\code{TRUE} Runs findROIlocations to locate putative inversions
#' @param ROI If\code{TRUE} Expects ROI list, if\code{FALSE} Expects chromosome table
#' @param minDepth The minimum number of reads/Mb required to analyze the chromosome
#' @param minReads The minimum number of reads within the ROI required for inclusion
#' @author Ashley D. Sanders, Mark Hills
#' @export
runInvertR <- function(regionTable, binSize=50, WCcutoff=0.75, dataDirectory='./InvertR_analysis/', gapfile=0, type='bed', dup=TRUE, qual=10, padding=0, minDepth=20, minReads=20, verbose=TRUE, png=TRUE, strand=TRUE, ROI=FALSE, genotype=TRUE, findROIs=T)
{#
  options(warn=-1)
  if(type == 'bam') {library(Rsamtools)}
  
  dir.create(dataDirectory)
  fileDestination <- dataDirectory
  
  #for every chromosome...
  for(i in seq(1,nrow(regionTable)))
  {##
    ch <- regionTable[i,1]
    chr <- paste('chr', regionTable[i,1], sep="")
    startLoc <- regionTable[i,2]
    endLoc <- regionTable[i,3]
    if (chr == 'chrY') { WCcutoff = 0}
    ## NOTE change WCcutoff=0 if chrY (b.c. cannot have a WC chr, and may have large inversions (e.g. cad11) which would be missed if wcCutoff high)
    
    if (ROI == TRUE)
    {
      ROIname <- paste('ROINo.', i, sep="")
      dir.create(paste(fileDestination, ROIname, sep=""))
      chrfileDestination <- (paste(fileDestination, ROIname, '/', sep=""))
      dir.create(paste(fileDestination, 'WCLibs', sep=""))
      padding<- round((endLoc-startLoc)*0.33, digits=0)
    }else{
      ROIname <- 'wholeChr'
      dir.create(paste(fileDestination, chr, sep=""))
      chrfileDestination <- (paste(fileDestination, chr, '/', sep=""))
      dir.create(paste(fileDestination, 'WCLibs', sep=""))
      padding<- 0
    }
    
    
    pattern <- paste('.', type,'$', sep="")
    fileList <- list.files(path='.', pattern=pattern, full.names=TRUE)
    fileLength <- length(fileList)
    indexCounter <- 0
    options(scipen=20)
    #for reading in multiple files at a particular location
    allFrequencies <- data.frame(index=vector(), rname=vector(), pos=vector(), strand=vector(), mapq=vector(), WCratio=vector(), chrState=vector())
    #allROIlocationTable <- data.frame(index=vector(), chr=vector(), ROIstart=vector(), ROIend=vector(), deltaWC=vector(), roiReadDepth=vector())	
    allROIlocationTable <- data.frame(index=vector(), chr=vector(), callingTh=vector(), ROIstart=vector(), ROIend=vector(), ROIsize= vector(), deltaWC=vector(), roiReads=vector())
    wcLibraries <- data.frame(index=vector(), chr=vector(), wcCall=vector())
    
    #for every filename....
    for(fileName in fileList)
    { ###
      indexCounter <- indexCounter + 1
      if(verbose==T){message(paste('** RUNNING ', fileName, ' [lib.No ', indexCounter,'/', fileLength, '], ', ' chromosome [', ch, '/', nrow(regionTable), '] **', sep=""))}	
      index <- basename(fileName)
      
      #read in files; either using processBed or processBam
      if(type == 'bed')
      {
        tempFile <- processBed(startLoc, endLoc, chr, fileName, qual=qual, rmdup=dup, padding=padding, verbose=verbose)   
        chrState <- tempFile[[2]]
        if(chrState >= WCcutoff) {chrState<-'ww'}else if(chrState <= -WCcutoff){chrState<-'cc'}else{chrState<-'wc'}
        # if chrState is Negative chr is CRICK/CRICK
        processFile <- tempFile[[1]]
      }else if(type == 'bam') {
        tempFile <- processBam(startLoc, endLoc, chr, fileName, qual=qual, rmdup=dup, padding=padding, verbose=verbose)
        chrState <- tempFile[[2]]
        if(chrState >= WCcutoff) {chrState<-'ww'}else if(chrState <= -WCcutoff){chrState<-'cc'}else{chrState<-'wc'}
        # if chrState is Negative chr is CRICK/CRICK
        processFile <- tempFile[[1]]
        processFile<- cbind(chr, processFile[2:length(processFile)]) # pastes chr instead of ch to file
        processFile<- processFile[!duplicated(processFile[2]),]
      }
      if(verbose==T){message(paste('-> bedFile generated for ', index, ', chromosome ', ch, sep=""))}
    
      if(length(processFile[[1]]) > 1)
      {####        
        #Filters out low (< minDepth) read depth libraries. If enough reads are in this library, proceed...
        if(length(processFile[[1]]) > 1 && nrow(processFile)/((endLoc-startLoc)/1000000) > minDepth)
        {
          #calculate the ration of _ to + reads (i.e. the wcCall)
          if (ROIname == 'wholeChr'){
            wcCall <- round((table(processFile$strand)[2]-table(processFile$strand)[1])/nrow(processFile), digits=3)
            #wcCall <- chrState
          }else{
            tempFile <- processFile[which(processFile$pos < startLoc),]
            tempFile <- rbind(tempFile, processFile[which(processFile$pos > endLoc),])
            #filters reads that flank the ROI to calculate the wcCall of these surrounding reads (since an inversion at the ROI will impact the wcCall)
            wcCall <- round(( table(tempFile$strand)[2]-table(tempFile$strand)[1] ) /nrow(tempFile), digits=3)
          }			
          
          if( is.na(wcCall)) {wcCall <- 1}
          
          ######
          # wcCall can become NA if 100% of reads are + or -
          if(wcCall != 'NaN')
          {
            if(wcCall <= -WCcutoff | wcCall >= WCcutoff)
            { ##This is a pure (WW or CC) library...
              
              # calculates A WCratio value FOR EACH READ based on the proportion of W and C in (binSize #) succeeding reads 
              fileFrequencies <- countFreqs(processFile, checkNum=binSize, verbose=FALSE)
              if(verbose==T){message(paste('-> fileFrequencies counted for ', index, ' file, chromosome ', chr, sep=""))}
              if(length(fileFrequencies) > 1 && nrow(fileFrequencies) > 2)
              {			     
                # reduces the table size by identifying only the locations where the WCratio values change
                outputFile <- collapseStrands(fileFrequencies, index, asBedgraph=TRUE)				
                if(verbose==T){message(paste('-> strands collapsed for ', index, ' file, chromosome ', chr, sep=""))}
                
                #find the location of ROIs that dip below threshold level
                fileFrequencies <- cbind(fileName, fileFrequencies) 
                fileFrequencies<- cbind(fileFrequencies, chrState)
                allFrequencies <- rbind(allFrequencies, fileFrequencies)
                
                if(findROIs==T){ # if true then run findROILocation script, else ROIlocationTable =1
                  #minReads specifies the minimum number of reads within the roi that are required to include it in the ROIlocationTable list
                  locationFile <- findROILocation(outputFile, fileFrequencies, chrState=chrState, verbose=verbose, baselineThreshold=0.8, minReads=minReads)
                  ROIlocationTable <- locationFile[[1]]
                  # ROIlocationTable <- data.frame(index=vector(), chr=vector(), callingThreshold=vector(), ROIstart=vector(), ROIend=vector(), ROIsize= vector(), deltaWC=vector(), roiReads=vector())
                  Th <- locationFile[[2]]
                }else{ROIlocationTable = 1
                      Th<- 1}
                
                if(length(ROIlocationTable) != 1)
                {
                  if(verbose==T){message(paste('-> Total of: ', nrow(ROIlocationTable), ' ROIs found for ', index, ' file, chromosome ', chr, sep=""))}
                  allROIlocationTable <- rbind(allROIlocationTable, ROIlocationTable)
                  deltaWC<-  ROIlocationTable[1,7] 
                  
                } else { deltaWC <- 0 }
                #if(verbose==T){message(paste('-> NO ROIs found for ', index, ' file, chromosome ', chr, sep=""))} 
                
                
                ### generates a png  of the single Ss library with gaps and ROI locations highlighted, also calculates the number of reads in the region -> table
                plotROIFrequencies(chrfileDestination, index, fileFrequencies, chr, binSize, startLoc, endLoc, ROIname=ROIname, gapfile=gapfile, callingThreshold=Th, ROI=ROIlocationTable, padding=padding, strand=strand, png=png)
              
                ###### write a bedfile of the reads, and then append the bedgraph
                bedfile<- cbind(chr,  processFile$pos, processFile$pos+100, index, processFile$mapq, as.character(processFile$strand))
                head<- paste('track name=', index, '_reads_', chr, ' visibility=1 colorByStrand="103,139,139 243,165,97"', sep="")
                write.table(head, file=paste(chrfileDestination, index, '_', chr, '_(b=', binSize, ', t=', Th, ').bedgraph', sep=""), row.names=FALSE, col.names=F, quote=F, append=F)   
                write.table(bedfile, file=paste(chrfileDestination, index, '_', chr, '_(b=', binSize, ', t=', Th, ').bedgraph', sep=""), row.names=FALSE, col.names=F, quote=F, append=T)
                
                ####  bedgraph of the collapsed WCratios for the single Ss library -> can be uploaded onto the ucsc genome browser
                write.table(outputFile, file=paste(chrfileDestination, index, '_', chr, '_(b=', binSize, ', t=', Th, ').bedgraph', sep=""), row.names=FALSE, quote=FALSE, col.names=FALSE, append=T)
              }
            } else {
              message(paste('~>  ', index, ' is WC for ', chr, ' - moving on to next lib', sep=""))
              wcLibrary <- cbind(index, chr, round(wcCall, digits=3))
              wcLibraries <- rbind(wcLibraries, wcLibrary)
            }
          } else { if(verbose==T){message(paste('~>  ', index, ' wcCall = NaN for ', chr, ' - moving on to next lib', sep="")) } }
          ######
        } else { if(verbose==T){message(paste('~~>  ', index,  'read count below minReads, moving on to next lib', sep="")) } }
        ####	
      } else { if(verbose==T){message(paste('~~>  ', index, ' has no reads in the region - moving on to next lib', sep="")) } }
    }###
    
    if(findROIs==T){
    write.table(allROIlocationTable, file=paste(fileDestination, 'ROI_locations_Table_b', binSize, '_', chr, '.txt', sep=""), row.names=FALSE, quote=FALSE, append=FALSE)  }
    
    write.table(allFrequencies, file=paste(fileDestination, 'allFrequencies_b', binSize, '_', chr, '.txt', sep=""), row.names=FALSE, quote=FALSE, append=FALSE) 
    write.table(wcLibraries, file=paste(fileDestination, '/WCLibs/wcLibraries_', chr, '.txt', sep=""), row.names=FALSE, quote=FALSE, append=FALSE)  
    
    # calculate Stats for the plots:
    if (nrow(allROIlocationTable) != 0){
      AvTh <- round(mean(allROIlocationTable[,3]), digits=2)
    }else{
      AvTh <-0
      allROIlocationTable<-0}
    
    if(nrow(allFrequencies) < 1){ allFrequencies <- data.frame("fileName", chr, startLoc, endLoc, '+', 0, 1, 'cc' )}
    plotROIFrequencies(fileDestination, 'overlay', allFrequencies, chr, binSize, startLoc, endLoc, ROIname=ROIname, gapfile=gapfile, callingThreshold=AvTh, ROI=allROIlocationTable, padding=padding, strand=FALSE, png=png)
    
    if(verbose==T){message(paste(' ~~ Overlaid plot generated for ', chr, ' *YIPEE!* moving on to next chromosome...', sep=""))}
  }##
  
}#




