#' Processes input bamfile
#'
#' Generates a dataframe of reads from a bam file for the specified genomic coordinate
#' 
#' @param startLoc Genomic coordinate start position
#' @param endLoc Genomic coordinate end position
#' @param chr Genomic coordinate chromosome
#' @param fileName Name of bamfile to process
#' @param padding Number of bases to extend beyond the startLoc and endLoc
#' @param qual Filter reads based on specified qulaity score 
#' @param rmdup if\code{TRUE} remove duplicate reads
#' @param verbose Verbose messages if\code{TRUE} 
#' @author Ashley D. Sanders, Mark Hills
#' @export
processBam <- function(startLoc, endLoc, chr, fileName, padding=100000, qual=0, rmdup=TRUE, verbose=TRUE)
{
  #library(Rsamtools)
  #This function only takes reads from region of interest
	if(verbose){message(paste('-> Creating processFile for index ', fileName, sep=""))}
	
	#read in first reads & create dataframe, then convert to table
	temp.1.bam <- scanBam(fileName, param=ScanBamParam(flag=scanBamFlag(isFirstMateRead=TRUE, isDuplicate=F, isProperPair=T), what=c("rname","pos","strand","mapq")))
	temp.1.dataframe <- do.call("DataFrame", temp.1.bam)
	temp.1.dataframe <- temp.1.dataframe[which(temp.1.dataframe$rname == as.character(chr) & temp.1.dataframe$pos >= startLoc-padding & temp.1.dataframe$pos <= endLoc+padding),]
	if(qual != 0)
	{
		temp.1.dataframe <- temp.1.dataframe[which(temp.1.dataframe$mapq >= qual),]
	}
	if(rmdup)
	{
		temp.1.dataframe <- temp.1.dataframe[!duplicated(paste(temp.1.dataframe$rname, temp.1.dataframe$pos)),]
	}

	#read in second pair reads, create dataframe then convert to table
	temp.2.bam <- scanBam(fileName, param=ScanBamParam(flag=scanBamFlag(isSecondMateRead=TRUE,isDuplicate=F, isProperPair=T), what=c("rname","pos","strand","mapq")))
	temp.2.dataframe <- do.call("DataFrame", temp.2.bam)
	temp.2.dataframe <- temp.2.dataframe[which(temp.2.dataframe$rname == as.character(chr) & temp.2.dataframe$pos >= startLoc-padding & temp.2.dataframe$pos <= endLoc+padding),]
	if(qual != 0)
	{
		temp.2.dataframe <- temp.2.dataframe[which(temp.2.dataframe$mapq >= qual),]
	}
	if(rmdup)
	{
		temp.2.dataframe <- temp.2.dataframe[!duplicated(paste(temp.2.dataframe$rname, temp.2.dataframe$pos)),]
	}
	
  #Now, flip all reverse complement reads to the different strand
	#add a new level to the factor
	levels(temp.1.dataframe[,2]) <- c(levels(temp.1.dataframe[,2]), "1")
	temp.1.dataframe$strand <- replace(temp.1.dataframe$strand, temp.1.dataframe$strand == '+', '1')
	temp.1.dataframe$strand <- replace(temp.1.dataframe$strand, temp.1.dataframe$strand == '-', '+')
	temp.1.dataframe$strand <- replace(temp.1.dataframe$strand, temp.1.dataframe$strand == '1', '-')  
  
	bamDataframe <- as.data.frame(rbind(temp.1.dataframe, temp.2.dataframe))
	bamDataframe <- bamDataframe[order(bamDataframe$pos),]
	bamDataframe <- data.frame(rname=bamDataframe[,1], pos=bamDataframe[,3], strand=bamDataframe[,2], mapq=bamDataframe[,4]) 

	chrState<- round((table(bamDataframe$strand)[2]-table(bamDataframe$strand)[1])/nrow(bamDataframe), digits=3)

  return(list(bamDataframe, chrState))

}