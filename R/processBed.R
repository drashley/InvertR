#' Processes input befile
#'
#' Generates a dataframe of reads from a bedfile for the specified genomic coordinate
#' 
#' @param startLoc Genomic coordinate start position
#' @param endLoc Genomic coordinate end position
#' @param chr Genomic coordinate chromosome
#' @param fileName Name of bamfile to process
#' @param skip Number of rows in header of bedfile
#' @param padding Number of bases to extend beyond the startLoc and endLoc
#' @param qual Filter reads based on specified qulaity score 
#' @param rmdup if\code{TRUE} remove duplicate reads
#' @param verbose Verbose messages if\code{TRUE} 
#' @author Ashley D. Sanders, Mark Hills
#' @export
processBed <- function(startLoc, endLoc, chr, fileName, skip=3, padding=100000, qual=0, rmdup=TRUE, verbose=TRUE)
{
	if(verbose){message(paste('-> Creating processFile for index ', fileName, sep=""))}
	bedFile <- read.table(fileName, skip=skip)
	temp.dataframe <- data.frame(rname=bedFile[,1], pos=bedFile[,2], strand=bedFile[,6], mapq=bedFile[,5]) 
	
	if(qual != 0)
	{
		temp.dataframe <- temp.dataframe[which(temp.dataframe$mapq >= qual),]
	}
	if(rmdup == TRUE)
	{
		temp.dataframe <- temp.dataframe[!duplicated(paste(temp.dataframe$rname, temp.dataframe$pos)),]
	}

	chrState<- round((table(temp.dataframe$strand)[1]-table(temp.dataframe$strand)[2])/nrow(temp.dataframe), digits=3)

	temp.dataframe <- temp.dataframe[which(temp.dataframe$rname == as.character(chr) & temp.dataframe$pos >= startLoc-padding & temp.dataframe$pos <= endLoc+padding),]

	bedDataframe <- temp.dataframe[order(temp.dataframe$pos),]
	return(list(bedDataframe, chrState))
}