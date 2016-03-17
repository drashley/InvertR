#' Collapse ratioFile and generate bedGraph 
#'
#' Collapse space between values in ratioFile that are consecutive, with bedGraph format option
#' 
#' @param ratioFile A dataframe with wcRatio coordinates.
#' @param fileName Identifier used in bedgraph header
#' @param verbose Verbose messages if\code{TRUE} 
#' @param asBedgraph if\code{TRUE} Formates output file as bedgraph
#' @author Ashley D. Sanders, Mark Hills
#' @export
collapseStrands <- function(ratioFile, fileName, verbose=TRUE, asBedgraph=FALSE)
{
 	options(scipen=999)
 	startSize <- nrow(ratioFile)
	ratioFile$tmp <- c(3, ratioFile[1:nrow(ratioFile)-1,5])
	ratioFile <- ratioFile[which(ratioFile[,5] != ratioFile$tmp),]
	ratioFile$tmp <- NULL


	if(asBedgraph)
	{
		bedGraphHead <- data.frame(chromosome='track type=bedGraph', start=paste('name=', fileName, '_', ratioFile[1,1],sep=""), end=paste('description=', fileName, sep=""), ratio='visibility=full color=0,0,0 autoScale=off graphType=points viewLimits=0:100', stringsAsFactors=FALSE)
		outputData <- data.frame(chromosome=ratioFile[,1], start=ratioFile[,2], end=c(ratioFile[2:nrow(ratioFile),2]-1, tail(ratioFile$pos,1)), ratio=ratioFile[,5]*100, stringsAsFactors=FALSE)
		outputData<- outputData[which(outputData$start <= outputData$end),] # corrects for any occurances where end is less than start (errors on UCSC Browser)
		outputData <- rbind(bedGraphHead, outputData)
	} else {

		outputData <- data.frame(chromosome=ratioFile[,1], start=ratioFile[,2], end=c(ratioFile[2:nrow(ratioFile),2]-1, tail(ratioFile$pos,1)), ratio=ratioFile[,5])
	}
	if(verbose){message(paste("Bed file collapsed from ", startSize , " lines to ", nrow(outputData) , sep=""))}
	return(outputData)
}
