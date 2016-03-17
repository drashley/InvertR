#' Calculate wcRatios for strandTable
#'
#' Calculate wcRatios for input bedfile using a read-based sliding window of checkNum size
#'
#' @param strandTable a dataframe object of sequencing reads
#' @param checkNum Window size used to calculate wcRatios
#' @param verbose Verbose messages if \code{TRUE}.
#' @author Ashley D. Sanders, Mark Hills
#' @export
countFreqs <- function(strandTable, checkNum=10, verbose=TRUE)
{
	if(verbose){message(paste('-> Creating WC frequencies for ', index, " [", indexCounter, "/", bamFileLength, "]", sep=""))}
	actualRow <- checkNum-1
	if (nrow(strandTable) <= actualRow)
	{
		message("Not enough reads for BINSIZE! consider changing criteria")
		return(0)
	} else {
		freqFrame <- matrix(nrow=nrow(strandTable), ncol=1)
		#calculate the average ratio of W and C across the entire region. This allows plotting ratio to default to 1 irrespecitve of template
		firstRead <- abs(as.numeric(table(strandTable[1:nrow(strandTable),3]))/nrow(strandTable))
		firstRead <- firstRead[1]

		if(firstRead > 0.5)
		{
			for(row in seq(1:(nrow(strandTable)-actualRow)))
			{
				#finds the numerical difference between number of + and - between the row and row+9
				frequency <- abs(as.numeric(table(strandTable[row:(row+actualRow),3]))/checkNum)
				frequency <- frequency[1]
				freqFrame[row,] <- frequency
			}

			for( row in rev(seq((nrow(strandTable)-checkNum+1), nrow(strandTable))) )
			{
				frequency <- abs(as.numeric(table(strandTable[row:(row+actualRow),3]))/checkNum)
				frequency <- frequency[1]
				freqFrame[row,] <- frequency
			}

		} else {
			for(row in seq(1:(nrow(strandTable)-actualRow)))
			{
				#finds the numerical difference between number of + and - between the row and row+9
				frequency <- 1-abs(as.numeric(table(strandTable[row:(row+actualRow),3]))/checkNum)	
				frequency <- frequency[1]				
				freqFrame[row,] <- frequency
			}

			for( row in rev(seq((nrow(strandTable)-checkNum+1), nrow(strandTable))) )
			{
				frequency <- 1-abs(as.numeric(table(strandTable[row:(row+actualRow),3]))/checkNum)
				frequency <- frequency[1]
				freqFrame[row,] <- frequency
			}					
		}

		colnames(freqFrame) <- "WCratio"	
		strandTable <- cbind(strandTable[1:nrow(strandTable),], freqFrame)
		strandTable[(nrow(strandTable)-checkNum):nrow(strandTable), 5] <-1 # replaces last WCratios to 1
		return(strandTable)
	}
}