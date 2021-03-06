% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plotROIFrequencies.R
\name{plotROIFrequencies}
\alias{plotROIFrequencies}
\title{Generates plot of wcRatios}
\usage{
plotROIFrequencies(fileLocation, fileName, fileFrequencies, chr, binSize,
  startLoc, endLoc, ROIname = "wholeChr", gapfile = 0,
  callingThreshold = 0.66, ROI = 0, padding = 0, strand = TRUE,
  png = TRUE)
}
\arguments{
\item{fileLocation}{Specifies directory to save file}

\item{fileName}{Unique identifier of input file}

\item{fileFrequencies}{Dataframe of wcRatio used to generate histrogram}

\item{chr}{Genomic coordinate chromosome}

\item{binSize}{The number of reads in each bin used to calculate wcRatio}

\item{startLoc}{Genomic coordinate start position}

\item{endLoc}{Genomic coordinate end position}

\item{ROIname}{Genomic coordinate number}

\item{gapfile}{Input file of gaps to plot}

\item{callingThreshold}{Threshold line used to find ROIs}

\item{ROI}{ROI location table}

\item{padding}{Number of bases to extend beyond genomic coordinates listed in regionTable  e.g. [10000]}

\item{strand}{if\code{TRUE} plot reads above histrogram}

\item{png}{if\code{TRUE} Save png}
}
\description{
Generates single plot of wcRatios for fileFrequencies
}
\author{
Ashley D. Sanders, Mark Hills
}

