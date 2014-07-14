\name{enviPickbatch}
\alias{enviPickbatch}

\title{
Batch peak picking for several .mzXML input files
}

\description{
Given a folder of .mzXML input files and subsequent processing, .csv peak tables with picked peaks are written to an output folder
}

\usage{
enviPickbatch(folderin, folderout, MSlevel=c(1), dmzgap=15, dmzdens=4, 
ppm=TRUE, drtgap=500, drtsmall=20, drtdens=250, drtfill=10, drttotal=200, 
minpeak=4, recurs=10, weight=2, SB=3, SN=2, minint=1E5, maxint=1E7, 
ended=2, progbar=FALSE)
}

\arguments{
	\item{folderin}{Folder containing .mzXML input files}
	\item{folderout}{Destination folder for .csv peak tables}
	\item{MSlevel}{numeric 1 (MS) or 2 (MS-MS), \code{\link[enviPick]{readMSdata}}}
	\item{dmzgap}{m/z gap width for partitioning, \code{\link[enviPick]{mzpart}}}
	\item{dmzdens}{Maximum measurement deviation (+/-) of m/z within an EIC, \code{\link[enviPick]{mzclust}}}
	\item{ppm}{\code{dmzgap},  given in ppm (TRUE) or as absolute value (FALSE)}
	\item{drtgap}{RT gap width for partitioning, \code{\link[enviPick]{mzpart}}}  
    \item{drtsmall}{Peak definition - RT window of a peak; cp. \code{minpeak} argument}   
	\item{drtdens}{Maximum length of EICs, \code{\link[enviPick]{mzclust}}} 
	\item{drtfill}{RT gap length to be filled, \code{\link[enviPick]{mzpick}}}   
	\item{drttotal}{Maximum RT length of a single peak, \code{\link[enviPick]{mzpick}}}
	\item{minpeak}{Peak definition - minimum number of measurements per peak to found in windows of \code{drtsmall} and \code{dmzsmall}, \code{\link[enviPick]{mzpick}}} 
	\item{recurs}{Maximum number of peaks within one EIC, \code{\link[enviPick]{mzpick}}}
	\item{weight}{Weight for assigning measurements to a peak, \code{\link[enviPick]{mzpick}}}   
	\item{SB}{Signal-to-base ratio, \code{\link[enviPick]{mzpick}}}
	\item{SN}{Signal-to-noise ratio, \code{\link[enviPick]{mzpick}}}   
	\item{minint}{Minimum intensity of a peak, \code{\link[enviPick]{mzpick}}}
	\item{maxint}{Peaks above this intensity are alwas included, regardless of checks for \code{SB} or \code{SN}, \code{\link[enviPick]{mzpick}}}
	\item{ended}{Within the peak detection recursion set by argument \code{recurs}, how often can a peak detection fail, to end the recursion? \code{\link[enviPick]{mzpick}}}
	\item{progbar}{Show a progress bar (TRUE or FALSE)? May only work under Windows OS.}
  }
  
  
\details{
For further details on the parameter settings, please refer to the arguments of the underlying functions \code{\link[enviPick]{readMSdata}}, 
\code{\link[enviPick]{mzagglom}}, \code{\link[enviPick]{mzclust}} and \code{\link[enviPick]{mzpick}}. 
}

\value{
Folder with .csv peak tables, each containing columns with:
		m/z (mean m/z of peak measurements), 
		var_m/z (m/z variation of peak measurements),
		max_int (base-line corrected maximum intensity),
		sum_int (sum of all base-line corrected peak measurement intensities),
		RT (retention time at maximum intensity),
		minRT (start peak RT),
		maxRT (end peak RT),
		peak# (peak ID number),
		EIC# (EIC ID number),
		Score (not yet implemented)
}

\author{Martin  Loos}

\seealso{
	\code{\link[enviPick]{readMSdata}} 
	\code{\link[enviPick]{mzagglom}} 
	\code{\link[enviPick]{mzclust}}
	\code{\link[enviPick]{mzpick}}
	\code{\link[enviPick]{plotMSlist}}	
}
