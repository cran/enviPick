\name{enviPickwrap}
\alias{enviPickwrap}

\title{
Peak detection for LC-HRMS .mzML and .mzXML data
}

\description{
A wrapper combining data upload, partitioning, EIC clustering and EIC peak detection from  \code{\link[enviPick]{readMSdata}}, 
\code{\link[enviPick]{mzagglom}}, \code{\link[enviPick]{mzclust}} and \code{\link[enviPick]{mzpick}}.
}

\usage{
enviPickwrap(filepath.mzXML, MSlevel=c(1), dmzgap=15, dmzdens=4, ppm=TRUE, 
drtgap=500, drtsmall=20, drtdens=250, drtfill=10, drttotal=200, minpeak=4,  
recurs=3, weight=2, SB=3, SN=2, minint=1E5, maxint=1E7, ended=2, progbar=FALSE)
}

\arguments{

	\item{filepath.mzXML}{Path to the .mzXML file to be read, \code{\link[enviPick]{readMSdata}}}  
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
	\item{maxint}{Peaks above this intensity are alwas included, regardless of checks for \code{SB} or \code{SN},\code{\link[enviPick]{mzpick}}}
	\item{ended}{Within the peak detection recursion set by argument \code{recurs}, how often can a peak detection fail, to end the recursion? \code{\link[enviPick]{mzpick}}}
	\item{progbar}{Show a progress bar (TRUE or FALSE)? May only work under Windows OS.}
  
}

\details{
For further details on the parameter settings, please refer to the arguments of the underlying functions \code{\link[enviPick]{readMSdata}}, 
\code{\link[enviPick]{mzagglom}}, \code{\link[enviPick]{mzclust}} and \code{\link[enviPick]{mzpick}}. 
}

\value{
MSlist
\item{State}{MSlist[[1]]: tags the individual steps the MSlist has undergone for peak picking.}
\item{Parameters}{MSlist[[2]]: saves the parameter settings.}
\item{Results}{MSlist[[3]]: saves the result summary values}
\item{Scans}{MSlist[[4]]: matrix with raw measurements (m/z, intensity, RT) and tags for partitions, EIC cluster or individual peaks.}
\item{Partition_Index}{MSlist[[5]]: Index assigning partitions to sections in the raw measurment of MSlist[[4]]. Required for fast (random) access, e.g., plotting.}
\item{EIC_index}{MSlist[[6]]: Index assigning EIC clusters to sections in the raw measurment of MSlist[[4]]. Required for fast access.}
\item{Peak_index}{MSlist[[7]]: Index assigning picked peaks to sections in the raw measurment of MSlist[[4]]. Required for fast access.}
\item{Peaklist}{MSlist[[8]]: Final peak list, cp. \code{\link[enviPick]{mzpick}}}
}

\author{Martin Loos}

\seealso{
	\code{\link[enviPick]{readMSdata}} 
	\code{\link[enviPick]{mzagglom}} 
	\code{\link[enviPick]{mzclust}} 
	\code{\link[enviPick]{mzpick}}
	\code{\link[enviPick]{plotMSlist}}	
}
