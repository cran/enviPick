\name{mzagglom}
\alias{mzagglom}

\title{
Agglomerative partitioning of raw LC-HRMS measurements
}

\description{
	Agglomerative partitioning of LC-HRMS measurements.
	Preparatory step for \code{\link[enviPick]{mzclust}} and \code{\link[enviPick]{mzpick}}. 
	Requires an MSlist initilialized by \code{\link[enviPick]{readMSdata}} as input.  
}

\usage{
	mzagglom(MSlist, dmzgap = 10,  ppm = TRUE, drtgap = 500, minpeak = 4, 
	maxint=1E7, progbar=FALSE)
}

\arguments{  
	\item{MSlist}{MSlist generated by \code{\link[enviPick]{readMSdata}}}
	\item{dmzgap}{m/z gap width for partitioning}
	\item{ppm}{\code{dmzgap} given in ppm (TRUE) or as absolute value (FALSE)?}
	\item{drtgap}{RT gap width for partitioning}
	\item{minpeak}{Minimum number of measurements in a partition}
	\item{maxint}{Measurements equal or above this intensity will be retained even if ranging below \code{minpeak}}
	\item{progbar}{For debugging, ignore}
}



\details{
Partitioning of the full set of measurements into subsets is necessary to speed up the clustering procedure of \code{\link[enviPick]{mzclust}}. 
To this end, an agglomerative partitioning approach is used, combining measurements that are linked by values smaller than \code{drtgap} and 
\code{dmzgap} into single subsets. No measurements of two different subsets can be closer than \code{drtgap} and \code{dmzgap} to each other.
}

\section{Imbecile}{
Do not set \code{minpeak} bigger than its counterpart in \code{\link[enviPick]{mzclust}} or \code{\link[enviPick]{mzpick}}. 
Too complicated? Then rather use \code{\link[enviPick]{enviPickwrap}} for adjusting all function arguments.
}

\value{
	Returns the argument MSlist, with entries made:

\item{Parameters}{MSlist[[2]]: saves the parameter settings.}
\item{Scans}{MSlist[[4]]: matrix with raw measurements and tags resorted for partitions.}
\item{Partition_Index}{MSlist[[5]]: Index assigning partitions to sections in the raw measurement of MSlist[[4]]; required for fast (random) access.}

}

\author{Martin  Loos}

\section{Warning}{
Despite optimized code, this function has a potential to run for a intolerable long time or out of memory if (a) the parameters are set wrongly, (b) the .mzML/.mzXML-file was not centroided or
(c) the underlying data is inadequate for this peak picker. 
With regards to (a), do not assume gaps being larger than actually present. Instead, use \code{\link[enviPick]{plotMSlist}} to have a look at your 
data contained in MSlist after upload with \code{\link[enviPick]{readMSdata}}.
}

\seealso{\code{\link[enviPick]{mzclust}},\code{\link[enviPick]{mzpart}}}








