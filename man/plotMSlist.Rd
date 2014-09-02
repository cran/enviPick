\name{plotMSlist}
\alias{plotMSlist}

\title{Interactive raw data and result viewer}

\description{
View your centroided LC-HRMS data and partitioning / clustering / peak-picking results. 
Monitor what peak-picking produces and if it fails; get a grip on optimal parameter settings from comparison with underlying raw data.
}

\usage{
plotMSlist(MSlist, RTlimit = FALSE, mzlimit = FALSE, shiny = FALSE, ppmbar = 8)
}

\arguments{
  \item{MSlist}{An MSlist returned by \code{\link[enviPick]{readMSdata}}, \code{\link[enviPick]{mzpart}}, \code{\link[enviPick]{mzclust}} or \code{\link[enviPick]{mzpick}}}
  \item{RTlimit}{Initialize plot: two-element vector of lower and upper RT plot limits. Set to FALSE to view full data.}
  \item{mzlimit}{Initialize plot: two-element vector of lower and upper m/z limits. Set to FALSE to view full data.}
  \item{shiny}{For debugging. Ignore.}
  \item{ppmbar}{Size of m/z bar (in ppm) shown at large zoom}  
}

\details{For more help, use the help button in the interactive plot. Based on low-level R plot functionality!}

\author{Martin Loos}

\note{\code{MSlist} may contain a lot of data; rendering of measurements for plotting may thus decrease in speed when zooming out or using the full-view mode.}

\seealso{
	\code{\link[enviPick]{readMSdata}}
	\code{\link[enviPick]{mzagglom}}
	\code{\link[enviPick]{mzclust}}	
	\code{\link[enviPick]{mzpick}}	
}

