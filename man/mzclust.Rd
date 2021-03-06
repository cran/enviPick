\name{mzclust}
\alias{mzclust}

\title{
EIC-clustering of LC-HRMS data.
}

\description{
Based on the measurement partitions generated by \code{\link[enviPick]{mzagglom}}, 
extracted ion chromatograms (EICs) are assigned by a clustering procedure. Preparatory step for \code{\link[enviPick]{mzpick}}.  
}

\usage{
mzclust(MSlist,dmzdens=10,ppm=TRUE,drtdens=60,minpeak=4,maxint=1E6,
progbar=FALSE,merged=TRUE,from=FALSE,to=FALSE  )
}

\arguments{
  \item{MSlist}{MSlist returned by \code{\link[enviPick]{mzpart}}}
  \item{dmzdens}{Maximum measurement deviation (+/-) of m/z from its mean within each EIC}
  \item{ppm}{\code{dmzdens} given in ppm (TRUE) or as absolute value (FALSE)?}
  \item{drtdens}{Retention time (RT) tolerance for clustering; defined as (+/-) time units relative to the lowest and highest RT value in each cluster}
  \item{minpeak}{Minimum number of measurements expected in an EIC}
  \item{maxint}{EIC cluster with measurements above this intensity are kept, even if they do not fulfill \code{minpeak}}
  \item{progbar}{For debugging, ignore}
  \item{merged}{Merge EIC cluster of comparable m/z (TRUE or FALSE)?}
  \item{from}{For debugging, ignore}
  \item{to}{For debugging, ignore}
}

\details{
Within individual partitions calculated by \code{\link[enviPick]{mzagglom}}, an unsupervised clustering of measurements to individual ion chromatograms (EICs) is performed.
For this purpose, a first EIC cluster is initialized with the most intense measurement, given an m/z uncertainty of 2*\code{dmzdens}. 
Along decreasing intensities, all other measurements are then sequentially either assigned to this cluster or used to define new clusters. 
For assignment, measurements must range both within the current tolerances of \code{dmzdens} and \code{drtdens} of an existing cluster. 
If several cluster are eligible for assignment, the one with the smallest mass difference between measurement m/z and cluster mean m/z will be used. 
Each time a new assignment to an existing cluster is made, its m/z estimate can be improved, i.e., the \code{dmzdens} tolerance around its mean m/z gradually 
shrinks from 2*\code{dmzdens} to \code{dmzdens}. In addition, \code{dmzdens} is used to update the RT tolerance of a cluster at each assignment.

With no measurements left, EIC clusters nested in m/z are then merged, relative to the m/z boundaries of the most intense cluster and stepwise along increasing mean m/z differences.
Finally, EIC cluster are filtered to fulfill either \code{minpeak} or \code{maxint}.   
}

\value{
	Returns the argument MSlist, with entries made:
	
	\item{Parameters}{MSlist[[2]]: saves the parameter settings.}
	\item{Scans}{MSlist[[4]]: matrix with raw measurements and tags resorted for EIC clusters within the partition subsets.}
	\item{EIC_index}{MSlist[[6]]: Index assigning EIC cluster to sections in the raw measurement of MSlist[[4]]; required for fast (random) access.}
}

\references{
Loos, M. (XXXX). Extraction of ion chromatograms by unsupervised clustering of high-resolution mass spectrometry data. Some Journal. Sometime.
}

\author{Martin Loos}

\note{Too small values for \code{dmzdens} or too large \code{drtdens} may cause erratic EICs cluster}

\seealso{\code{\link[enviPick]{mzpick}}}


