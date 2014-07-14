\name{readMSdata}
\alias{readMSdata}

\title{
Read LC-HRMS measurement .mzXML data
}

\description{
Initiates an MSlist object and reads LC-HRMS measurement data from .mzXML files into it. 
}

\usage{
readMSdata(filepath.mzXML, MSlevel=c(1), progbar=FALSE, minRT=FALSE, 
maxRT=FALSE, minmz=FALSE, maxmz=FALSE)
}

\arguments{
  \item{filepath.mzXML}{Path to the .mzXML file to be read}
  \item{MSlevel}{numeric 1 (MS) or 2 (MS-MS)}
  \item{progbar}{Debugging, ignore. Show a progress bar (TRUE or FALSE)? Might only work in Windows OS}
  \item{minRT}{Filter for measurements with retention time >= minRT. Otherwise set to FALSE.}
  \item{maxRT}{Filter for measurements with retention time <= maxRT. Otherwise set to FALSE.}
  \item{minmz}{Filter for measurements with m/z >= minmz. Otherwise set to FALSE.}
  \item{maxmz}{Filter for measurements with m/z <= maxmz. Otherwise set to FALSE.}  
}

\details{
The return value, a so-called MSlist object, is a simple R list object that contains (a) the raw measurement data, (b) intermediate/final 
results of the peak picking procedure and (c) indices for random access, to be passed among functions. 
Peaks are nested in EIC clusters which in turn are nested in partitions which in turn are subsets of measurements; 
MSlist[[4]] is resorted accordingly during all peak picking steps.

Setting \code{minRT}, \code{maxRT}, \code{minmz} or \code{maxmz} allows you to filter your .mzML data. 
On the one hand, this may be very useful if only being interested in certain ranges of an experiment.
On the other hand, this allows you to upload subset data of an experiment too large to be loaded into R at once.
}

\value{
MSlist
\item{State}{MSlist[[1]]: tags the individual steps the MSlist has undergone.}
\item{Parameters}{MSlist[[2]]: saves parameter settings.}
\item{Results}{MSlist[[3]]: saves a result summary.}
\item{Scans}{MSlist[[4]]: matrix with raw measurements (m/z, intensity, RT) and tags for partitions, EIC cluster and individual peaks.}
\item{Partition_Index}{MSlist[[5]]: Index assigning partitions to sections in the raw measurment of MSlist[[4]]. Needed for fast (random) access during, e.g., plotting.}
\item{EIC_index}{MSlist[[6]]: Index assigning EIC clusters to sections in the raw measurment of MSlist[[4]]. Required for fast access.}
\item{Peak_index}{MSlist[[7]]: Index assigning picked peaks to sections in the raw measurment of MSlist[[4]]. Required for fast access.}
\item{Peaklist}{MSlist[[8]]: Final peak list, cp. \code{\link[enviPick]{mzpick}}}
}

\section{Warning}{
Use \code{\link[enviPick]{plotMSlist}} to check your data in MSlist for consistency at an early stage before further processing.
}

\author{
Martin Loos
}

\note{
It is your responsibility to ensure your input files are centroided. If not, R may freeze and the peak picker will not return valid results.
}

\seealso{
\code{\link[enviPick]{mzagglom}}
}

