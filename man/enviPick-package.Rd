\name{enviPick-package}
\alias{enviPick-package}
\alias{enviPick}
\docType{package}
\title{
Peak picking for high resolution liquid chromatography - mass spectrometry (LC-HRMS) data
}
\description{
Peak picking for centroided .mzXML high-resolution LC-MS data. 
Built on a three-step approach of (1) data partitioning, (2) unsupervised clustering of extracted ion chromatograms (EIC) of and 
(3) shape-independent peak detection wihtin individual EICs. Interactive plot access to all results and the underlying raw measurements.
Browser UI for non-R users. Batch processing.
}


\details{
\tabular{ll}{
Package: \tab enviPick\cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 2014-14-07\cr
License: \tab GPL-2\cr
}
After initial upload of an .mzXML standard format file with \code{\link[enviPick]{readMSdata}}, above steps (1) to (3) are 
calculated by \code{\link[enviPick]{mzagglom}}, \code{\link[enviPick]{mzclust}} and \code{\link[enviPick]{mzpick}}, respectively. 
The wrapper for joint upload and processing is \code{\link[enviPick]{enviPickwrap}}.
The raw data and the results of each step, a so-called MSlist object, can be viewed by \code{\link[enviPick]{plotMSlist}}, producing an 
interactive plot that conveniently offers zoom, drag and select functionality and easy identification of individual partitions, EIC cluster or
peaks. Batch processing can be done via \code{\link[enviPick]{enviPickbatch}}. 
To export a peak list from an MSlist object, use \code{\link[enviPick]{writePeaklist}}.

For converting LC-HRMS measurement files from various vendor formats or .mzML to centroided .mzXML we strongly recommend the MSConvert tool from
ProteoWizard; for centroidization choose Filters -> Peak Picking -> Prefer Vendor -> Add.
}

\section{Warning}{
This package has only been tested on HIGH-RESOLUTION Thermo Orbitrap and QExactive measurements processed (centroided) with ProteoWizard`s MSConvert. 
It may not give satisfying results for chromatograms affected by mass shifts from centroid-centroid interferences prevalent at low resolutions.
}

\note{
In the package context, peak picking refers to extracting individual ion chromatograms (EICs) from centroided data and identifying peaks 
in these EICs. In the ProteoWizard MSConvert context, peak picking refers to identifying individual peaks within single HRMS scans, alias centroidization. 
}

\author{
Martin Loos
Maintainer: Martin Loos <Martin.Loos@eawag.ch>
}

\references{
Loos, M. (XXXX). Extraction of ion chromatograms by unsupervised clustering of high-resolution mass spectrometry data. Some Journal. Sometime.

ProteoWizard: Open Source Software for Rapid Proteomics Tools Development Darren Kessner; Matt Chambers; Robert Burke; David Agus; Parag Mallick Bioinformatics 2008
\url{http://proteowizard.sourceforge.net/}
}

\seealso{

	\code{\link[enviPick]{readMSdata}}
	\code{\link[enviPick]{mzagglom}}
	\code{\link[enviPick]{mzpart}}
	\code{\link[enviPick]{mzclust}}
	\code{\link[enviPick]{mzpick}}
	\code{\link[enviPick]{plotMSlist}}
	\code{\link[enviPick]{enviPickwrap}}	
	\code{\link[enviPick]{enviPickbatch}}	
	
}

\examples{
\dontrun{

##################################################
# (1) Define path to an LC-HRMS .mzML file (not provided with package):

filepath.mzML<-"C:/.../2012_07_01.mzML"

# (2) Initialize an MSlist object and load this .mzML file into it:

MSlist<-readMSdata(filepath.mzML, MSlevel=c(1))

# (3) Partition the measurements now available in MSlist:

MSlist<-mzagglom(MSlist,dmzgap=10,ppm=TRUE,drtgap=500,minpeak=4,maxint=1E7)

# (4) EIC clustering of the partitions now available in MSlist:

MSlist<-mzclust(MSlist,dmzdens=5,ppm=TRUE,drtdens=120,minpeak=4)      

# (5) Peak picking within the EIC clusters now available in MSlist:

MSlist<-mzpick(MSlist, minpeak = 4, drtsmall = 50, drtfill = 10,  drttotal = 200, recurs = 4, 
	weight = 2, SB = 3, SN=2, minint = 1E4, maxint = 1e+07, ended = 2)

# (6) Export a peak list now available in MSlist:

writePeaklist(MSlist,"directory","filename")

# (7) View your partitioning / EIC clustering / peak picking results:

plotMSlist(MSlist,ppmbar=10);
##################################################

##################################################
# Do above steps (1) to (5) in one wrap, then export a peak list:

MSlist<-enviPickwrap( filepath.mzML,
                      MSlevel=c(1),
                      dmzgap=10,  
                      dmzdens=5,       
                      ppm=TRUE,
                      drtgap=1000, 
                      drtsmall=20,
                      drtdens=250,
                      drtfill=10,
                      drttotal=200,
                      minpeak=4,
                      recurs=10,
                      weight=2,
                      SB=3,
                      SN=2,
                      minint=10E4,
                      maxint=10E6,
                      ended=2,
                      progbar=TRUE
                       )
					   
writePeaklist(MSlist,"directory","filename")
##################################################

}
}

\keyword{ package }
