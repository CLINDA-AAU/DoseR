\name{DoseR-package}
\alias{DoseR-package}
\alias{DoseR}
\docType{package}
\title{
Analysis of dose response experiments
}
\description{
Analysis of dose response experiments using the procedure described in Falgreen et. al (2013). 
Data analysis of dose response experiments. The data created during the analysis is saved in an RData file such that cell line data already analysed is not analysed again. The analysis is completed using the following steps:

\enumerate{
    \item \code{createMetaData}
    \item \code{readDBFData}
    \item \code{drugColorCorrection}
    \item \code{bgModel}
    \item \code{bootstrap}
    \item \code{DoseRModel}
  }
}


\details{
\tabular{ll}{
Package: \tab DoseR\cr
Type: \tab Package\cr
Version: \tab 0.1\cr
Date: \tab 2013-10-21\cr
License: \tab What license is it under?\cr
}
}
\author{
The package was written at Department of haematology, Aalborg University Hospital.

Maintainer: Steffen Falgreen Larsen to <sfl@rn.dk>
}
\references{
~~ Literature or other references for background information ~~
}
\keyword{Dose Response Experiments }

\examples{
## See the Sweave document for an example of the workflow.
}
