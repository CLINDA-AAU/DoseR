\name{CI}
\alias{CI}
\alias{CI.summary.DRdata}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Print information regarding absorbance data object.
}
\description{
Prints information regarding drugs and cell lines to be analysed 
}
\usage{
CI(A.data, ...)

\method{CI}{summary.DRdata}(A.data, ..., model = "G",
   types =  c("GI50", "TGI", "LC50", "AUC"),
   conf  = c(2.5, 97.5),
   time  = "48",
   dose.scale  = "mol/l",
   dose.logfun = "log10")
}
\arguments{
\item{A.data}{A.data object created by the function summary.DRdata.}
\item{...}{Further arguments passed to other methods.}
\item{model}{character, the dose response model to get the summary satistices from.}
\item{types}{vector of characters with the names of the summary statistics.}
\item{conf}{Numeric vector of legth two with the limits of the confidence interval.}
\item{time}{The time used for the summary statistic LC50 for the G-model. E.g. setting \code{t=48} indicates that the LC50 for the G model is actually the LC48.}
\item{dose.scale}{Character indicating the scale used for concentrations when estimating the isotonic regression and the summary statistcs. The unit is written as e.g. ug/ml to indicate micro grams per milli litre and defaults to \code{mol/l}. The unit is written as ug/ml The current implementations for multiples are:
 \tabular{lccccccccccc}{
    Name  \tab  \tab deca  \tab hecto \tab kilo \tab mega \tab giga \tab tera \tab peta \tab exa \tab zetta \tab yotta \cr
    Prefix  \tab  \tab da  \tab h \tab k \tab M \tab G \tab T \tab P \tab E \tab Z \tab Y \cr
    Factor \tab 10^0 \tab 10^1  \tab 10^2 \tab 10^3 \tab 10^6 \tab 10^9 \tab 10^12 \tab 10^15 \tab 10^18 \tab 10^21 \tab 10^24 
  }
  The current implementations for fractions are: 
\tabular{lcccccccccc}{
    Name   \tab deci  \tab centi \tab milli \tab micro \tab nano \tab pico \tab femto \tab atto \tab zepto \tab yocto \cr
    Prefix  \tab d  \tab c \tab m \tab u \tab n \tab p \tab f \tab a \tab z \tab y \cr
    Factor  \tab 10^-1  \tab 10^-2 \tab 10^-3 \tab 10^-6 \tab 10^-9 \tab 10^-12 \tab 10^-15 \tab 10^-18 \tab 10^-21 \tab 10^-24 
  }
}
\item{dose.logfun}{Character indicating if the concentrations should be log transformed. The possible inputs are \code{nolog} for non log transformation, \code{log10, log2, log} for log transforming with base 10, 2, and e, respectively.}
}
\value{
Confidence intervals for the chosen summary statistics and drugs. \code{T0} corresponds to the estimated doubling time. 
}
\references{
Steffen Falgreen et al. Exposure time independent summary statistics for assessment of drug dependent cell line growth inhibition (2013)
}
\author{
The function was written at department of haematology, Aalborg University Hospital and maintained by Steffen Falgreen.
}
\note{
When the dose response models have been estimated it is possible to make various plots of the data. The estimated summary statistics are obtainable using the function CI. See the examples for further usage of the A.data object.
}
%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{pdfCI}},\code{\link{tableCI}},\code{\link{DRdataBoxplot}},\code{\link{plot.DRdata}},\code{\link{plot.growthModel}},\code{\link{plotGrid}}
}
\examples{
require(DoseR)
## load Dose Response data
data(A.data)
CI(A.data)

}
