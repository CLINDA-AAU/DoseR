\name{tableCI}
\alias{tableCI}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
LaTeX code for table of summary statistics.
}
\description{
Creates LaTeX code for table of summary statistics.
}
\usage{
tableCI(A.data, 
           model       = "G",
           dec         = 2,
           conf        = c(2.5, 97.5), 
           types       = c("T0", "GI50", "TGI","LCt", "AUCq"),
           type.order  = "GI50",
           dose.scale  = "mol/l", 
           dose.logfun = "log10",
           drugs       = 1,
           splitvar    = NULL)
}
\arguments{
\item{A.data}{A.data object created by the function summary.DRdata.}
\item{model}{character, the dose response model to get the summary satistices from.}
\item{dec}{Numeric indicator for the number of decimals.}
\item{conf}{Numeric vector of legth two with the limits of the confidence interval.}
\item{types}{vector of characters with the names of the summary statistics.}
\item{type.order}{Character indicating the summary statistic used for sorting.}
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
\item{drugs}{Character vector, names of the drugs the table should be based on.}
\item{splitvar}{Character naming a variable in the dataset by which the cell lines should be grouped.}
}
\value{
Confidence intervals for the chosen summary statistics and drugs. \code{T0} corresponds to the estimated doubling time. The results will be converted to a tex table and compiled to a pdf.
}
\references{
Steffen Falgreen et al. Exposure time independent summary statistics for assessment of drug dependent cell line growth inhibition (2013)
}
\author{
The function was written at department of haematology, Aalborg University Hospital and maintained by Steffen Falgreen.
}
%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{CI}},\code{\link{pdfCI}},\code{\link{DRdataBoxplot}},\code{\link{plot.DRdata}},\code{\link{plot.growthModel}},\code{\link{plotGrid}}
}
\examples{
require(DoseR)
## load Dose Response data
data(A.data)
tableCI(A.data)
}
