\name{DeltaCI}
\alias{DeltaCI}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Computation of confidence intervals for difference in summary statistics.
}
\description{
This function is used for calculation of confidence intervals for the difference in summary statistics when the same cell line have been treated with a drug in two conditions. 
}
\usage{
DeltaCI(A.data, 
        conf     = c(2.5, 97.5), 
        diff.var = 1, 
        drug     = 1, 
        model    = "G", 
        type     = "AUC", 
        time     = "48",
        dose.scale  = "mol/l",
        dose.logfun = "log10",
        reverse     = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{A.data}{
The A.data object created by the function \code{doseResponseModel}.
}
  \item{conf}{
Numeric vector indicating the limits of the confidence intercal defaults to \code{c(2.5, 97.5)}.
}
  \item{diff.var}{
Character indicating the variable used for separating the experiments for which the difference is to be calculated.
}
  \item{drug}{
The drug for which the difference should be calculated.
}
  \item{type}{
The summary statistic for which the difference is to be calculated.
}
  \item{time}{
  The time point where the difference is to be calculated.
}
  \item{dose.scale}{Character indicating the scale used for the summary statistcs. The unit is written as e.g. ug/ml to indicate micro grams per milli litre and defaults to \code{mol/l}. The unit is written as ug/ml The current implementations for multiples are:
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
  \item{dose.logfun}{
Character indicating if the concentrations should be log transformed. The possible inputs are \code{nolog} for non log transformation, \code{log10, log2, log} for log transforming with base 10, 2, and e, respectively.
}
  \item{reverse}{
Reverse the calculation of difference.
}
}
\value{
Confidence interval on differences for calculated summary statistics.
}
\references{
Steffen Falgreen et al. Exposure time independent summary statistics for assessment of drug dependent cell line growth inhibition (2013)
}
\author{
The function was written at department of haematology, Aalborg University Hospital and maintained by Steffen Falgreen.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{doseResponseModel}},\code{\link{DeltaPlot}}
}

\examples{
require(DoseR)
data(A.data)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{~CI}
\keyword{~difference}

