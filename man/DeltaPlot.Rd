\name{DeltaPlot}
\alias{DeltaPlot}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Plot of the dose response curves with each cell line presented in a panel.
}
\description{
Plot of the dose response curves with each cell line presented in a panel.
}
\usage{
DeltaPlot(A.data, 
          diff.var    = 1, 
          drug        = 1, 
          model       = "G", 
          type        = "AUC", 
          time        = "48",
          names       = NULL,
          dose.scale  = "mol/l",
          dose.logfun = "log10",
          col         = "#4F6E9F",
          reverse     = FALSE, 
          las         = 2, 
          ylab        = NULL, 
          main        = "Difference", 
          ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{A.data}{
The A.data object created by the function bgModel.
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
  \item{names}{
  The names of the cell lines that are to plotted.
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
  \item{col}{
Character vector of colours. 
}
  \item{reverse}{
Logical indicater for whether or not the the number of plots per figure. 
}
  \item{las}{
Logical indicater for whether or not the the number of plots per figure. 
}
  \item{ylab}{
Logical indicater for whether or not the the number of plots per figure. 
}
  \item{main}{
Logical indicater for whether or not the the number of plots per figure. 
}
  \item{\dots}{
Further arguments passed on to the plot. 
}
}
\value{
plot of absorbance values plotted aginst time with the result of the G-model superimposed.
}
\references{
Steffen Falgreen et al. Exposure time independent summary statistics for assessment of drug dependent cell line growth inhibition (2013)
}
\author{
The function was written at department of haematology, Aalborg University Hospital and maintained by Steffen Falgreen.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{plotbgModelresid}}
}

\examples{
require(DoseR)
data(A.data)


}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{~diagnostics}
\keyword{~pre-processing}
\keyword{~absorbance}
