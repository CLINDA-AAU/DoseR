\name{DRdataBoxplot}
\alias{DRdataBoxplot}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Box plots of summary statistics
}
\description{
Box plots of chosen summary statistics according to the specified dose reponse model.
}
\usage{
DRdataBoxplot(A.data, splitvar = NULL, splitvar.col = NULL, col.all = NULL, 
              model = "G", type = "AUC", names = NULL, drug = 1, time = NULL, 
              main = NULL, ylab = NULL, clean.names = NULL, plot.order = NULL, 
              dose.scale = "mol/l", dose.logfun = "log10", ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{A.data}{
An A.data object created by the function \code{doseResponseModel}.
}
  \item{splitvar}{
Character specifying a variable according to which the boxplots should be separated e.g. a column specifying the disease origin of the cell lines.
}
  \item{splitvar.col}{
Specify the colour used for each element of splitvar.
}
  \item{col.all}{
Vector of colours. The colours are used from left to right.
}
  \item{model}{
Character specifying which dose response model the boxplot should be produced for. Currently the follwing models are implemented \code{G}, \code{Gres}, \code{R}, \code{D}, \code{DG}, and \code{RG}. Only models estimated by the function \code{doseResponseModel} can be used.
}
  \item{type}{
Character indicating the summary statistic used in the plot. This can be either of \code{GI50}, \code{TGI}, \code{LC50}, or \code{AUC}.
}
  \item{names}{
Character vector specifying the names of the cell lines that are to be plotted.
}
  \item{drug}{
Character specifying the drug to be plotted.
}
  \item{time}{
character or numeric indicator for which time point the plot is based upon. Not applicable for the \code{model = "G"}-model.
}
  \item{main}{
The main title for the plot
}
  \item{ylab}{
The label for the y-axis
}
  \item{clean.names}{
character specifying a repetitive part of cell line names that should be remove. e.g. " - inHS" if the cell lines are named as "cell.line.name - inHS".
}
  \item{plot.order}{
character vector of cell line names specifying the order of the boxplots.
}
  \item{dose.scale}{
Character indicating the scale used for concentrations when estimating the isotonic regression and the summary statistcs. The unit is written as e.g. ug/ml to indicate micro grams per milli litre and defaults to \code{mol/l}. The unit is written as ug/ml The current implementations for multiples are:
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
  \item{\dots}{
Further arguments passed along to boxplot.
}
}
\value{
Box plot of bootstrapped summary statistics.
}
\references{
Steffen Falgreen et al. Exposure time independent summary statistics for assessment of drug dependent cell line growth inhibition (2013)
}
\author{
The function was written at department of haematology, Aalborg University Hospital and maintained by Steffen Falgreen.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{CI}},\code{\link{DRdataBoxplot}},\code{\link{plot.DRdata}},\code{\link{plot.growthModel}},\code{\link{plotGrid}}
}

\examples{
require(DoseR)
data(A.data)

DRdataBoxplot(A.data, type = "AUC", model = "G", 
         splitvar = "disease", col.all = c("#71965A", "#4F6E9F"),
         drug = "Rituximab",
         time = 48)

DRdataBoxplot(A.data, type = "AUC", model = "G", 
         splitvar = "disease", col.all = c("#71965A", "#4F6E9F"),
         drug = "Doxorubicin",
         time = 48)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{~boxplot}
