\name{plotdrugColorCorrection}
\alias{plotdrugColorCorrection}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Plot the result of drug color correction.
}
\description{
Three types of plots can be generated using the function specified by the \code{plot.type}. A plot of the estimated absorbance measures is produced when \code{plot.type = "raw"}. Finally, when \code{plot.type = "output"} the drug specific absorbance associated with each concentration is plotted against the concentration.
}
\usage{
plotdrugColorCorrection(A.data, file = file.path(getwd(), "drugColorCorrection.pdf"), 
                         drugs = NULL, pdfit = FALSE, width = 7, height = 3.5, 
                         type.plot = c("output", "raw"), xlab = NULL, ylab = NULL,
                         main = NULL, xlim = NULL, ylim = NULL, 
                         col = c("#71965A", "#4F6E9F", "#9F9692", "#9D2441",
                                 "#333333", "#662D91", "#71DEC0", "#F7931E"),
                         col.plate = TRUE, set.par = FALSE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{A.data}{
An A.data object created by the function \code{drugColorCorrection}.
}
  \item{file}{
Character specifying the pdf file for saving the resulting plots.
}
  \item{drugs}{
Character vector specifying the drugs for which the plot should be made. Defeaults to \code{NULL} where plots for all drugs are created.
}
  \item{pdfit}{
Logical indicating whether or not the plots should be saved in a pdf. 
}
  \item{width}{
Numeric value specifying the width of the pdf in inches.
}
  \item{height}{
Numeric value specifying the height of the pdf in inches.
}
  \item{type.plot}{
Character specifying the type of plot to be made. Currently the implemented types are 1) \code{"output"} which produces a plot of the  estimated drug specific absorbances and 2 ) \code{"raw"} which produces a plot of the normalised absorbace measures.
}
  \item{xlab}{
Character specifying the the x-axis label.
}
  \item{ylab}{
Character specifying the the y-axis label. For \code{plot.type = "output"} this should be a vector of length three.
}
  \item{main}{
The main title of the plots. 
}
  \item{xlim}{
Numeric vector for the range of the x-axis.
}
  \item{ylim}{
Numeric vector for the range of the y-axis. 
}
  \item{col}{
Vector of colours. 
}
  \item{col.plate}{
With \code{type.plot = c("raw")} should the colour set per plate. Defeault = \code{TRUE}.
}
 \item{set.par}{
Logical indicating whether or not options to par should be set by the function.
} 
 \item{\dots}{
Additional arguments passed along to plot. 
}
}
\value{
plot of residual absorbance values after pre-processing. The lefthand plot is for the model that does not account for the heteroschedastic variance whereas the model for the righthand plot does.
}
\references{
Steffen Falgreen et al. Exposure time independent summary statistics for assessment of drug dependent cell line growth inhibition (2013)
}
\author{
The function was written at department of haematology, Aalborg University Hospital and maintained by Steffen Falgreen.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{plotbgModel}}
}

\examples{
require(DoseR)
data(A.data)

## Plot of the drug specific Absorbance measures with an 95\% confidence interval.
plotdrugColorCorrection(A.data=A.data, type.plot = c("output"))

## Plot of the absorbance values. Use this plot for dagnostics of each dose response 
## experiment. Plates that are not fitted adequately should be removed a plate is found 
## In the plot to the left, the raw data is shown, the middle plot shows the absorbance values
## corrected in the conventional way, and the right plot shows the model corrected absorbance values.

par(mfrow = c(2,3))
plotdrugColorCorrection(A.data    = A.data, 
                        drugs     = "Doxorubicin",
                        set.par   = FALSE,
                        type.plot = c("raw"))
    
plotdrugColorCorrection(A.data    = A.data, 
                        drugs     = "Rituximab",
                        set.par   = FALSE,
                        type.plot = c("raw"))

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{~diagnostics}
\keyword{~pre-processing}