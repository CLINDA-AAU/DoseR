\name{plotbgModel}
\alias{plotbgModel}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Plot of the normalised absorbance values.
}
\description{
Plot of the normalised absorbance values. 
}
\usage{
plotbgModel(A.data, names = NULL, drugs = NULL, times = NULL, 
            pdfit = FALSE, unit = "mol/l", logfun = "log10", 
            plots = c("Raw", "Color Corrected", 
            "Simple Correction", "Model Correction"), 
            ylab = "Absorbance", xlab = NULL, main = NULL, 
            title.1 = NULL, 
            title.2 = c("Raw", "Color Corrected", 
            "Simple Correction", "Model Correction"), 
            figure.output = file.path(getwd(), "Normalisation"), 
            col = c("#4F6E9F", "#71965A", "#9F9692", "#9D2441", 
                    "#333333", "#662D91", "#71DEC0", "#F7931E"), 
                    width = 6.69, height = NULL, pointsize = 5, 
                    cex = 1, set.par = TRUE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{A.data}{
The A.data object created by the function bgModel.
}
  \item{names}{
Character vector of the cell lines that are to be plotted.
}
  \item{drugs}{
Character vector of the drugs that are to be plotted.
}
  \item{times}{
Numeric vector of the times that are to be plotted.
}
  \item{pdfit}{
Logical indicating whether or not the plots should be saved in a pdf.
}
  \item{unit}{Character indicating the scale used for concentrations when estimating the isotonic regression and the summary statistcs. The unit is written as e.g. ug/ml to indicate micro grams per milli litre and defaults to \code{mol/l}. The unit is written as ug/ml The current implementations for multiples are:
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
  \item{logfun}{
Character indicating if the concentrations should be log transformed. The possible inputs are \code{nolog} for non log transformation, \code{log10, log2, log} for log transforming with base 10, 2, and e, respectively.
}
  \item{plots}{
character vector of the plots made. The value \code{Raw} plots of the raw absorbance values. The value \code{Color Corrected} plots the color corrected absorbance values. The value \code{Simple Correction} is the color and background corrected. The value \code{Model Correction} plots the model corrected absorbance values. Defeaults to c("Raw", "Color Corrected", "Simple Correction", "Model Correction").
}
  \item{ylab}{
Label for the y-axis.
}
  \item{xlab}{
Vector of labels for the x-axis
}
  \item{main}{
Vector of titles for he plot. 
}
  \item{title.1}{
%%     ~~Describe \code{title.1} here~~
}
  \item{title.2}{
\code{c("Raw", "Color Corrected", "Simple Correction", "Model Correction")}
}
  \item{figure.output}{
Directory for storing the pdf files.
}
  \item{col}{
Vector of colors used for the absorbance values. 
}
  \item{width}{
Numeric value indicating the width of the pdf.
}
  \item{height}{
Numeric value indicating the height of the pdf.
}
  \item{pointsize}{
Numeric value indicating the pointsize used in the pdf.
}
  \item{cex}{
Numeric value indicating the cex value used in the plot.
}
  \item{set.par}{
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

# Diagnostics plot of for cell line "RPMI-8226" for the drug 
# RPMI-8226" at time point 48. From left to right the plots are 
# for Raw data, colour corrected, simple background correction, 
# and model-based correction.

plotbgModel(A.data  = A.data, 
            pdfit   = FALSE,
            names   = "RPMI-8226",
            drugs   = "Doxorubicin",
            times   = 48,
            set.par = TRUE)

# The following code chunk saves the diagnostics plot for all cell lines,
# all drugs, and all times in the folder normalisation.

# plotbgModel(A.data=A.data, 
#             figure.output="Normalisation", 
#             pdfit = TRUE, 
#             pointsize = 8, 
#             set.par=TRUE)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{~diagnostics}
\keyword{~pre-processing}
\keyword{~absorbance}
