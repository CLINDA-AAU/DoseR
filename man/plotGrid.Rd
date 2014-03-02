\name{plotGrid}
\alias{plotGrid}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Plot of all dose reponse curves for a single drug in a panel.
}
\description{
The investigation of dose response curves individually is much easier if they are plotted in individual panels with bootstrapped confidence interval shown.
}
\usage{
plotGrid(A.data = A.data, model = "G", type = "AUC", times = NULL, 
         dose.scale = "mol/l", dose.logfun = "log10", drug = 1, 
         name = 2, names = NULL, plot.order = NULL, conc.names = NULL, 
         ylim = NULL, xlim = NULL, nrows = NULL, ncols = NULL, ylab = NULL, 
         xlab = NULL, main = paste("Dose Response Curves for", drug), 
         absorbance.CI = FALSE, absorbance.CI.col = "#C6C6C5", 
         absorbance.CI.alpha = 80, bootstrap.conf = TRUE, 
         barcol = "#71965A", bar.height = 1.3, 
         plotgrid = TRUE, grid.col = "#C6C6C5", grid.lty = 1, grid.lwd = 1, 
         bs.col = c("#333333", "#9D2441"), bs.lty = 1, bs.lwd = 0.5, 
         bs.alpha = 50, line.col = c("#333333", "#9D2441"), 
         line.lty = rep(1, 8), line.lwd = rep(1, 8), line.alpha = "", 
         col.by.identifier = TRUE, col.points = c("#333333", "#9D2441"), 
         pch = 1, plot.data = FALSE, log = "", pdfit = FALSE,
         pdf.width = 6.6929, pdf.height = 6.6929, pointsize = 8)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{A.data}{
An A.data object created by the function \code{doseResponseModel}.
}
  \item{model}{
Character specifying which dose response model the boxplot should be produced for. Currently the follwing models are implemented \code{G}, \code{Gres}, \code{R}, \code{D}, \code{DG}, and \code{RG}. Only models estimated by the function \code{doseResponseModel} can be used.
}
  \item{type}{
Character indicating the summary statistic used to sort the pqnels. This can be either of \code{GI50}, \code{TGI}, \code{LC50}, or \code{AUC}.
}
  \item{times}{
Numeric indicator of the time points that are to be plotted.
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
  \item{drug}{
Character vector of the drugs that are to be plotted.
}
  \item{name}{
When set to 2, only the the name part of the cell line is used to separate the plots. 
}
  \item{names}{
The names of the cell lines to be plotted.
}
  \item{plot.order}{
The order the panels are produced. This overwrights the argument set by \code{type}.
}
  \item{conc.names}{
Names of the concentrations used in the plot. This is the names specified in the protocols.
}
  \item{ylim}{
Limit of the y-axis.
}
  \item{xlim}{
Limit of the x-axis.
}
  \item{nrows}{
The number of rows.
}
  \item{ncols}{
The number of cols.
}
  \item{ylab}{
The name of the x-axis.
}
  \item{xlab}{
The name of the x-axis.
}
  \item{main}{
The main title of the plot.
}
  \item{absorbance.CI}{
Should 95\% confidence intervals for the absorbance values be shown.
}
  \item{absorbance.CI.col}{
must be specified as Hex rgb, i.e. something like "#334455".
}
  \item{absorbance.CI.alpha}{
Numeric value between 10 and 99 specifying the opacity of the confidence interval.
}
  \item{bootstrap.conf}{
Should the bootstrapped curves be shown.
}
  \item{barcol}{
Colour of the bar above each plot.
}
  \item{bar.height}{
The height of the abovementioned bar.
}
  \item{plotgrid}{
Should a grid be plotted. 
}
  \item{grid.col}{
The colour of the grid.
}
  \item{grid.lty}{
Line type of the grid.
}
  \item{grid.lwd}{
Line width og the grid.
}
  \item{bs.col}{
colour of bootstrapped curves must be specified as Hex rgb, i.e. something like "#334455".
}
  \item{bs.lty}{
Line type for the bootstrapped curves.
}
  \item{bs.lwd}{
Line width of the bootstrapped curves
}
  \item{bs.alpha}{
Numeric value between 10 and 99 specifying the opacity of the bootstrap curves.
}
 \item{line.col}{
Colour of the the curves must be specified as Hex rgb, i.e. something like "#334455".
}
  \item{line.lty}{
Line type for the curves.
}
  \item{line.lwd}{
Line width for  the curves.
}
  \item{line.alpha}{
Numeric value between 10 and 99 specifying the opacity of the curves.
}
  \item{col.by.identifier}{
Should the data be colourd in accordace with the specifier.
}
   \item{col.points}{
Colour of the points.
}
  \item{pch}{
The point type for the data.
}
  \item{plot.data}{
Should the data be plotted.
}
  \item{log}{
Argument passed along to the plot. if \code{log = "y"} the y-axis is log-tranformed.
}
  \item{pdfit}{
Should the created plots be save in pdfs.
}
  \item{pdf.width}{
Width of the pdf.
}
  \item{pdf.height}{
Height of te pdf.
}
  \item{pointsize}{
Pointsize of the pdf.
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
\code{\link{CI}},\code{\link{DRdataBoxplot}},\code{\link{plot.DRdata}},\code{\link{plot.growthModel}}
}

\examples{
require(DoseR)
data(A.data)

# plot a grid of the dose response curves with the bootstrapped results shown 
plotGrid(A.data=A.data, barcol="#33333350",   
         ncol = 2, drug = "Doxorubicin")
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{~dose}
\keyword{~response}
\keyword{~curves}