\name{plot.growthModel}
\alias{plot.growthModel}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Plot of the growth curves for the dose response experiments.
}
\description{
This function plots the absorbance values against time. The bootstrapped curves can be plotted together with the curves. Since a plot like this shold be made for each cell line for each drug it is possible to save the plots in pdf files for later inspection. The folder for the pdfs is specified in the argument \code{figure.output}.
}
\usage{
\method{plot}{growthModel}(x, ..., time.points.used = c("all"), 
     drugs = NULL, names = NULL, conc.names = NULL, 
     ylim = NULL, xlim = NULL, nrows = NULL, 
     ncols = NULL, xlab = "Time (Hours)", ylab = "Absorbance", 
     main = "Dose Response Growth Curves", 
     absorbance.CI = FALSE, absorbance.CI.col = "#C6C6C5", 
     absorbance.CI.alpha = 80, bootstrap.conf = TRUE, 
     barcol = "grey", bar.height = 1.3, 
     plotgrid = FALSE, grid.col = "#C6C6C5", 
     grid.lty = 1, grid.lwd = 1, bs.col = "#C6C6C5", 
     bs.lty = 1, bs.lwd = 0.5, bs.alpha = 80, 
     line.col = "#662D91", line.lty = 1, 
     line.lwd = 1, line.alpha = "", plot.data = TRUE, 
     col.by.identifier = TRUE, col.points = "#9F9692", 
     pch.points.outlier = 4, pch.points = 1, plot.all = TRUE, 
     plot.data.all = FALSE, line.col.all = "#9F9692", line.lty.all = 1, 
     line.lwd.all = 0.8, line.col.C0 = "#71965A", line.col.GI50 = "#4F6E9F", 
     line.col.TGI = "#9D2441", line.col.LC48 = "#333333", log = "", 
     pdfit = FALSE, figure.output = getwd(), pdf.width = 6.6929, pdf.height = 6.6929, pointsize = 8)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
An A.data object created by the function \code{doseResponseModel}.
}
  \item{\dots}{
Further arguments passed to plot.
}
  \item{time.points.used}{
The G model have been calculated based on the all and the individual time points. Here you can specify which of the established models should be plotted. defeaults to \code{time.points.used = "all"} indicating that the fit based on all time points is plotted.
}
  \item{drugs}{
Character vector of the drugs for which the plot should be made.
}
  \item{names}{
Character vector of cell lines for which the plot should be made.
}
  \item{conc.names}{
Names of the concentrations used in the plot. These are the names specified in the protocols.
}
  \item{ylim}{
Limit of the y-axis.
}
  \item{xlim}{
Limit of the x-axis.
}
  \item{nrows}{
Number of rows.
}
  \item{ncols}{
Number of columns
}
  \item{xlab}{
Name of the x-axis.
}
  \item{ylab}{
Name of the y-axis.
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
  \item{plot.data}{
Should the data be plotted.
}
  \item{col.by.identifier}{
Should the data be colourd in accordace with the specifier.
}
  \item{col.points}{
Colour of the points.
}
  \item{pch.points.outlier}{
pch used to indicate outliers.
}
  \item{pch.points}{
pch used to indicate non-outliers.
}
  \item{plot.all}{
Should the data be combined into a single plot for all concentrations. 
}
  \item{plot.data.all}{
Should the abovementioned panel be plotted with data.
}
  \item{line.col.all}{
Colour of the lines in the all data panel.
}
  \item{line.lty.all}{
Line type for the lines in the all data panel.
}
  \item{line.lwd.all}{
Line width for the lines in the all data panel.
}
  \item{line.col.C0}{
Colour specifying the untreated control.
}
  \item{line.col.GI50}{
Colour specifying the GI50 value.
}
  \item{line.col.TGI}{
Colour specifying the TGI value.
}
  \item{line.col.LC48}{
Colour specifying the LC48 value.
}
  \item{log}{
Argument passed along to the plot. if \code{log = "y"} the y-axis is log-tranformed.
}
  \item{pdfit}{
Should the created plots be save in pdfs.
}
  \item{figure.output}{
The folder where the plots should be stored
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
\code{\link{CI}},\code{\link{DRdataBoxplot}},\code{\link{plot.DRdata}},\code{\link{plotGrid}}
}

\examples{
require(DoseR)
data(A.data)
# plot of the fitted models with the bootstrapped results shown.
# This plot can be used for investigation of the  fitted model.
plot.growthModel(x = A.data, 
                 drugs = "Rituximab",
                 names = "OCI-Ly7",
                 time.points.used=c("all", 48),
                 bootstrap.conf=FALSE,
                 line.col= c("#662D91",  "#F7931E"),
                 absorbance.CI=FALSE,
                 conc.names = paste("C", c(0, 1:17), sep = ""),
                 pointsize = 8,
                 pdfit = FALSE,
                 plotgrid = TRUE,  nrows=4, ncols=5,
                 col.by.identifier = FALSE)

plot.growthModel(x = A.data, 
                 drugs = "Doxorubicin",
                 names = "SU-DHL-4",
                 time.points.used=c("all", 48),
                 bootstrap.conf=TRUE,
                 line.col= c("#662D91",  "#F7931E"),
                 absorbance.CI=FALSE,
                 conc.names = paste("C", c(0, 1:17), sep = ""),
                 pointsize = 8,
                   pdfit = FALSE,
                 plotgrid = TRUE,  nrows=4, ncols=5,
                 col.by.identifier = FALSE)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{~dose}
\keyword{~response}
\keyword{~curves}
