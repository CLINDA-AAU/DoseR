\name{plotbgModelresid}
\alias{plotbgModelresid}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Function for plotting the residuals of the fitted absorbance model. 
}
\description{
This function is used to test that the fit of the normalisation model is adequate and that the heteroschedastic variance is fitted correctly.
}
\usage{
plotbgModelresid(A.data, names = NULL, times = NULL, drugs = NULL, 
                  pdfit = FALSE, file = file.path(getwd(), "resid.plots.pdf"), 
                  pointsize = 12, xlab = "Fitted values", ylab = "Residuals", 
                  main = NULL, col = 1:8, width = 3.34, height = 2.23097112)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{A.data}{
An A.data object created by either the function \code{bgModel}.
}
  \item{names}{
Character vector with the names of the cell lines for which the plot should be made. Defeaults to NULL where a plot is generated for all cell lines.
}
  \item{times}{
Numeric vector with the time.points for which the plot should be made. Defeaults to NULL where a plot is generated for all time points.
}
  \item{drugs}{
Character vector with the drugs for which the plot should be made. Defeaults to NULL where a plot is generated for all drugs.
}
  \item{pdfit}{
Logical. if \code{TRUE} all plots are saved in a pdf with name given in the \code{file} argument.
}
  \item{file}{
Character giving the name of the created pdf file.
}
  \item{pointsize}{
The pointsize used in the created pdf.
}
  \item{xlab}{
Character specifying the name for the x-axis.
}
  \item{ylab}{
Character specifying the name for the x-axis.
}
  \item{main}{
Character specifying the title for the plot.
}
  \item{col}{
Vector with colours for each experiment.  
}
  \item{width}{
Numeric value specifying the width of the pdf.
}
  \item{height}{
Numeric value specifying the width of the pdf.
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

# Plot the residuals to check whether any of the fits are 
# very poor. In order to do the plot for the cell line RPMI-8226 
# at the 48 hour time.point for the drug Doxorubidin the 
# following code is used:

plotbgModelresid(A.data, 
                 names = "RPMI-8226",
                 times = 48,
                 drugs = "Doxorubicin",
                  col = c("#71965A", "#4F6E9F", "#9F9692", "#9D2441",
                         "#333333", "#662D91", "#71DEC0", "#F7931E"))
      
# The following code chunk saves the residual plot for all cell lines,
# all drugs, and all times in the pdf file "resid.plots.pdf".                     
                         
# plotbgModelresid(A.data, names = NULL,
#                 times = NULL,
#                 drugs = NULL,
#                 pdfit = TRUE,
#                 file=  "resid.plots.pdf",
#                 col = c("#71965A", "#4F6E9F", "#9F9692", "#9D2441",
#                         "#333333", "#662D91", "#71DEC0", "#F7931E"),
#                 pointsize = 6,
#                 width = 13/2.54, height =(13/2.54)/1.5)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{~diagnostics}
\keyword{~pre-processing}
\keyword{~absorbance}
\keyword{~Residual}
