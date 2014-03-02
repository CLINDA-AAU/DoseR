\name{drugColorCorrection}
\alias{drugColorCorrection}
\alias{drugColorCorrectionShiny}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Correction for absorbance caused by the drug.
}
\description{
This functionction is used to estimate the absorbance values caused by the drug. Next, the estimates are used to correct the absorbance measures for all experiments using the drug. 
}
\usage{
drugColorCorrection(A.data = A.data, weights = "fitted", fitted.a = FALSE, 
                    contr = c("sum", "helmert", "treatment"), outlier.test = 3, 
                    outlier.iter = 2, varpower.min = 0.001, varpower.iter = 50, 
                    parametrisation = "unrestricted", 
                    update = TRUE, save = TRUE, progressbar = "text", 
                    shiny.input = NULL, session = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{A.data}{
An A.data object created by the function \code{readDBFData}.
}
  \item{weights}{
  Either \code{fitted} or NULL. The variance of the absorbance measures are usually of a heteroschedastic form e.g. \eqn{|\delta\alpha + \beta|^{2\xi}\sigma^2}. If the value of \code{weights} is not equal to NULL the power \eqn{\xi} is fitted.
}
  \item{fitted.a}{
 If FALSE the heteroschedastic variance is fitted according \eqn{|\delta\alpha + \beta|^{2\xi}\sigma_{\epsilon}^2} and if TRUE it is fitted according to \eqn{|\alpha|^{2\xi}\sigma^2}
}
  \item{contr}{
  The contrast used when fitting the normalisation model. This can be be either of \code{sum}, \code{helmert}, or \code{treatment}. We strongly urge the use of the sum contrast which is used as defeault. 
}
  \item{outlier.test}{
Numeric value indicating the number of standard deviations an obserbance measure should be from the estimated mean in order to be deemed an outlier. Defaults to 3.
}
  \item{outlier.iter}{
Numeric value indicating the number of iterations the outlier.test should be run. Defaults to \code{2}.
}
  \item{varpower.min}{
Numeric value indicating the convergence criterion for fitted weight for the variance. Defeaults to \code{0.001}
}
  \item{varpower.iter}{
Numeric value determining the number of allowed iterations for fitting the weight.
}
  \item{parametrisation}{
The parametrisation used for fitting the model can be either of \code{none}, \code{square}, \code{abs}, \code{exp}, or \code{expinv}. when set to \code{none} there is no restriction on the fitted absorbance values. On the contrary, if any of the others are chosen the absorbance is restricted to be positive. It is also possible to give the arguments as a vector. Then the optimization routine will circle through each element in the vector until it converges.
}
  \item{update}{
Should the analysis be updated or run from scratch. When set to \code{TRUE} only new experiments will be normalised. When \code{FALSE} The normalisation is conducted for all experiments.
}
  \item{save}{
Should the data be saved. Defeaults to \code{TRUE}.
}
  \item{progressbar}{
The type of progress bar used to show how far along the function is. Can be either \code{"window"}, \code{text} or \code{none}.
}
  \item{shiny.input}{
Used for the shiny server.
}
  \item{session}{
Used for the shiny server.
}
}
\details{
This function is used to correct for absorbance measures caused by the drug itself. This is the case for drugs with colour, and drugs that interacts with the MTS assy. To gain precise estimates it is necessary to account for the elevated absorbance values when estimating the cell count. To use this function for the purpose it is necessary to setup a number of plates where each concentration used in the experiment is presented but no cells are added. The name used in place of a cell line name must be specified in the function \code{createMetaData} using the argument \code{correctionname}.
}
\value{
The ouput of the function is an A.data object of class \code{drugColorCorrection}. This is a list with the following component  
 \item{meta.list}{This is a list of meta data objects.}
 \item{call}{A list containing information regarding the call to the function.}
 \item{auxiliary}{List of auxiliary data used by other functions.}
 \item{data}{List of data frames. The colour corrected Absorbance data are stored in \code{raw.data} and the protocols are stored in \code{protocols}}
 \item{drug.color.correct}{Contains the results of the fitted dose reponse experiments for colour correction.}


In the \code{data.frame} \code{raw.data} the column \code{absorbance.nc} corresponds to the uncorrected absorbance measurements and the the column  \code{absorbance} contains the corrected ones.

The element \code{background.list} within \code{drug.color.correct} contains the fitted abosorbance values for the various concentrations:
 \item{Concentration}{The concentrations}
 \item{BC}{The absorbance associated with the drug at the concentration.}
 \item{lo}{The lower limit of the 95\% CI for \code{BC}}
 \item{up}{The upper limit of the 95\% CI for \code{BC}}
 \item{BC2}{The absorbance associated with the drug at each concentration in relation to the control.}
 \item{Conc}{The concentrations on a molar}
}
\references{
Steffen Falgreen et al. Exposure time independent summary statistics for assessment of drug dependent cell line growth inhibition (2013)
}
\author{
The function was written at department of haematology, Aalborg University Hospital and maintained by Steffen Falgreen.
}
\note{
After the data have been read into R the need to be preprocessed. If any drug colour correction plates have been established continue to the funktion \code{\link{drugColorCorrection}} otherwise continue to the function \code{\link{bgModel}}.
}
%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{createMetaData}},\code{\link{readDBFData}},\code{\link{bgModel}},\code{\link{plotdrugColorCorrection}}
}
\examples{
require(DoseR)
## load Dose Response data
data(A.data)

A.data <- 
      drugColorCorrection(A.data = A.data, 
                          progressbar = "text", save = FALSE)

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
\keyword{~colourcorrection}
\keyword{~colorcorrection}% __ONLY ONE__ keyword per line
