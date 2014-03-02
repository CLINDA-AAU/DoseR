\name{bgModel}
\alias{bgModel}
\alias{bgModelShiny}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Normalisation of dose response data according to the model  \eqn{\delta\alpha + \beta|^{2\xi}}.
}
\description{
This function normalises the data according to the model \eqn{\delta\alpha + \beta|^{2\xi}}, where \eqn{\delta} is multiplicative error associated with each experiment, \eqn{\beta} is the plate specific background, and finally \eqn{\alpha} is the mean absorbance. 
}

\usage{
bgModel(A.data, 
        update          = TRUE,
        parametrisation = "unrestricted", 
        outlier.test    = 3,
        outlier.iter    = 2,
        weights         = "fitted", 
        fitted.a        = FALSE,
        varpower.min    = 10e-4,
        varpower.iter   = 50,
        contr           = c("sum", "helmert", "treatment"),
        progressbar     = "text",
        verbose         = FALSE,
        save            = TRUE,
        shiny.input     = NULL,
        session         = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{A.data}{
An A.data object created by either the function \code{readDBFData} or \code{drugColorCorrection}.
}
  \item{update}{
Should the analysis be updated or run from scratch. When set to \code{TRUE} only new experiments will be normalised. When \code{FALSE} The normalisation is conducted for all experiments.
}
  \item{parametrisation}{
The parametrisation used for fitting the model can be either of \code{unrestricted} or \code{restricted}. when set to \code{unrestricted} there is no restriction on the fitted absorbance values. On the contrary, if it set equal to \code{restricted} the absorbance is restricted to be positive. 
}
 \item{outlier.test}{
Numeric value indicating the number of standard deviations an obserbance measure should be from the estimated mean in order to be deemed an outlier. Defaults to 3.
}
  \item{outlier.iter}{
Numeric value indicating the number of iterations the outlier.test should be run. Defaults to \code{2}.
}
  \item{weights}{
  Either \code{fitted} or NULL. The variance of the absorbance measures are usually of a heteroschedastic form e.g. \eqn{|\delta\alpha + \beta|^{2\xi}\sigma^2}. If the value of \code{weights} is not equal to NULL the power \eqn{\xi} is fitted.
}
  \item{fitted.a}{
 If FALSE the heteroschedastic variance is fitted according \eqn{|\delta\alpha + \beta|^{2\xi}\sigma_{\epsilon}^2} and if TRUE it is fitted according to \eqn{|\alpha|^{2\xi}\sigma^2}
}
   \item{varpower.min}{
Numeric value indicating the convergence criterion for fitted weight for the variance. Defeaults to \code{0.001}
}
  \item{varpower.iter}{
Numeric value determining the number of allowed iterations for fitting the weight.
}
  \item{contr}{
  The contrast used when fitting the normalisation model. This can be be either of \code{sum}, \code{helmert}, or \code{treatment}. We strongly urge the use of the sum contrast which is used as defeault. 
}
  \item{progressbar}{
The type of progress bar used to show how far along the function is. Can be either \code{"window"}, \code{text} or \code{none}.
}
  \item{verbose}{
Should the function indicate what experiment it is currently reading. Defaults to \code{FALSE}.
}
  \item{save}{
Should the results saved (TRUE) or not (FALSE). Internally used for bootstrapping and should be left as TRUE.
}
  \item{shiny.input}{
Used for the shiny server.
}
  \item{session}{
Used for the shiny server.
}
}
\details{
Dose response experiments are usually done in triplicates to gain accuracy in the resulting absorbance measures. Normally one simply subtracts the background. However this does not take into account various errors that occur during seeding of the cell lines. This function normalises the data according to the model \eqn{\delta\alpha + \beta|^{2\xi}}, where \eqn{\delta} is multiplicative error associated with each experiment, \eqn{\beta} is the plate specific background, and finally \eqn{\alpha} is the mean absorbance. 
}
\value{
The ouput of the function is an A.data object of class \code{bgModel}. This is a list with the following components:
 \item{meta.list}{This is a list of meta data objects.}
 \item{call}{A list containing information regarding the call to the function.}
 \item{auxiliary}{List of auxiliary data used by other functions.}
 \item{data}{List of data frames. The normalised data for all wells are stored in the element \code{bc.data}. The estimated absorbance values are stored in the element \code{bc.mean}.}
 \item{drug.color.correct}{Contains the results of the fitted dose reponse experiments for colour correction.}
 \item{fits}{List of the fitted objects. The fitted objects for model corrections is stored in element \code{bgModel}}
 
The cotribution to the abovementioned data from this function is:

1) \code{data.frame} \code{bc.data} within the element \code{data} which contains:
 \item{absorbance}{The colour corrected absorbance values.}
 \item{absorbance.nc}{The original absorbation values.}
 \item{para}{The parametrisation used for fitting the model.}
 \item{NB}{The colour corrected absorbance values subtracted the plate specific background.}
 \item{outlier}{Was the well deemed an outlier.}
 \item{BC2}{The model corrected absorbance values corrected for drug colour, background, and variation between cell seeding.}
 \item{BC}{The estimated absorbance value for. This is the value used in further analysis.}
 \item{Std.Error}{The estimated standard deviation of \code{BC}}
 \item{lo}{The lower limit of the 95\% CI for \code{BC}}
 \item{up}{The upper limit of the 95\% CI for \code{BC}}

2) \code{data.frame} \code{bc.data} within the element \code{data} which contains:
 \item{absorbance}{The mean of colour corrected absorbance values.}
 \item{para}{The parametrisation used for fitting the model.}
 \item{NB}{The mean of colour corrected absorbance values subtracted the plate specific background.}
 \item{BC}{The estimated absorbance value for. This is the value used in further analysis.}
 \item{Std.Error}{The estimated standard deviation of \code{BC}}
 \item{lo}{The lower limit of the 95\% CI for \code{BC}}
 \item{up}{The upper limit of the 95\% CI for \code{BC}}
 
3) The list of fitted \code{gnls} objects.
 
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

## perfor the model based pre-processing:
A.data <- bgModel(A.data          = A.data, 
                  weights         = "fitted", 
                  fitted.a        = FALSE,
                  contr           = c("sum"),
                  outlier.test    = 3,
                  outlier.iter    = 2,
                  varpower.min    = 10e-4,
                  varpower.iter   = 50,
                  parametrisation = "unrestricted",
                  save            = FALSE,
                  update          = TRUE)

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
\keyword{pre-processing}
\keyword{normalisation}% __ONLY ONE__ keyword per line
\keyword{background}
\keyword{correction}
