\name{doseResponseModel}
\alias{doseResponseModel}
\alias{doseResponseModelShiny}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Calculation of dose response models.
}
\description{
This function is used to estimate the dose response models: \code{G}: the G model is calculated, \code{R}: the relative model is calculated, \code{D}: the difference based model is calculated, \code{DG}: the G model is calculated based on the D model, \code{RG} the G model is calculated based on the R model. When the models are calculated isotonic regression is used to fit the dose response curves while taking the singulatity of the \code{G} and \code{D} models into account. Finally, the summary statistics \code{GI50}, \code{TGI}, \code{LC50}, and \code{AUC0} is calculated.
}
\usage{
doseResponseModel(A.data, models = c("G", "R", "D", "DG", "RG"), 
                  dose.scale = "mol/l", dose.logfun.from = "nolog", 
                  dose.logfun.to = "log10", parametrisation = "restricted", 
                  AUC.q = 0, t = 48, doublingvar = NULL, cut = 0.025, 
                  update = TRUE, verbose = FALSE, 
                  progressbar = "text", save = TRUE, 
                  shiny.input = NULL, session = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{A.data}{
An A.data object created by the function \code{createMetaData}.
}
  \item{models}{
Character vector with indicating what models should be calculated. Currently the follwing models are implemented \code{G}: the G model is calculated, \code{R}: the relative model is calculated, \code{D}: the difference based model is calculated, \code{DG}: the G model is calculated based on the D model, \code{RG} the G model is calculated based on the R model.

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
  \item{dose.logfun.from}{
Character indicating if the concentrations given in the protocols are log transformed. The possible inputs are \code{nolog} for non logtransformed concentrations, \code{log10, log2, log} for log transformed with base 10, 2, and e, respectively.
}
  \item{dose.logfun.to}{
Character indicating if the concentrations should be log transformed. The possible inputs are \code{nolog} for non log transformation, \code{log10, log2, log} for log transforming with base 10, 2, and e, respectively.
}
  \item{parametrisation}{
Character indicating whether or not the estimated growth of the treated cells should be restricted such that it is always slower than untreated. Can be either \code{restricted} or \code{unrestricted}.
}
  \item{AUC.q}{
Numeric value determining the minimum for the summary statistic AUC.q. When zero this is the area under the dose response curve above zero. Defeaults to 0.
}
  \item{t}{
Numeric value indicating the time spand LC50 should be based upon for the G model. e.g. when \code{t = 48} the LC50 values corresponds to the concentration where the halving time is 48 hours. 
}
  \item{doublingvar}{
Character giving the name of the column contaning doubling time for each cell line (Optional). When only one time point is used for a dose response experiment it is possible to convert the R model to the G model if the doubling time for each cell line is supplied. 
}
  \item{cut}{
Numeric value determining the lowest allowable value for absorbance value defeaults to 0.025. Since all absorbance values must be positive all negative ones need to be converted to positive values. If the value is set to 0.025 all values below 0.025 are truncated at 0.025.
}
  \item{update}{
Logical value. Should the analysis be updated or run from scratch. Defaults to \code{TRUE}.
}
  \item{verbose}{
Logical value. Should the function indicate what experiment it is currently reading. Defaults to \code{FALSE}. 
}
  \item{progressbar}{
The type of progress bar used to show how far along the function is. Can be either \code{"window"}, \code{text} or \code{none}.
}
  \item{save}{
Should the data be saved.
}
  \item{shiny.input}{
Used for the shiny server.
}
  \item{session}{
Used for the shiny server.
}
}
\value{
The ouput of the function is an A.data object of class \code{bgModel}. This is a list with the following components:
 \item{meta.list}{This is a list of meta data objects.}
 \item{call}{A list containing information regarding the call to the function.}
 \item{auxiliary}{List of auxiliary data used by other functions.}
 \item{data}{List of data frames. The added elements are described below.}
 \item{drug.color.correct}{Contains the results of the fitted dose reponse experiments for colour correction.}
 \item{fits}{List of the fitted objects. The fitted objects for model corrections is stored in element \code{bgModel} The fitted objects for the G-model is stored in element \code{growthModel}}
 \item{summary}{Summary statistics for the dose response data. The data consists of nested lists as:
 summary$drug$model$summarystatistic}
  
 
 The element data is expanded with the following elements:
 \item{GM.mean}{The estimated G model for the dose response data including all bootstrapped data}
 \item{T0.list}{list of \code{data.frame}s. For each drug the estimated doubling times are stored in a \code{data.frame} named after the drug}
 \item{GI.mean}{The estimated D and R models for the dose response data including all bootstrapped data}
 \item{DR.data}{Combination of \code{GI.mean} and \code{GM.mean}.}
 \item{iso.fits}{List with the fitted curves for each dose response experiment. The data inculding the fitted curves are stored as: drug$model. The data is available for all bootstrapped datasets.}

See the examples for usage of the A.data object. 
}
\references{
Steffen Falgreen et al. Exposure time independent summary statistics for assessment of drug dependent cell line growth inhibition (2013)
}
\author{
The function was written at department of haematology, Aalborg University Hospital and maintained by Steffen Falgreen.
}
\note{
When the dose response models have been estimated it is possible to make various plots of the data. The estimated summary statistics are obtainable using the function CI. See the examples for further usage of the A.data object.
}
%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{CI}},\code{\link{DRdataBoxplot}},\code{\link{plot.DRdata}},\code{\link{plot.growthModel}},\code{\link{plotGrid}}
}
\examples{
require(DoseR)
data(A.data)
# Not run
if(FALSE)
A.data <- doseResponseModel(A.data = A.data,
                            models  = c("G", "R", "RG",  "D", "DG"),
                            t      = 48,
                            AUC.q  = 0,
                            dose.scale = "mol/l",
                            dose.logfun.to = "log10",
                            parametrisation = "restricted",
                            cut = 0.001,   
                            verbose = FALSE,
                            update = FALSE)

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

# plot a grid of the dose response curves with the bootstrapped results shown 
plotGrid(A.data=A.data, barcol="#33333350", 
         ncol = 2, drug = "Doxorubicin")


# plot of the dose response curve based on the D-model for the different time points
par(mfrow = c(1,2))
plot.DRdata(x = A.data,#xlim = c(-7.5, -4.8),
            model = "D",
            col.scheme =  c("#71965A", "#4F6E9F"),
             drug = "Doxorubicin",legend = TRUE,
            n.columns = 1, plot.data = TRUE,
            legend.cex = 1, times = c(12, 24, 36, 48))

plot.DRdata(x = A.data,#xlim = c(-7.5, -4.8),
            model = "D",
            col.scheme =  c("#71965A", "#4F6E9F"),
             drug = "Rituximab",legend = TRUE,
            n.columns = 1, plot.data = TRUE,
            legend.cex = 1, times = c(12, 24, 36, 48))

par(mfrow = c(2,2))

plot.DRdata(x = A.data,
            model = "G", 
            drug = "Rituximab",
            col.scheme =  c("#71965A", "#4F6E9F"),
            n.columns = 1,
            plot.data = TRUE,legend = TRUE,
            legend.cex = 1)

plot.DRdata(x = A.data,
            model = "G", 
            drug = "Doxorubicin",
            col.scheme =  c("#71965A", "#4F6E9F"),
            n.columns = 1,
            plot.data = TRUE,legend = TRUE,
            legend.cex = 1)


DRdataBoxplot(A.data, type = "AUC", model = "G", 
         splitvar = "disease", col.all = c("#71965A", "#4F6E9F"),
         drug = "Rituximab",
         time = 48)

DRdataBoxplot(A.data, type = "AUC", model = "G", 
         splitvar = "disease", col.all = c("#71965A", "#4F6E9F"),
         drug = "Doxorubicin",
         time = 48)
         
         
# calculate summary statistics based on the G-model
CI(A.data, model="G")
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
