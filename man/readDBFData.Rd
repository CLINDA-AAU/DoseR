\name{readDBFData}
\alias{readDBFData}
\alias{readDBFDataShiny}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Function for reading dbf files into R in accordance with specified protocols.
}
\description{
This function is used read dBase files into R according to a specified protocol.
}
\usage{
readDBFData(A.data, update = TRUE, well.id = "WELLNUM", absorbance.id  = "M1", 
            discard.lines  = c(1, 2), dosevar = "Concentration",
            additivevar = "Additive", controlval = "Control",
            backgroundval = "Background", mistakeval = "X", 
            remove.rows = c("A","H"), remove.cols = c(1, 12), 
            progressbar = "text", verbose= FALSE, save = TRUE, 
            shiny.input = NULL, session = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{A.data}{
An A.data object created by the function \code{createMetaData}.
}
  \item{update}{
Should the analysis be updated or run from scratch. Defaults to \code{TRUE}.
}
  \item{well.id}{
Name of the coloumn in the .dbf file containg well id's. Defeaults to "WELLNUM".
}
  \item{absorbance.id}{
Name of the coloumn in the .dbf file containg absorbance values. Defeaults to "M1".
}
\item{discard.lines}{
A vector of linenumber that are to be discarded from the .dbf fiel. Defeaults to c(1, 2)
}
 \item{dosevar}{
Name of the column that contains concentration values in protocol. 
}
  \item{additivevar}{
Name of the column that contains information regarding what each well  values in protocol. 
}
  \item{controlval}{
Value indicating that the cell line is untreated for the column defined in \code{additivevar}. Defeaults to \code{Control}
}
  \item{backgroundval}{
Value indicating that the well is used to measure background for the column defined in \code{additivevar}. Defeaults to \code{Background}
}
  \item{mistakeval}{
Value indicating that the particular well should not be used for analysis. This is specified in the setup part of the protocol where the default value \code{X} denote a mistake. 
}
  \item{remove.rows}{
Vector of row indices that should be discarded according to the plate. Defeaults to c("A", "H") indicating that the border wells should be removed for a 96 well plate. 
}
  \item{remove.cols}{
Vector of column indices that should be discarded according to the plate. Defeaults to c("1", "12") indicating that the border wells should be removed for a 96 well plate. 
}
  \item{progressbar}{
The type of progress bar used to show how far along the function is. Can be either \code{"window"}, \code{text} or \code{none}.
}
  \item{verbose}{
Should the function indicate what experiment it is currently reading. Defaults to \code{FALSE}.
}
  \item{save}{
Should the data be saved. Defeaults to TRUE.
}
  \item{shiny.input}{
Used for the shiny server.
}
  \item{session}{
Used for the shiny server.
}
}
\details{
This function is used for reading the protocols defining the setup for each plate and the associated dbf files storing the absorbance measures into R. Because of the nature of the function no examples of function usage is presented here. Instead you should look at the Sweave document accompanying the R package.
}
\value{
The ouput of the function is an A.data object of class \code{Absorbance}. This is a list with the following component  
 \item{meta.list}{This is a list of meta data objects.}
 \item{call}{A list containing information regarding the call to the function.}
 \item{auxiliary}{List of auxiliary data used by other functions.}
 \item{data}{List of data frames. The Absorbance data are stored in \code{raw.data} and the protocols are stored in \code{protocols}}
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
\seealso{
\code{\link{drugColorCorrection}, \link{bgModel}}
}
\examples{
##  Since this function is for reading data into R there are no examples of the usage
## See the associated Sweave document for detailed usage.
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{absorbance}

