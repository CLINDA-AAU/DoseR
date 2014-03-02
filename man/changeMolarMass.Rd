\name{changeMolarMass}
\alias{changeMolarMass}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Change or set the molar mass of a drug.
}
\description{
The function \code{readDBFData} automatically downlads information regarding molar mass of the drugs. The information is e.g. used for changing unit in a plot.
}
\usage{
changeMolarMass(A.data, drug, mol.mass, save = TRUE, shiny.input = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{A.data}{
An A.data object created by the function \code{readDBFData}
}
  \item{drug}{
Character indicating the drug for which the molar mass is to be changed.
}
  \item{mol.mass}{
Numeric value specifying the molar mass of the drug
}
  \item{save}{
Should the data be saved for later use. Defeaults to TRUE.
}
  \item{shiny.input}{
Used for the shiny webtool
}
}
\details{
When it is not possible to find information of mol.mass it is necessary to add it yourself.
}
\value{
An A.data object with the molar mass changed for specified drug.
}
\references{
Steffen Falgreen et al. Exposure time independent summary statistics for assessment of drug dependent cell line growth inhibition (2013)
}
\author{
The function was written at department of haematology, Aalborg University Hospital and maintained by Steffen Falgreen.
}
\examples{
require(DoseR)
## load Dose Response data
data(A.data)

## print of the A.data object
print.Absorbance(A.data)

## Say the molar mass of Rituximab was wrong and you wanted to set it 
## to another value e.g. 100000.

A.data <- 
   changeMolarMass(A.data, drug = "Rituximab", 
   mol.mass = 100000, save = FALSE) 
   
# Setting save to TRUE ensures that the molar mass is not overwritten 
# when using the function readDBFdata with update set to TRUE.

}

