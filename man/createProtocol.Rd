\name{createProtocol}
\alias{createProtocol}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Creates a protocol
}
\description{
This function is used for creating protocols in the shiny webserver. 
}
\usage{
createProtocol(n.doses = 16, n.backgrounds = 6, n.controls = 6, 
               fold = 2, max = 60, remove.edges = TRUE, drug = "drug")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{n.doses}{
The number of differnet doses.
}
  \item{n.backgrounds}{
Number of background probes.
}
  \item{n.controls}{
Number of controls
}
  \item{fold}{
The dilution fold,
}
  \item{max}{
The maximum dilution
}
  \item{remove.edges}{
Should the edges be removed.
}
  \item{drug}{
Name of the drug.
}
}
\note{
This function is used within the shiny server.
}
\seealso{
 \code{\link{editProtocol}},\code{\link{editDBFProtocol}},
 \code{\link{readProtocol}},\code{\link{useExistingProtocol}},
 \code{\link{useNewProtocol}},\code{\link{useTempProtocol}}
}
\examples{
## An example of a protocol with 16 doses created by a two fold dilution.
createProtocol(n.doses = 16, n.backgrounds = 6, n.controls = 6, fold = 2, max = 60, remove.edges = TRUE, drug = "Doxorubicin")
}
