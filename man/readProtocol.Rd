\name{readProtocol}
\alias{readProtocol}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Reads a protocol into R.
}
\description{
This function reads the specified protocol into R in an appropriate format.
}
\usage{
readProtocol(xls.file, mistakeval = "X")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{xls.file}{
file path to the Excel file.
}
  \item{mistakeval}{
The value used in the setup to specify an error.
}
}
\note{
This function is used within the shiny server.
}
\seealso{
 \code{\link{createProtocol}},\code{\link{editProtocol}},
 \code{\link{editDBFProtocol}},\code{\link{useExistingProtocol}},
 \code{\link{useNewProtocol}},\code{\link{useTempProtocol}}
}