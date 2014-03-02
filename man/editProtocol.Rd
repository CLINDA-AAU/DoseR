\name{editProtocol}
\alias{editProtocol}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Saves a protocol as an Excel file.
}
\description{
Saves the protocol as a temporary Excel file and returns the path to the file.
}
\usage{
editProtocol(protocol)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{protocol}{
A protocol created by  \code{createProtocol}.
}
}
\note{
This function is used within the shiny server.
}
\seealso{
 \code{\link{createProtocol}},\code{\link{editDBFProtocol}},
 \code{\link{readProtocol}},\code{\link{useExistingProtocol}},
 \code{\link{useNewProtocol}},\code{\link{useTempProtocol}}
}