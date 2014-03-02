\name{useExistingProtocol}
\alias{useExistingProtocol}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Use an existing protocol
}
\description{
Use an already established protocol with a new .dbf file.
}
\usage{
useExistingProtocol(A.data, protocolname, dbf.files)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{A.data}{
An \code{A.data} object
}
  \item{protocolname}{
Name of the protocol
}
  \item{dbf.files}{
Names of the dbf files the protocol should be assigned to
}
}
\note{
This function is used within the shiny server.
}
\seealso{
 \code{\link{createProtocol}},\code{\link{editProtocol}},
 \code{\link{editDBFProtocol}},\code{\link{readProtocol}},
 \code{\link{useTempProtocol}},\code{\link{useNewProtocol}}
}