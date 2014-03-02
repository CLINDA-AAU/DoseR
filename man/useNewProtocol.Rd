\name{useNewProtocol}
\alias{useNewProtocol}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Use newly established protocol
}
\description{
This function assigns a newly established protocol to the specified dbf.files.
}
\usage{
useNewProtocol(A.data, protocol = NULL, protocolname, dbf.files)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{A.data}{
An A.data object.
}
  \item{protocol}{
The protocol that is to be used with the specified dbf file
}
  \item{protocolname}{
The name of the protocol
}
  \item{dbf.files}{
Character vector of the .dbf files that should use the protocol.
}
}
\note{
This function is used within the shiny server.
}
\seealso{
 \code{\link{createProtocol}},\code{\link{editProtocol}},
 \code{\link{editDBFProtocol}},\code{\link{readProtocol}},
 \code{\link{useExistingProtocol}},\code{\link{useTempProtocol}}
}