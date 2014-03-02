\name{useTempProtocol}
\alias{useTempProtocol}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Use a protocol for specified dbf.files. 
}
\description{
Assign a protocol to the specified dbf.files.
}
\usage{
useTempProtocol(A.data, protocol = NULL, protocolname, dbf.files)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{A.data}{
An \code{A.data} object.
}
  \item{protocol}{
A protocol as that created by \code{createProtocol}. If left as NULL the A.data object is searched for a protocol with the specified protocolname.
}
  \item{protocolname}{
The protocol name to be used.
}
  \item{dbf.files}{
The .dbf files that are to use the protol
}
}
\note{
This function is used within the shiny server.
}
\seealso{
 \code{\link{createProtocol}},\code{\link{editProtocol}},
 \code{\link{editDBFProtocol}},\code{\link{readProtocol}},
 \code{\link{useExistingProtocol}},\code{\link{useNewProtocol}}
}