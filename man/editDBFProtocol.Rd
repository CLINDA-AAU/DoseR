\name{editDBFProtocol}
\alias{editDBFProtocol}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Edit protocol associated with an .dbf file
}
\description{
This function edits a protocol associated with a certain dbf file. The protocol element of the dbf file will be renamed in accordance with the newly established name. The protocol opens in Excel with a newly assigned name. The new name must not be replaced.
}
\usage{
editDBFProtocol(A.data, dbf.file = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{A.data}{
An \code{A.data} created by \code{createMetadata}.
}
  \item{dbf.file}{
The name of the .dbf file for which the protocol should be altered.
}
}
\note{
This function is used within the shiny server.
}
\seealso{
 \code{\link{createProtocol}},\code{\link{editProtocol}},
\code{\link{readProtocol}},\code{\link{useExistingProtocol}},
 \code{\link{useNewProtocol}},\code{\link{useTempProtocol}}
}
