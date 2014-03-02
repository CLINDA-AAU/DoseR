\name{createMetaData}
\alias{createMetaData}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Create an A.data object for storing the absorbance data.
}
\description{
This function creates an A.data object. This object stores all data related the dose response experiments which is saved in an .RData file. This ensures that the analyses is only conducted on new data including the time consuming bootstrap analysis. The most convinient approach for using the function is to create a metedata object directly from the .dbf filenames. Other possibilities include: 1) having an Excel document with experiment information, 2) A data.frame with experiment information, and 3) a list with dbf filenames and protocol filenames.  
}
\usage{
createMetaData(data = NULL, data.file = file.path(getwd(), "Absorbance"), 
               save = TRUE, 
               namevar = "Cellline", drugvar = "chemo", protocolvar = "R.protocol", 
               identifier = "identifier", timevar = "Hour", correctionname = "Control", 
               incubationvar = "incubation", doublingvar = NULL, format = "long", 
               dbf.path = getwd(), protocol.path = getwd(), dbf.files = NULL, 
               file.extension = ".dbf", protocol.files = NULL,  are.paths.full = TRUE,
               colnames = c("namevar", "drugvar", "protocolvar", 
                            "identifier", "timevar"), 
               sep = c(";"), namesep = " ", identifiersep = "_", 
               date.format = NULL, unit = "ug/ml", additional.metadata = NULL, 
               show.warnings = TRUE, update = TRUE, idcomb = NULL, idvar = "sampleid", 
               namecols = NULL, dbf.file.var = "dbf.file.name", shiny.input = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
This may be a data frame containing the metadata for the dose response experiments or the path to an Excel file containg the metadat. Not used when the metadata is inferred from the filename.
}
  \item{data.file}{
Character giving the path to the Excel metadata sheet created for the project. This is only created when the metadata is inferred from the filename.
}
  \item{save}{
Should the data be saved. 
}
  \item{namevar}{
The column name containg the name of the cell line. Defeaults to \code{Cellline}. 
}
  \item{drugvar}{
Character denoting the name of the column specifying the drug used in the experiment. Defeaults to \code{chemo}. 
}
  \item{protocolvar}{
Character denoting the name of the column specifying the protocol used in the experiment. Defeaults to \code{R.protocol}.
}
  \item{identifier}{
  Character denoting the name of the column specifying an id for the experiment. Defeaults to \code{identifier} and in combination with \code{namevar}, \code{drugvar}, and \code{timevar} it is expected to be a unique identifier. It is used to indicate what experiments are run together, so it should be identical for separate time points. The date where the drug is added may be a good value to use. 
}
  \item{timevar}{
The column name specifying the number of hours the cell line was exposed to the drug. Defeaults to \code{Hour}.
}
  \item{correctionname}{
The name used to indicate that this plate is to be used for correction of drugcolor. The name is to be included within the column \code{namevar} containg the names of the cell lines'. Defeaults to \code{Control}.
}
  \item{incubationvar}{
The column name specifying the number of hours the cell line incubated with the MTS assay. Defeaults to \code{NULL} for which 2 hours is assummed for all experiments.
}
  \item{doublingvar}{
Name of the column containg information of cell line doubling time when not treated with drug. This is only used if the experiment only contains 1 time point.
}
  \item{format}{
Format of the metadata. Can be either long or wide format the created Excel sheet will also be in this format. Defaults to \code{long}. 
}
  \item{dbf.path}{
Path to the directory storing the dbf files.
}
  \item{protocol.path}{
Path to the directory storing the protocol files.
}
  \item{dbf.files}{
The path to each individual .dbf file.
}
  \item{file.extension}{
The extension for the dBase files. defaults to \code{.dbf}
}
  \item{protocol.files}{
The path to each individual protocol file.
}    
 \item{are.paths.full}{
Logical. if \code{TRUE} the paths are assummed to full.
}
  \item{colnames}{
The order of the colnames specified by the filename for the .dbf file. Defeaults to \code{c("namevar", "drugvar", "protocolvar","identifier", "timevar")} which means that part 1) includes the name of the cell line, 2) specifies what drug was used in the experiment, 3) specifies the name of the protocol, 4) specify the unique identifier, and 5) the number of hours the cell line were exposed to a drug. Thus if the the cell line OCI-Ly7 was treated with Doxorubicin according to protocol DOX18 on the November 30th 2013, and exposed to the drug for 48 hours, the file name would be: OCI-Ly7;Doxorubicin;DOX18;20131130;48. The separator ";" is explained below.
}
  \item{sep}{
The character used to specify column separation in the filename. Defealts to ";" rsulting in a filename as described above.
}
  \item{namesep}{
Aditional separator for cell line name, defeaults to "_". The namecharater in the file name can further be used to split up experiments. E.g. the filenames 
OCI-Ly7_cond1;Doxorubicin;DOX18;20131130;48 and OCI-Ly7_cond2;Doxorubicin;DOX18;20131130;48 says that cell line OCI-Ly7 is treated with Doxorubicin under two different conditions.
}
  \item{identifiersep}{
Aditional separator for identifier, defeaults to "_". The indentifier can be the date used for the setup, however if the same cell line is treated with the same drug on the same date the filename will not be unique. Then we may use an additional seperator in the identifier name. E.g.\ 
OCI-Ly7;Doxorubicin;DOX18;20131130_1;48 and OCI-Ly7;Doxorubicin;DOX18;20131130_2;48 indicate that the two different experiments.
}
  \item{date.format}{
If the identifier indicates the setup date for the experiment the format of the date is specified.
Defeaults to \code{NULL}.
}
  \item{unit}{
The unit for the drug concentration. Defeaults to \code{ug/ml}. where u is used as substitute for \eqn{\mu}. If a column of the data supplied to argument \code{data} is named \code{unit} this column will be used as the unit.
}
  \item{additional.metadata}{
An optinal additional metadata set containg information regarding the cell lines. The first column must specify the name of the cell line.
}
  \item{show.warnings}{
Should warnings be displayed.
}
  \item{update}{
Do you want to update an already made metasheet or create a new.
}
  \item{idcomb}{
Combination of column names that identifies an experiment id (Optional)
}
  \item{idvar}{
The variable name used for id (Optional)
}
  \item{namecols}{
If the name part of the file name is separated using \code{namesep} the resulting additional column names are given here.
}
  \item{dbf.file.var}{
Name of the column in the created metadata containg the filename (optional)
}
  \item{shiny.input}{
Used for the shiny web application.
}
}
\details{
When conducting dose response experiments the amount of .dbf becomes large thus metadata is normally stored in Excel sheets. This function converts such data into an \code{A.data} object of class \code{createMetaData}. Since many mistakes occur when copy pasting filenames into Excel the function can also be used to create metadata from filenames provided the filenames are of a certain structure, such as: OCI-Ly7;Doxorubicin;DOX18;20131130;48.dbf.
}
\value{
The ouput of the function is an A.data object of class \code{createMetaData}. This is a list with the following component  
 \item{meta.list}{This is a list of meta data objects.}
 \item{call}{A list containing information regarding the call to the function.}
 \item{auxiliary}{List of auxiliary data used by other functions.}

The meta.list output contains
\item{metadata.full}{Which is the metadata for all dose response experiment in eihter long or wid format.}
\item{metadata.correction}{Which is the metadata for dose response experiment used for correction of drug colour.}
\item{metadata.new}{Metadata for all experiments not used to correct for drug colour.}
\item{additional}{The metadata sheet supplied to \code{additional.metadata}.}
\item{metadata.long}{Metadata for all dose response experiments sorted in long format.}
}
\references{
Steffen Falgreen et al. Exposure time independent summary statistics for assessment of drug dependent cell line growth inhibition (2013)
}
\author{
The function was written at department of haematology, Aalborg University Hospital and maintained by Steffen Falgreen.
}
\note{
Following the creation of the A.data object the function \code{\link{readDBFData}} is used to read the dBase files into R.
}
\seealso{
\code{\link{readDBFData}}
}
\examples{
## Example of creating the metadata based on file names

## A.data <- createMetaData(
##  dbf.path       = data.path,     # data.path is a predifined directory storing the dbf files
##  protocol.path  = protocol.path, # protocol.path is a predifined directory storing the Excel protocol files
##  colnames       = c("namevar", "drugvar", "protocolvar",   ## The order of colnames specified by the file names:
##                     "identifier", "timevar"),              ## OCI-Ly7;Doxorubicin;DOX18;20131130;48
##  sep            = ";",          # The separator used in the filename as noted above
##  namevar        = "name",     
##  drugvar        = "chemo",
##  protocolvar    = "R.protocol",
##  identifier     = "identifier",
##  timevar        = "Hour",
##  namecols       = "serum",     
##  date.format    = "%d%m%y",
##  correctionname = "Control",    # The name specified in the column "namevar" for drug colour correction plates  
##  unit           = "ug/ml",
##  additional.metadata = file.path(BCell.ext.dir, "cell.line.metadata.xls"),
##  show.warnings  = TRUE,
##  update         = TRUE,
##  format         = "long",
##  data.file      = file.path(BCell.gen.dir, "Absorbance.test"))
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
