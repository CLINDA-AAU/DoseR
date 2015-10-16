DoseR: Dose response analysis in R using the time independent G-model. 
=================================

The R-package `DoseR` is a set of tools used for dose response analysis.

The packages is features (or, is planned to feature):
* Pre-processing of dose response data
* Estimation of the G-model


Installation
------------
If you wish to install the latest version of `DoseR` directly from the master branch here at GitHub, run 

```R
# Install necessary packages 

# First from bioconductor
source("http://bioconductor.org/biocLite.R")
biocLite("prada")

# Then from CRAN
install.packages(c("foreign", "gdata", "plyr", "WriteXLS", "nlme", "RJSONIO",
                   "RCurl", "XML", "RColorBrewer", "Hmisc", "animation"))
                   
#Then install the rvest package version 0.2.0 from CRAN
package <- "https://cran.r-project.org/src/contrib/Archive/rvest/rvest_0.2.0.tar.gz"
install.packages(package, repos = NULL, type = "source")

# From GitHub 
install.packages("devtools")
devtools::install_github("rstudio/shiny-incubator")
devtools::install_github("oncoclass/DoseR", dependencies = TRUE)
```

`DoseR` is still under development and should be considered unstable. Be sure that you have the [package development prerequisites](http://www.rstudio.com/ide/docs/packages/prerequisites) if you wish to install the package from the source.

**Note:** The interface and function names may still see significant changes and
modifications!


Using DoseR
----------
A small tutorial goes here!

### Conventions and interface
The supplied expression matrices to be converted are arrange in the tall 
configuration. I.e. we have features in the rows and samples in columns.


References
----------
References goes here!
