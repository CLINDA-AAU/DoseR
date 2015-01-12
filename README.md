DoseR: Dose response analysis in R using the time independent G-model. 
=================================

The R-package `DoseR` is a set of tools used for dose response analysis.

The packages is features (or, is planned to feature):
* Pre-processing of dose response dataœ“
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
install.packages("foreign")
install.packages("gdata")
install.packages("plyr")
install.packages("WriteXLS")
install.packages("nlme")
install.packages("RJSONIO")
install.packages("RCurl")
install.packages("XML")
install.packages("RColorBrewer")
install.packages("shiny")
install.packages("Hmisc")
install.packages("animation")

# From GitHub 
install.packages("devtools")  # Uncomment if devtools is not installed

devtools::install_github("rstudio/shiny-incubator")

devtools::install_github("Falgreen/DoseR", 
                         auth_token = "...", 
                         dependencies = TRUE)
```

`DoseR` is still under heavy development and should be considered unstable. Be sure that you have the [package development prerequisites](http://www.rstudio.com/ide/docs/packages/prerequisites) if you wish to install the package from the source.

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
