print.createMetaData <- function(x,...){
  
  drugvar <- x$auxiliary$passed.var$drugvar
  namevar <- x$auxiliary$passed.var$namevar
  correctionname <- x$auxiliary$passed.var$correctionname
  
  
  meta <- x$meta.list$metadata.full
  
  cat("Dose response experiment conducted on the drugs:\n")
  cat(unique(meta[, drugvar]), "\n")
  
  for(drug in unique(meta[, drugvar])){
    cat("\n")
    cat(length(unique(meta[meta[, drugvar] == drug,namevar] %w/o% correctionname)), 
        "cell lines have been exposed to", paste(drug, ":", sep = ""),  "\n") 
    y <- unique(meta[meta[, drugvar] == drug,namevar] %w/o% correctionname)
    y <- y[order(y)]
    if(length(y) > 6){
      cat(y[1:3], "...", y[(length(y)-3):length(y)], "\n")
    }else{
      cat(y, "\n")
    }   
  }
  cat("\n")
  
  y <- unique(meta[meta[,namevar] == correctionname, drugvar])
  
  if(length(y)==1){
    cat("The following drug is associated with a drug color correction plate\n")
    cat(y, "\n")
  }
  if(length(y)>1){
    cat("The following drugs are associated with drug color correction plates\n")
    cat(y, "\n")
  }
  if(length(y)==0){
    cat("No drugs are associated with drug color correction plates\n")
  }
  cat("\n")
  cat("Continue to function readDBFData. \n")
}

print.Absorbance <- function(x,...){
  drugvar <- x$auxiliary$passed.var$drugvar
  namevar <- x$auxiliary$passed.var$namevar
  correctionname <- x$auxiliary$passed.var$correctionname
  
  
  
  meta <- x$meta.list$metadata.full
  
  cat("Dose response experiment conducted on the drugs:\n")
  cat(unique(meta[, drugvar]), "\n")
  
  for(drug in unique(meta[, drugvar])){
    cat("\n")
    cat(length(unique(meta[meta[, drugvar] == drug,namevar] %w/o% correctionname)), 
        "cell lines have been exposed to", paste(drug, ":", sep = ""),  "\n") 
    y <- unique(meta[meta[, drugvar] == drug,namevar] %w/o% correctionname)
    y <- y[order(y)]
    if(length(y) > 6){
      cat(y[1:3], "...", y[(length(y)-3):length(y)], "\n")
    }else{
      cat(y, "\n")
    }   
  }
  cat("\n")
  
  cat("The molar masses of the drugs are:\n")
  print(x$auxiliary$mol.data)
  cat("If these are not correct you can change them using the function changeMolarMass")
  cat("\n\n")
  y <- unique(meta[meta[,namevar] == correctionname, drugvar])
  
  if(length(y)==1){
    cat("The following drug is associated with a drug color correction plate\n")
    cat(y, "\n\n")
    cat("Continue to function drugColorCorrection \n")
  }
  if(length(y)>1){
    cat("The following drugs are associated with drug color correction plates\n")
    cat(y, "\n\n")
    cat("Continue to function drugColorCorrection \n")
  }
  if(length(y)==0){
    cat("No drugs are associated with drug color correction plates\n\n")
    cat("Continue to function bgModel to perform normalisation \n")
  }
}


print.drugColorCorrection <- function(x,...){
  drugvar <- x$auxiliary$passed.var$drugvar
  namevar <- x$auxiliary$passed.var$namevar
  correctionname <- x$auxiliary$passed.var$correctionname
  
  
  
  meta <- x$meta.list$metadata.full
  
  cat("Dose response experiment conducted on the drugs:\n")
  cat(unique(meta[, drugvar]), "\n")
  
  for(drug in unique(meta[, drugvar])){
    cat("\n")
    cat(length(unique(meta[meta[, drugvar] == drug,namevar] %w/o% correctionname)), 
        "cell lines have been exposed to", paste(drug, ":", sep = ""),  "\n") 
    y <- unique(meta[meta[, drugvar] == drug,namevar] %w/o% correctionname)
    y <- y[order(y)]
    if(length(y) > 6){
      cat(y[1:3], "...", y[(length(y)-3):length(y)], "\n")
    }else{
      cat(y, "\n")
    }   
  }
  cat("\n")
  
  cat("Continue to function bgModel to perform normalisation \n")
}



print.bgModel <- function(x,...){
  drugvar <- x$auxiliary$passed.var$drugvar
  namevar <- x$auxiliary$passed.var$namevar
  correctionname <- x$auxiliary$passed.var$correctionname
  
  
  
  meta <- x$meta.list$metadata.full
  
  cat("Dose response experiment conducted on the drugs:\n")
  cat(unique(meta[, drugvar]), "\n")
  
  for(drug in unique(meta[, drugvar])){
    cat("\n")
    cat(length(unique(meta[meta[, drugvar] == drug,namevar] %w/o% correctionname)), 
        "cell lines have been exposed to", paste(drug, ":", sep = ""),  "\n") 
    y <- unique(meta[meta[, drugvar] == drug,namevar] %w/o% correctionname)
    y <- y[order(y)]
    if(length(y) > 6){
      cat(y[1:3], "...", y[(length(y)-3):length(y)], "\n")
    }else{
      cat(y, "\n")
    }   
  }
  cat("\n")
  
  cat("Continue to function bootstrap to simulate datasets used for confidence intervals \n")
}


print.bootstrap <- function(x,...){
  drugvar <- x$auxiliary$passed.var$drugvar
  namevar <- x$auxiliary$passed.var$namevar
  correctionname <- x$auxiliary$passed.var$correctionname
  
  
  
  meta <- x$meta.list$metadata.full
  
  cat("Dose response experiment conducted on the drugs:\n")
  cat(unique(meta[, drugvar]), "\n")
  
  for(drug in unique(meta[, drugvar])){
    cat("\n")
    cat(length(unique(meta[meta[, drugvar] == drug,namevar] %w/o% correctionname)), 
        "cell lines have been exposed to", paste(drug, ":", sep = ""),  "\n") 
    y <- unique(meta[meta[, drugvar] == drug,namevar] %w/o% correctionname)
    y <- y[order(y)]
    if(length(y) > 6){
      cat(y[1:3], "...", y[(length(y)-3):length(y)], "\n")
    }else{
      cat(y, "\n")
    }   
  }
  cat("\n")
  
  cat("Continue to function doseResponseModel to fit the dose response models. \n")
}



print.doseResponseModel <- function(x,...){
  drugvar <- x$auxiliary$passed.var$drugvar
  namevar <- x$auxiliary$passed.var$namevar
  correctionname <- x$auxiliary$passed.var$correctionname
  
  meta <- x$meta.list$metadata.full
  
  cat("Dose response experiment conducted on the drugs:\n")
  cat(unique(meta[, drugvar]), "\n")
  
  for(drug in unique(meta[, drugvar])){
    cat("\n")
    cat(length(unique(meta[meta[, drugvar] == drug,namevar] %w/o% correctionname)), 
        "cell lines have been exposed to", paste(drug, ":", sep = ""),  "\n") 
    y <- unique(meta[meta[, drugvar] == drug,namevar] %w/o% correctionname)
    y <- y[order(y)]
    if(length(y) > 6){
      cat(y[1:3], "...", y[(length(y)-3):length(y)], "\n")
    }else{
      cat(y, "\n")
    }   
  }
  cat("\n")
  
  #cat("Continue to function doseResponseModel to simulate datasets used for confidence intervals \n")
}


print.summary.DRdata <- function(x,...){
  drugvar <- x$auxiliary$passed.var$drugvar
  namevar <- x$auxiliary$passed.var$namevar
  correctionname <- x$auxiliary$passed.var$correctionname
  
  
  
  meta <- x$meta.list$metadata.full
  
  cat("Dose response experiment conducted on the drugs:\n")
  cat(unique(meta[, drugvar]), "\n")
  
  for(drug in unique(meta[, drugvar])){
    cat("\n")
    cat(length(unique(meta[meta[, drugvar] == drug,namevar] %w/o% correctionname)), 
        "cell lines have been exposed to", paste(drug, ":", sep = ""),  "\n") 
    y <- unique(meta[meta[, drugvar] == drug,namevar] %w/o% correctionname)
    y <- y[order(y)]
    if(length(y) > 6){
      cat(y[1:3], "...", y[(length(y)-3):length(y)], "\n")
    }else{
      cat(y, "\n")
    }   
  }
  cat("\n")
  
  #cat("Continue to function doseResponseModel to simulate datasets used for confidence intervals \n")
}

