
##
##
## Function for printing the analysis flow
##
##


jsonToList <- function(json) {
  jsonAsList <- fromJSON(json, simplify = StrictCharacter, nullValue = "")
  return(jsonAsList)
}

jsonToDataFrame <- function(json) {
  jsonAsList <- jsonToList(json)
  jsonAsDataFrame <- ldply(jsonAsList,function(x) as.data.frame(t(unlist(as.numeric(x)))))
  return(jsonAsDataFrame)
}


printAnalysisFlow <- function(A.data){
  if(!"call" %in% names(A.data))
    return()
  funs <- names( A.data$call) %w/o% 
    c("record", "createMetaDataDir", "createMetaDataXls", 
      "createMetaDataList","createMetaDataListDataFrame",
      "growthModel", "createGIData", "combineGIdata", "isoreg.DRdata", "summary.DRdata")
  
  
  funs <- c("createMetaData", "readDBFData", "drugColorCorrection", "bgModel", 
    "bootstrap", "doseResponseModel")
  funs <- funs[funs %in% names( A.data$call)]
  
  for(names in funs){
    cat("A.data <- \n     ")
    if(names == "createMetaData"){
      A.data$auxiliary$shiny.calls$createMetaData
   # print(A.data$auxiliary$shiny.calls[[names]])
    }else{
      cat(A.data$call[[names]][["shiny.call"]])
    }
    
    cat("\n\n")
  }
}
# choose.dir <- function(old.dir = getwd()){
#   dbf.file <- try(file.choose(), silent=TRUE)
#   if(is(dbf.file, 'try-error')){
#     old.dir
#   }else{
#     paste(strsplit(dbf.file, "/")[[1]][1:(length(strsplit(dbf.file, "/")[[1]])-1)], collapse = "/")
#   }
# }


choose.dir <- function(old.dir = getwd(), caption = "Select directory"){
  dbf.file <- tk_choose.dir(default = old.dir, caption = caption)
  if(is.na(dbf.file)){
    old.dir
  }else{
    dbf.file
  }
}


# my.file.choose <- function(old.file = getwd()){
#   file <- try(file.choose(), silent=TRUE)
#   if(is(file, 'try-error')){
#     old.file
#   }else{
#     file
#   }
# }

my.file.choose <- function(old.file = getwd()){
  file  <- tk_choose.files()#try(file.choose(), silent=TRUE)
  
  #if(is(analysisfile, 'try-error')){
  if(length(file) == 0 ){
    old.file
  }else{
    file
  }
}

"%w/o%" <- function(x, value) x[!x %in% value]

##
##
## Function for numbering the 
## different calls, not yet implemented
##
##

callNumbering <- function(fun = NULL, record = NULL){
  if(is.null(record)){
    functions  = c("createMetaData","createMetaDataDir", 
                   "createMetaDataXls", "createMetaDataDataFrame",
                   "createMetaDataList",
                   "readDBFData", "drugColorCorrection",
                   "bgModel", "bootstrap", "growthModel", "createGIData", 
                   "combineGIdata", "summary.DRdata")
    record <- data.frame(function.name  = functions,
                         numbering = c(rep(1, 4), 1:(length(functions)-4)),
                         last.visited = as.POSIXct("2011-04-04 14:18:58", tz="GB"))
    row.names(record) <- record$function.name
  }
  
  record[fun, "last.visited"] <- Sys.time()
  record[,"function.name" ] <- row.names(record)
  return(record)
}


##
##
## Function for checking dependency of the 
## different calls, not yet implemented
##
##

dependencyCheck <- function(A.data, call = "summary.DRdata"){
  not.visited = as.POSIXct("2011-04-04 14:18:58", tz="GB")
  dependency <- list("createMetaData" = NULL, 
                     "createMetaDataDir" = NULL, 
                     "createMetaDataXls" = NULL, 
                     "createMetaDataDataFrame" = NULL,
                     "createMetaDataList" = NULL,
                     "readDBFData" = c("createMetaData"), 
                     "drugColorCorrection" = c("createMetaData", "readDBFData"),
                     "bgModel" =  c("createMetaData", "readDBFData", "drugColorCorrection"), 
                     "bootstrap" =  c("createMetaData", "readDBFData", 
                                      "drugColorCorrection", "bgModel"), 
                     "growthModel" = c("createMetaData", "readDBFData", 
                                       "drugColorCorrection", "bgModel",  
                                       "bootstrap" ), 
                     "createGIData" = c("createMetaData", "readDBFData", 
                                        "drugColorCorrection", "bgModel",  
                                        "bootstrap" ), 
                     "combineGIdata" = c("createMetaData", "readDBFData", 
                                         "drugColorCorrection", "bgModel",  
                                         "bootstrap" ), 
                     "summary.DRdata" = c("createMetaData", "readDBFData", 
                                          "drugColorCorrection", "bgModel",  
                                          "bootstrap" ))
  
  fun.vec <- dependency[[call]]
  call.record <- A.data$call$record
  missing.fun.calls <- fun.vec %w/o% names(A.data$call)
  call.record[fun.vec, ]
  call.flow <- ""
  
  record <- A.data$auxiliary$record
  
  
  
  fun.vec2 <- fun.vec %w/o% "createMetaData"
  record[,fun.vec2]
  colnames(record)
}

##
##
## Help functions for handling extensions
##
##

fileExt <- function(x) {
  db <- grepl("\\.[^.]+\\.(gz|bz2|xz)$", x)
  ans <- sub(".*\\.", "", x)
  ans[db] <- sub(".*\\.([^.]+\\.)(gz|bz2|xz)$", "\\1\\2", 
                 x[db])
  ans
}

no.extension <- function(astring) { 
  if(fileExt(astring) == astring){
    return(astring)
  }else{
    
    if (substr(astring, nchar(astring), nchar(astring))==".") { 
      return(substr(astring, 1, nchar(astring)-1)) 
    } else { 
      no.extension(substr(astring, 1, nchar(astring)-1)) 
    }
  } 
} 

no.extension.vec <- function(x){
  vec <- vector()
  for(i in 1:length(x))
    vec[i] <- no.extension(x[i])
  return(vec)
} 


##
##
## Function that extracts information 
## regarding a drug from wikipedia
##
##
# drugInfo <- function(drug){
# 
#   the_url <- paste0("http://en.wikipedia.org/wiki/", drug)
#   table <- NULL
#   if(url.exists(the_url)) table <- readHTMLTable(doc = rvest::html(the_url))[[1]]
# 
#   if(!is.null(table)){
# 
#     #flyttet dette {
#     for(i in 1:ncol(table))
#       table[,i] <- as.character(table[,i])
# 
#     table <- table[!is.na(table[,2]), ]
#     table <- table[table[,1] != "", ]
#     rownames(table) <- table[,1]
# 
#     table[,2] <- gsub("????\u0080\u0093", "-", table[,2])
#     table[,2] <- gsub("????\u0084\u009e", "Rx", table[,2])
#     # } og hertil, ind i if sentenc
# 
#     chem.info <- list()
# 
#     if(any(grepl("IUPAC", table[,1]))){
#       wh<- which(grepl("IUPAC", table[,1])) +1
#       chem.info$"Systematic (IUPAC) name" <- table[wh,1]
#     }
# 
#     if(any(grepl("Target", table[,1])))
#       chem.info$"Target" <- table["Target", 2]
# 
#     chem.info$"Clinical data" <- table[c("Trade names", "AHFS/Drugs.com", "MedlinePlus",
#                                          "Pregnancy cat.", "Legal status", "Routes"), 2, drop = FALSE]
#     chem.info$"Pharmacokinetic data" <- table[c("Bioavailability", "Metabolism", "Half-life", "Excretion"), 2, drop = FALSE]
#     chem.info$"Identifiers" <- table[c("CAS number", "ATC code", "PubChem",
#                                        "DrugBank", "ChemSpider", "UNII", "KEGG",
#                                        "ChEBI", "ChEMBL") , 2, drop = FALSE]
# 
#     chem.info$"Identifiers"[,1] <- gsub("\\s[YN]", "", chem.info$"Identifiers"[,1], perl = TRUE)
# 
#     chem.info$"Chemical data" <- table[c("Formula", "Molar mass") , 2, drop = FALSE]
# 
#     chem.info$"SMILES" <-
#       gsub("\\s[YN]", "", gsub("SMILES\n", "", table[grepl("SMILES", table[,1]),1]))
# 
#     InChI <-
#       gsub("\\s[YN]", "", gsub("InChI\n\nInChI=", "", table[grepl("InChI", table[,1]),1]))
# 
#     if(length(InChI) > 0) {
#       chem.info$"InChI" <-  strsplit(InChI, "\nKey:")[[1]]
#       names(chem.info$"InChI") <- c("InChI", "Key")
#     }
#   }else{
#     chem.info <- "Drug not found"
#   }
# 
#   return(chem.info)
# }

drugInfo <- function(drug){
  data(drug.mass)
  chem.info <- list()
  
  if(drug %in% drug.mass$drug_name){
    table <- drug.mass[drug.mass$drug_name == drug, ]
    chem.info$"Chemical data" <- data.frame(table$molar_mass, row.names = "Molar mass")
    }else{
    chem.info <- "Drug not found"
  }
  
  return(chem.info)
}



progressBar <- function(title = "progress bar", min = 0,
                        max = 100, width = 300, window = TRUE){
  operating.system <- sessionInfo()[1]$R.version$os
  glb.window.progrs <- window 
  if(!glb.window.progrs){
    txtProgressBar(min = min, max = max, style = 3)
  } else{
    if(grepl("darwin", operating.system )){
      tkProgressBar(title = title, min = min,
                    max = max, width = width)
    }else{
      winProgressBar(title = title, min = min,
                     max = max, width = width)
    }
  }
}

setProgressBar <- function(pb, i,label=NULL,
                           window = TRUE){
  if(is.null(label))
    label <- paste( round(i*100, 0),
           "% done")
  operating.system <- sessionInfo()[1]$R.version$os
  
  if(!window){
    setTxtProgressBar(pb, i)
  }else{
    if(grepl("darwin", operating.system )){
      setTkProgressBar(pb, i, label=label)
    }else{
      setWinProgressBar(pb, i, label=label)
    }
  }
}


##
##
## Create a call to display in shiny and 
## save information for later us
##
##

createCall <- function(call = match.call(), fun){
  
  new.call <- as.list(call)
  
  myfor <- formals(eval(new.call[[1]]))
  
  newnames <- names(myfor) %w/o% names(new.call) 
  for(i in newnames)
    new.call[[i]] <- myfor[[i]]
  
  if("shiny.input" %in% names(new.call))
  new.call$shiny.input <- NULL
  if("session" %in% names(new.call))
  new.call$session <- NULL
  
  new.call[[1]] <- NULL
  
  
  vec <- vector()
  for(i in 1:length(new.call))
    vec[i] <- ifelse(is.character(new.call[[i]]) & 
                       length(new.call[[i]]) == 1, '"',  "")
  
  
  this.shiny.call <- 
    paste(eval(fun), '(', paste(names(new.call), ' = ', vec, 
                                new.call, vec, collapse = ', ', sep =""), ')', 
          collapse = '', sep = '')
  
  
  new.call$shiny.call <- this.shiny.call
  new.call$call <- call
  return(new.call)
}




##
##
## Create metadata from a directory of files
##
##

createMetaDataDir <- function(
  dbf.files        = NULL,
  file.extension = ".dbf",
  protocol.files   = NULL,
  dbf.path         = getwd(),
  protocol.path    = getwd(),
  sep = c(";"),
  colnames = c("namevar", "drugvar", "protocolvar", 
               "identifier", "timevar"),
  namevar        = "Cellline",
  drugvar        = "chemo",
  protocolvar    = "R.protocol",
  identifier     = "identifier",
  timevar        = "Hour",
  incubationvar  = "incubation",
  doublingvar    = "T0",
  namesep        = " ",
  identifiersep = "_",
  namecols = "serum",
  date.format = "%d%m%y",
  correctionname = "Control",
  unit = "ug/ml",
  additional.metadata= NULL,
  idcomb   = NULL,#c(namevar, drugvar, identifier) # vars used to generate id variable
  idvar    = "sampleid", 
  dbf.file.var = "dbf.file.name",
  show.warnings = TRUE,
  update = TRUE,
  format = "long",
  save   = TRUE,
  data.file = file.path(getwd(), "Absorbance"),
  shiny.input = NULL)
{
  
  #############################################
  ##
  ##     Function starts 
  ##
  #############################################
  call2 <- match.call()
  this.call <- list()
  ## formals with default arguments
  myfor <- formals(createMetaDataDir) 
  for ( v in names(myfor)){
    if (!(v %in% names(call2))){
      this.call[[v]] <- myfor[[v]]
    }else{
      this.call[[v]] <- call2[[v]]
    }  ## if arg is missing I add it
  }
  
  this.call$call <- call2
  
  namevar2 <- paste(namevar, 2, sep = "")
  
  diff <- setdiff(c("namevar", "drugvar", "protocolvar", "identifier", "timevar"),
                  colnames) 
  
  if(length(diff) > 0)
    stop("The colnames need to include the variables ", diff)
  # check colnames
  
  # nameCheck <- c(namevar, drugvar, protocolvar, identifier, timevar)
  #  names(nameCheck) <- c("namevar", "drugvar", "protocolvar", "identifier", "timevar")
  
  # mat <- as.data.frame(matrix( nrow = 5, ncol = 4))
  #  colnames(mat) <- c("Variable", "Variable_name", "colnames")
  # mat[,1] <- c("namevar", "drugvar", "protocolvar", "identifier", "timevar")
  #  mat[,2] <- c(namevar, drugvar, protocolvar, identifier, timevar)
  # mat[,3] <- colnames
  #  mat[,4] <- mat[,2] == mat[,3]
  #  if(any(!mat[,4]))
  #   stop("colnames and varnames are not the same", print(mat))
  
  if(class(additional.metadata) == "character"){
    additional.metadata <- 
      data.frame(read_excel(additional.metadata))
  }
  
  if(!is.null(additional.metadata)){
    if(!(namevar %in% colnames(additional.metadata)))
      stop('the column name "', namevar, '" is not included in the supplied metadata')
  }
  
  if(is.null(idcomb))
    if(format == "wide"){
      idcomb  <- c(namevar, drugvar, identifier)
    }else{
      idcomb  <- c(namevar, drugvar, identifier, timevar)
    }
  
  
  idcombwide  <- c(namevar, drugvar, identifier)
  
  idcomblong  <- c(namevar, drugvar, identifier, timevar)
  
  
  #############################################
  ##
  ##     Read old data 
  ##
  #############################################
  
  xls.file  <- paste(data.file, ".xlsx",   sep = "")
  data.file <- paste(data.file, ".RData", sep = "")
  
  
  
  prior <- all(file.exists(data.file), 
               file.exists(xls.file),
               update)
  
  if(prior){
    load(data.file)
    this.call <- this.call
    old.call  <- old$call$createMetaDataDir
    if(length(intersect(names(old.call) , names(this.call))) != 
         length(names(old.call)) |
         length(intersect(names(old.call) , 
                          names(this.call))) != length(names(this.call) ))
      prior <- FALSE
    vec <- ""
    if(prior){
      call.names <- names(old.call) %w/o% c("call", "format", "update", "shiny.input")
      
      vec <- list()
      for(call.iter in call.names){
        vec[[call.iter]] <- this.call[[call.iter]] == old.call[[call.iter]]
      }
      if(!any(unlist(vec)))
        prior <- FALSE
    }
    if(!prior){
      
      excel.2 <- paste(no.extension(xls.file), " last changed ",
                        file.info(xls.file)$mtime, 
                        ".xlsx", sep = "")
      
      file.copy(from=xls.file, to=excel.2)
      
      data.file2 <- paste( no.extension(xls.file), " last changed ",
                           file.info(xls.file)$mtime, 
                           ".RData", sep = "")
      file.copy(from=data.file, to=data.file2)
      
      warning("The input variables have changed but the data is saved to an allready exisiting project without setting update to FALSE.\n\n", 
              "The input update is changed to FALSE. \n\n", 
              "The old Excel metasheet is backed up to:\n",
              excel.2,
              " \n\nThe old R.data file is backed up to:\n",
              data.file2)
      print(vec)
    }
  }
  
  meta.sheets <- list()
  if(prior){
    #load(data.file)
    
    sheet.names <- excel_sheets(xls.file)
    sheet.count <- length(sheet.names)
    
    for(i in 1:sheet.count) {
      meta.sheets[[i]] <- data.frame(read_excel(xls.file, sheet = i))
    }
    names(meta.sheets) <- excel_sheets(xls.file)[1:sheet.count]
    metadata.xls <- meta.sheets[["Data"]]
    metadata.xls[,"unit"] <- gsub("????", "", metadata.xls[,"unit"])
    time.points <- colnames(metadata.xls)[
      grepl("date.", colnames(metadata.xls)) & colnames(metadata.xls) != "setupdate" |
        colnames(metadata.xls) == "date" ]
    
    
    for(i in time.points){ 
      metadata.xls[, i] <- 
        as.POSIXct(strptime(metadata.xls[, i], "%Y-%m-%d %H:%M:%S") )
    }
    if(!is.null(date.format))
      metadata.xls[, "setupdate"] <- as.Date(metadata.xls[, "setupdate"])
    as.POSIXct(as.Date(metadata.xls[, "setupdate"]))    
    
    
    if(old$call$createMetaData$format == "wide" &
         format == "long"  ){
      if(idvar == "id")
        metadata.xls$id2 <- metadata.xls[,idvar]
      if(timevar == "time")
        metadata.xls$time2 <- metadata.xls[,timevar]
      
      
      metadata.xls.long <- 
        reshape(metadata.xls, direction = "long", 
                varying = list(
                  protocolvar = colnames(metadata.xls)[
                    grepl(paste(protocolvar,".",sep = ""), colnames(metadata.xls))], 
                  timevar = colnames(metadata.xls)[
                    grepl(paste(timevar,".",sep = ""), colnames(metadata.xls))],
                  date = colnames(metadata.xls)[
                    grepl("date.", colnames(metadata.xls))],
                  file.name = colnames(metadata.xls)[
                    grepl(dbf.file.var, colnames(metadata.xls))],
                  omit = colnames(metadata.xls)[
                    grepl("omit.", colnames(metadata.xls))],
                  reason = colnames(metadata.xls)[
                    grepl("reason.", colnames(metadata.xls))],
                  incubation = colnames(metadata.xls)[
                    grepl(paste(incubationvar, ".", sep = ""), colnames(metadata.xls))],
                  idlong = colnames(metadata.xls)[
                    grepl("idlong.", colnames(metadata.xls))] 
                ),
                v.names = c(protocolvar, timevar,"date", 
                            dbf.file.var, "omit", "reason", incubationvar,
                            "idlong"))
      
      if(idvar == "id"){
        metadata.xls.long$id <- metadata.xls.long$id2
        metadata.xls.long <- metadata.xls.long[,colnames(metadata.xls.long) != "id2"]
      }else{
        metadata.xls.long <- metadata.xls.long[,colnames(metadata.xls.long) != "id"]
      }
      if(timevar == "time"){
        metadata.xls.long$time <- metadata.xls.long$time2
        metadata.xls.long <- metadata.xls.long[,colnames(metadata.xls.long) != "time2"]
      }else{
        metadata.xls.long <- metadata.xls.long[,colnames(metadata.xls.long) != "time"]
      }
      
      metadata.xls <- metadata.xls.long[metadata.xls.long$idlong != "", ]  
      
      metadata.xls[, idvar] <- apply(metadata.xls[, idcomblong], 1, 
                                     function(x) paste(x, collapse = " "))
    }
    
    
    if(old$call$createMetaData$format == "long" &
         format == "wide"  ){
      
      metadata.xls[, idvar] <- apply(metadata.xls[, idcomb], 1, 
                                     function(x) paste(x, collapse = " "))
      
      metadata.xls <- 
        reshape(metadata.xls, 
                timevar = timevar, 
                idvar = "idwide", 
                direction = "wide",
                v.names =  c(timevar,  protocolvar,
                             "date", dbf.file.var, "omit" ,"reason", 
                             "idlong", incubationvar))
      
    }
  }
  if(is.null(additional.metadata)){
    if(prior){
      additional.metadata  <- meta.sheets[["additional.metadata"]]
    }else{
      additional.metadata <- NULL
    }
  }  
  
  #############################################
  ##
  ##     create a dataset with the files 
  ##
  #############################################
  
  if(is.null(dbf.files))    
    dbf.files <- dir(dbf.path, pattern = file.extension, full.names=TRUE)
  
  
  dbf.files4 <- basename(dbf.files)
  dbf.path2 <- dirname(dbf.files)
  
  dbf.files2 <- strsplit(dbf.files4,"\\.")
  #ext <- sapply(dbf.files2, function(x) x[length(x)])  
  #ext2 <- paste("\\.", ext, sep = "")
  
  ext2 <- file.extension
  for(i in ext2)
    dbf.files3 <- gsub(i, "", dbf.files4)
  
  spl <- strsplit(dbf.files3, sep)
  
  discard <- lapply(spl, length) != length(colnames)
  
  if(any(discard)){
    if(show.warnings)
      if(sum(discard) > 1){
        warning("The following entries in the file directory are discarded due to incorrect filenames: ", paste(dbf.files[discard], sep = "", colapse = ", ") )
      }else{
        warning("The following entry in the file directory is discarded due to incorrect filename: ", dbf.files[discard] )
      }
    spl    <- spl[!discard]
    dbf.files  <- dbf.files[!discard]
    dbf.files2 <- dbf.files2[!discard]
    dbf.files3 <- dbf.files3[!discard] 
    dbf.files4 <- dbf.files4[!discard]
    dbf.path2  <- dbf.path2[!discard]
  } 
  
  
  dbf.info <- data.frame(file        = paste(dbf.files3, file.extension, sep = ""), 
                         file.path   = dbf.path2, 
                         file.full   = dbf.files,
                         last.change = file.info(dbf.files)$mtime, 
                         stringsAsFactors=FALSE)
  
  row.names(dbf.info) <- dbf.files3
  
  if(is.null(protocol.files))    
    protocol.files <- dir(protocol.path, pattern = ".xls", full.names=TRUE)
  
  protocol.info <- data.frame(
    file        = basename(protocol.files),
    file.path   = dirname(protocol.files) ,                          
    file.full   = protocol.files, 
    last.change = file.info(protocol.files)$mtime, 
    stringsAsFactors=FALSE)
  
  row.names(protocol.info) <- no.extension.vec(basename(protocol.files))
  
  if(!prior){
    file.info <- list(dbf = dbf.info , protocol = protocol.info)
  }else{
    file.info <- old$auxiliary$file.info
    file.info$dbf <- dbf.info
    file.info$protocol = protocol.info
  }
  basename(protocol.files)
  
  met <- t(sapply(spl, t))
  
  colnames(met) <- colnames
  
  met <- as.data.frame(met)
  
  if("incubationvar" %in% colnames(met)){
    colnames(met)[match(c("namevar", "drugvar", 
                          "identifier", "timevar", "protocolvar",
                          "incubationvar"),  colnames(met) )] <- 
      c(namevar, drugvar, identifier, timevar, protocolvar, incubationvar)
  }else{
    colnames(met)[match(c("namevar", "drugvar", "identifier",
                          "timevar", "protocolvar"),  colnames(met) )] <- 
      c(namevar, drugvar, identifier, timevar, protocolvar)
  }
  
  if(!is.null(date.format)){
    if(!is.null(identifiersep)){
      met[, "setupdate"] <- 
        as.Date(unlist(lapply(strsplit(as.character(met[, identifier]), identifiersep), 
                              function(x) x[1])), format = date.format)
      met[, "replicate"]<- 
        unlist(lapply(strsplit(as.character(met[, identifier]), identifiersep), function(x) x[2]))
    }else{
      met[, "setupdate"] <- as.Date(as.character(met[, identifier]), format = date.format)
    }
  }
  
  
  #met[, "id"] <- paste(met[, namevar], met[, drugvar], met[, identifier], sep = " ")
  #if(identifier != "id")
  #  met <- met[, colnames(met) %w/o% identifier]
  
  ## the identifier is changed to id as this is used in following functions
  #identifier <- "id"
  met[, dbf.file.var]  <- dbf.files3
  met[, "date"] <- (file.info(dbf.files)$mtime)
  if(! unit %in% colnames(met))
    met[, "unit"] <- unit
  met[, "omit"] <- FALSE
  met[, "reason"] <- ""
  
  
  
  
  if(! incubationvar %in% colnames(met))
    met[, incubationvar] <- 2
  
  
  
  met[, idvar] <- apply(met[, idcomb], 1, function(x) paste(x, collapse = " "))
  
  met[, "idwide"] <- apply(met[, idcombwide], 1, function(x) paste(x, collapse = " "))
  met[, "idlong"] <- apply(met[, idcomblong], 1, function(x) paste(x, collapse = " "))
  
  if(format == "wide"){
    met.wide <- reshape(met, timevar = timevar, idvar = idvar, direction = "wide",
                        v.names =  c(timevar,  protocolvar, "date", dbf.file.var, "omit" ,"reason", "idlong", incubationvar))
    
    col.sort <- c(c(idvar, namevar, drugvar), colnames(met.wide) %w/o% 
                    c(idvar, namevar, drugvar))
    met.wide <- met.wide[, col.sort]
  }else{
    met.wide <- met
  }
  
  
  for(col in colnames(met.wide)){
    if(class(met.wide[, col])[1] == "factor")
      met.wide[, col] <- as.character(met.wide[, col])
  }
  
  if(!is.null(namesep)){
    spl <- strsplit(met.wide[, namevar], split = namesep) 
    cols <- c(namevar2, namecols)
    for(i in 1:length(cols))
      met.wide[, cols[i]] <- sapply(spl, function(x) x[i])
  }else{
    namevar2 <- namevar
  }
  if(!is.null(additional.metadata)){
    rownames(additional.metadata) <- additional.metadata[,namevar, drop = T]
    cols <-colnames(additional.metadata) %w/o% namevar
    met.wide[, cols] <-  additional.metadata[met.wide[, namevar2], cols]
  }
  if(format == "long"){
    met.wide[, timevar] <- as.numeric(met.wide[,timevar])
  }else{
    timevars <- colnames(met.wide)[grepl(timevar, colnames(met.wide))]
    met.wide[, timevars] <- apply(met.wide[,timevars ], 2, as.numeric)
  }
  #############################################
  ##
  ##     Merging the old and new dataset
  ##
  #############################################
  time.points.keep <- unique(as.character(met[, timevar]))
  
  if(prior){ 
    if(format == "wide"){
      old.times   <- colnames(metadata.xls)[grepl(dbf.file.var, 
                                                  colnames(metadata.xls))]
      old.times   <- gsub(dbf.file.var, "", old.times)
      rem.times   <- old.times %w/o% time.points.keep
      
      rep.vec <- c(timevar, protocolvar, "date", 
                   dbf.file.var, "omit" ,"reason", "idlong")
      rep.vec <- paste(rep(rep.vec, length(rem.times)), ".", 
                       rep(rem.times, each = length(rep.vec)), sep ="")
      
      
      metadata.xls <- metadata.xls[, colnames(metadata.xls) %w/o% rep.vec]
      
      
      row.names(met.wide) <- met.wide[, idvar]
      row.names(metadata.xls) <- metadata.xls[, idvar, drop = T]
      
      rem.id <- intersect(metadata.xls[, idvar, drop = T], met.wide[, idvar, drop = T])
      metadata.xls <- metadata.xls[rem.id, ]
      
      new.id <- met.wide[, idvar] %w/o% metadata.xls[, idvar]
      if(length(new.id))
        metadata.xls[new.id, idvar] <- new.id
      
      new.col <- colnames(met.wide)  %w/o% colnames(metadata.xls)
      for(col in new.col)
        metadata.xls[new.id, col] <- met.wide[new.id, col]
      
      int.cols <- intersect(colnames(metadata.xls), colnames(met.wide))
      if(length(new.id))
        metadata.xls[new.id, int.cols] <- met.wide[new.id, int.cols]
      
      omits <- colnames(metadata.xls)[grepl("omit.", colnames(metadata.xls))]
      reasons <- colnames(metadata.xls)[grepl("reason", colnames(metadata.xls))]
      incubations <- colnames(metadata.xls)[grepl(incubationvar, colnames(metadata.xls))]
      
      int.cols <- c(intersect(colnames(metadata.xls), colnames(met.wide)) %w/o%
                      c(omits , reasons, "unit"))
      
      
      if(!incubationvar %in% colnames)
        int.cols <- int.cols %w/o% incubations
      #for(col in int.cols) 
      #  metadata.xls[met.wide[, "idwide"] , col] <- met.wide[, col]
      
      
      metadata.xls[met.wide[, idvar] , int.cols] <- met.wide[, int.cols]
      
      time.points <- new.col[grepl("date.", new.col)]
      
      for(i in time.points){ 
        metadata.xls[, i] <- met.wide[metadata.xls[, "idwide"], i]
      }
      
      metadata.new <- metadata.full <- metadata.xls
      
      for(i in time.points.keep){
        row <- metadata.new[, paste("omit." , i, sep = "")]
        row[is.na(row)] <- FALSE
        
        col <- paste(dbf.file.var, i, sep = ".")
        
        metadata.new[row, col] <- ""
        
        metadata.new[, paste(timevar, "." , i, sep = "")] <- 
          as.numeric(metadata.new[, paste(timevar, "." , i, sep = "")])
        
      }
    }else{
      
      row.names(met.wide) <- met.wide[, idvar]
      row.names(metadata.xls) <- metadata.xls[, idvar, drop = T]
      
      rem.id <- intersect(metadata.xls[, idvar], met.wide[, idvar])
      metadata.xls <- metadata.xls[rem.id, ]
      
	    
      new.id <- met.wide[, idvar] %w/o% metadata.xls[, idvar]
	  
	  if(length(new.id))
        metadata.xls[new.id, idvar] <- new.id
      
      new.col <- colnames(met.wide)  %w/o% colnames(metadata.xls)
      for(col in new.col)
        metadata.xls[new.id, col] <- met.wide[new.id, col]
      
      int.cols <- intersect(colnames(metadata.xls), colnames(met.wide))
      
      if(length(new.id))
        metadata.xls[new.id, int.cols] <- met.wide[new.id, int.cols]
      
      int.cols <- c(intersect(colnames(metadata.xls), colnames(met.wide)) %w/o%
                      c("omit" , "reason", "unit"))
      
      if(!incubationvar %in% colnames)
        int.cols <- int.cols %w/o% incubationvar
      
      metadata.xls[met.wide[, idvar] , int.cols] <- met.wide[, int.cols]
      
      metadata.new <- metadata.full <- metadata.xls
      
      row <- metadata.new[,"omit"]
      row[is.na(row)] <- FALSE
      
      metadata.new[row, dbf.file.var] <- ""
      
      metadata.new[, timevar] <- as.numeric(metadata.new[, timevar])
    }
    
    metadata.correction <- metadata.new[metadata.new[, namevar] == correctionname, ]
    
    metadata.new <- metadata.new[ metadata.new[, namevar2] != correctionname, ]
    
  }else{
    metadata.new <- metadata.full <- met.wide
    
    if(format == "wide"){
      for(i in time.points.keep){
        row <- metadata.new[, paste("omit." , i, sep = "")]
        row[is.na(row)] <- FALSE
        
        col <- paste(dbf.file.var, i, sep = ".")
        
        metadata.new[row, col] <- ""
        
        metadata.new[, paste(timevar, "." , i, sep = "")] <- 
          as.numeric(metadata.new[, paste(timevar, "." , i, sep = "")])
        
      }
    }else{
      row <- metadata.new[,"omit"]
      row[is.na(row)] <- FALSE
      
      metadata.new[row, dbf.file.var] <- ""
      
      metadata.new[, timevar] <- as.numeric(metadata.new[, timevar])
    }
    
    metadata.correction <- metadata.new[metadata.new[, namevar] == correctionname, ]
    
    metadata.new <- metadata.new[ metadata.new[, namevar2] != correctionname, ]
  }
  
  
  #########################
  ##
  ##    Create file.info  
  ##
  ##########################
  
  dbf.files <- dir(dbf.path, pattern = file.extension, full.names=TRUE)
  
  dbf.info <- data.frame(
    file = basename(dbf.files),
    file.path = dirname(dbf.files) ,                          
    file.full = dbf.files,
    last.change = file.info(dbf.files)$mtime, 
    stringsAsFactors=FALSE)
  
  row.names(dbf.info) <- no.extension.vec(basename(dbf.files))
  
  protocol.files <- dir(protocol.path, pattern = ".xls", full.names=TRUE)
  
  protocol.info <- data.frame(
    file = basename(protocol.files),
    file.path = dirname(protocol.files) ,                          
    file.full = protocol.files, 
    last.change = file.info(protocol.files)$mtime, 
    stringsAsFactors=FALSE)
  
  row.names(protocol.info) <- no.extension.vec(basename(protocol.files))
  
  if(!prior){
    file.info <- list(dbf = dbf.info , protocol = protocol.info)
  }else{
    file.info <- old$auxiliary$file.info
    file.info$dbf <- dbf.info
    file.info$protocol = protocol.info
  }
  
  #############################################
  ##
  ##     Saving the data
  ##
  #############################################
  meta.list <- list(metadata.full       = metadata.full, 
                    metadata.correction = metadata.correction, 
                    metadata.new        = metadata.new,
                    additional          = additional.metadata)
  if(prior){
    old$call$record <- callNumbering("createMetaData", record = old$call$record)
    old$meta.list$metadata.full <- meta.list$metadata.full
    old$meta.list$metadata.correction <- meta.list$metadata.correction
    old$meta.list$metadata.new <- meta.list$metadata.new
    old$meta.list$additional <- meta.list$additional
    
    old$call[["createMetaDataDir"]] <- this.call
    old$auxiliary$file.info <- file.info
    old$auxiliary$shiny.input <- shiny.input
  }else{
    call <- list()
    call[["createMetaDataDir"]] <- this.call
    call$record <- callNumbering("createMetaDataDir")
    auxiliary <- list(file.info = file.info)
    auxiliary$shiny.input <- shiny.input
    old <- list(meta.list = meta.list, call = call,
                auxiliary = auxiliary)
  }
  
  if(is.null(additional.metadata)){
    additional.metadata <- data.frame(unique(metadata.full[, namevar]),
                                      stringsAsFactors=FALSE)
    colnames(additional.metadata) <- namevar
  }
  
  additional.metadata <- additional.metadata[order(additional.metadata[, namevar]),, drop = FALSE ]
  
  meta.sheets[["Data"]] <- metadata.full
  meta.sheets[["additional.metadata"]] <- additional.metadata
  write_xlsx(meta.sheets, 
			 path=  xls.file, 
             format_headers = TRUE,
		     col_names = TRUE)
  
  if(format == "wide"){
    data <- metadata.full
    if(idvar == "id")
      data$id2 <- data[,idvar]
    if(timevar == "time")
      data$time2 <- data[,timevar]
    
    
    data.long <- 
      reshape(data, direction = "long", 
              varying = list(
                protocolvar = colnames(data)[
                  grepl(paste(protocolvar,".",sep = ""), colnames(data))], 
                timevar = colnames(data)[
                  grepl(paste(timevar,".",sep = ""), colnames(data))],
                date = colnames(data)[
                  grepl("date.", colnames(data))],
                dbf.file.var = colnames(data)[
                  grepl(dbf.file.var, colnames(data))],
                omit = colnames(data)[
                  grepl("omit.", colnames(data))],
                reason = colnames(data)[
                  grepl("reason.", colnames(data))],
                incubation = colnames(data)[
                  grepl(paste(incubationvar, ".", sep = ""), colnames(data))],
                idlong = colnames(data)[
                  grepl("idlong.", colnames(data))] 
              ),
              v.names = c(protocolvar, timevar,"date", 
                          dbf.file.var, "omit", "reason",incubationvar,
                          "idlong"))
    
    
    
    data.long <- data.long[!is.na(data.long[dbf.file.var]),]
    
    if(idvar == "id"){
      data.long$id <- data.long$id2
      data.long <- data.long[,colnames(data.long) != "id2"]
    }else{
      data.long <- data.long[,colnames(data.long) != "id"]
    }
    if(timevar == "time"){
      data.long$time <- data.long$time2
      data.long <- data.long[,colnames(data.long) != "time2"]
    }else{
      data.long <- data.long[, colnames(data.long) != "time"]
    }
    
    
    idcombwide  <- c(namevar, drugvar, identifier)
    idcomb      <- c(namevar, drugvar, identifier)
    
    idcomblong  <- c(namevar, drugvar, identifier, timevar)
    
    # data <- data.long[data.long$idlong != "", ]  
    data.long[, idvar] <- apply(data.long[, idcomblong], 1, 
                                function(x) paste(x, collapse = " "))
    if(!("idwide" %in% colnames(data.long)))
      data.long[, "idwide"] <- apply(data.long[, idcombwide], 1, 
                                     function(x) paste(x, collapse = " "))
    if(!("idlong" %in% colnames(data.long)))
      data.long[, "idlong"] <- apply(data.long[, idcomblong], 1, 
                                     function(x) paste(x, collapse = " "))
    
    data <- data.long[data.long$idlong != "" &
                        !is.na(data.long$idlong ), ]  
    
    
    
    old$meta.list$metadata.long <- data.long
  }else{
    old$meta.list$metadata.long <- metadata.full
  }
  
  class(old) <- c("createMetaData", class(old) %w/o% "createMetaData")
  if(save)
  save(old, file = data.file)
  return(old)
}

##
##
## Create metadata from an Excel sheet
##


createMetaDataXls <- function(
  excel.file       = NULL,
  dbf.path         = getwd(),
  protocol.path    =  getwd(),
  namevar          = "Cellline",
  drugvar          = "chemo",
  protocolvar      = "R.protocol",
  dbf.file.var     = "dbf.file.name",
  identifier       = "identifier",
  timevar          = "Hour",
  additional.metadata   = NULL, # data.frame with extra information for the cell lines
  # either excel or data.frame
  incubationvar    = "incubation",
  doublingvar       = NULL,
  correctionname   = "Control",
  unit             = "ug/ml",
  idvar            = "sampleid",
  show.warnings    = TRUE,
  update           = TRUE,
  format           = "long",
  save            = TRUE,
  data.file      = file.path(getwd(), "Absorbance"),
  shiny.input      = NULL){
  
  #############################################
  ##
  ##     Function starts 
  ##
  #############################################
  call2 <- match.call()
  this.call <- list()
  
  myfor <- formals(createMetaDataXls)               ## formals with default arguments
  for ( v in names(myfor)){
    if (!(v %in% names(call2))){
      this.call[[v]] <- myfor[[v]]
    }else{
      this.call[[v]] <- call2[[v]]
    }  ## if arg is missing I add it
  }
  this.call$call <- call2
  
  
  ######### old data is loaded if available
  
  data.file   <- paste(data.file, ".RData", sep = "")
  prior <- all(file.exists(data.file), update)
  if(prior)
    load(data.file)
  
  if(is.null( excel.file ) & prior)
    if("createMetaData" %in% old$call)
      excel.file <- old$call[["createMetaData"]][["excel.file"]]
  if(is.null( excel.file ))
    excel.file <- file.choose()
  
  metadata <- 
    data.frame(read_excel(excel.file))
  
  
  if(class(additional.metadata) == "character"){
    additional.metadata <- 
      data.frame(read_excel(additional.metadata))
  }
  
  if(!is.null(additional.metadata)){
    if(!(namevar %in% colnames(additional.metadata)))
      stop('the column name "', namevar, 
           '" is not included in the supplied external metadata')
  }
  
  if(format == "wide"){
    idcomb  <- c(namevar, drugvar, identifier)
  }else{
    idcomb  <- c(namevar, drugvar, identifier, timevar)
  }
  
  idcombwide  <- c(namevar, drugvar, identifier)
  
  idcomblong  <- c(namevar, drugvar, identifier, timevar)
  
  
  metadata[, idvar] <- apply(metadata[, idcomb], 1, function(x) paste(x, collapse = " "))
  
  metadata[, "idwide"] <- apply(metadata[, idcombwide], 1, function(x) paste(x, collapse = " "))
  metadata[, "idlong"] <- apply(metadata[, idcomblong], 1, function(x) paste(x, collapse = " "))
  
  if(!is.null(additional.metadata)){
    rownames(additional.metadata) <- additional.metadata[,namevar]
    cols <-colnames(additional.metadata) %w/o% namevar
    metadata[, cols] <-  additional.metadata[metadata[, namevar], cols]
  }
  
  if(!"unit" %in% colnames(metadata)){
    metadata[, "unit"] <- unit
  }else{
    if(any(is.na( metadata[, "unit"])))
    metadata[is.na( metadata[, "unit"]), "unit"] <- unit
  }
  metadata[,timevar] <- as.numeric(metadata[,timevar])
  metadata.full <- metadata
  
  metadata.correction <- metadata[metadata[, namevar] == correctionname, ]
  
  metadata.new <- metadata[ metadata[, namevar] != correctionname, ]
  
  

  #############################################
  ##
  ##     Saving the data
  ##
  #############################################
  meta.list <- list(metadata.full       = metadata.full, 
                    metadata.correction = metadata.correction, 
                    metadata.new        = metadata.new,
                    additional          = additional.metadata)
  if(prior){
    old$call$record <- callNumbering("createMetaData", record = old$call$record)
    
    old$meta.list$metadata.full <- meta.list$metadata.full
    old$meta.list$metadata.correction <- meta.list$metadata.correction
    old$meta.list$metadata.new <- meta.list$metadata.new
    old$meta.list$additional <- meta.list$additional
    
    old$call[["createMetaData"]] <- this.call
    old$auxiliary$file.info <- file.info
    old$auxiliary$shiny.input <- shiny.input
  }else{
    call <- list()
    call[["createMetaData"]] <- this.call
    call$record <- callNumbering("createMetaData")
    auxiliary <- list(file.info = file.info)
    auxiliary$shiny.input <- shiny.input
    old <- list(meta.list = meta.list, call = call,
                auxiliary = auxiliary)
  }
  
  ###############################################
  ##
  ##   Creating the metadata.long
  ##
  #################################################
  
  if(format == "wide"){
    data <- metadata.full
    if(idvar == "id")
      data$id2 <- data[,idvar]
    if(timevar == "time")
      data$time2 <- data[,timevar]
    
    
    data.long <- 
      reshape(data, direction = "long", 
              varying = list(
                protocolvar = colnames(data)[
                  grepl(paste(protocolvar,".",sep = ""), colnames(data))], 
                timevar = colnames(data)[
                  grepl(paste(timevar,".",sep = ""), colnames(data))],
                date = colnames(data)[
                  grepl("date.", colnames(data))],
                dbf.file.var = colnames(data)[
                  grepl(dbf.file.var, colnames(data))],
                omit = colnames(data)[
                  grepl("omit.", colnames(data))],
                reason = colnames(data)[
                  grepl("reason.", colnames(data))],
                incubation = colnames(data)[
                  grepl(paste(incubationvar, ".", sep = ""), colnames(data))],
                idlong = colnames(data)[
                  grepl("idlong.", colnames(data))] 
              ),
              v.names = c(protocolvar, timevar,"date", 
                          dbf.file.var, "omit", "reason",incubationvar,
                          "idlong"))
    
    
    
    data.long <- data.long[!is.na(data.long[dbf.file.var]),]
    
    if(idvar == "id"){
      data.long$id <- data.long$id2
      data.long <- data.long[,colnames(data.long) != "id2"]
    }else{
      data.long <- data.long[,colnames(data.long) != "id"]
    }
    if(timevar == "time"){
      data.long$time <- data.long$time2
      data.long <- data.long[,colnames(data.long) != "time2"]
    }else{
      data.long <- data.long[, colnames(data.long) != "time"]
    }
    
    
    idcombwide  <- c(namevar, drugvar, identifier)
    idcomb      <- c(namevar, drugvar, identifier)
    
    idcomblong  <- c(namevar, drugvar, identifier, timevar)
    
    # data <- data.long[data.long$idlong != "", ]  
    data.long[, idvar] <- apply(data.long[, idcomblong], 1, 
                                function(x) paste(x, collapse = " "))
    if(!("idwide" %in% colnames(data.long)))
      data.long[, "idwide"] <- apply(data.long[, idcombwide], 1, 
                                     function(x) paste(x, collapse = " "))
    if(!("idlong" %in% colnames(data.long)))
      data.long[, "idlong"] <- apply(data.long[, idcomblong], 1, 
                                     function(x) paste(x, collapse = " "))
    
    data <- data.long[data.long$idlong != "" &
                        !is.na(data.long$idlong ), ]  
    
    
    
    old$meta.list$metadata.long <- data.long
  }else{
    old$meta.list$metadata.long <- metadata.full
  }
  
  
  class(old) <- c("createMetaData", class(old) %w/o% "createMetaData")
  if(save)
  save(old, file = data.file)
  return(old)
}



createMetaDataList <- function(
  list             = NULL,  # Absorbance, data.frame, list 
  dbf.files        = NULL,
  protocol.files   = NULL,
  are.paths.full    = TRUE, 
  dbf.path         = getwd(),
  protocol.path    = getwd(),
  namevar          = "Cellline",
  drugvar          = "chemo",
  protocolvar      = "R.protocol",
  dbf.file.var     = "dbf.file.name",
  identifier       = "identifier",
  timevar          = "Hour",
  additional.metadata = NULL, # data.frame with extra information for the cell line.
  incubationvar    = "incubation",
  doublingvar       = NULL,
  correctionname   = "Control",
  unit             = "ug/ml",
  idvar            = "sampleid",
  format           = "long",
  file.extension   = "dbf",
  show.warnings    = TRUE,
  update           = TRUE,
  save             = TRUE,
  data.file      = file.path(getwd(), "Absorbance"),
  shiny.input      = NULL)
{
  
  #############################################
  ##
  ##     Function starts 
  ##
  #############################################
  call2 <- match.call()
  this.call <- list()
  
  myfor <- formals(createMetaDataList)               ## formals with default arguments
  for ( v in names(myfor)){
    if (!(v %in% names(call2))){
      this.call[[v]] <- myfor[[v]]
    }else{
      this.call[[v]] <- call2[[v]]
    }  ## if arg is missing I add it
  }
  this.call$call <- call2  
  
  ######### old data is loaded if available
  
  data.file   <- paste(data.file, ".RData", sep = "")
  prior <- all(file.exists(data.file), update)
  if(prior)
    load(data.file)
  
  
  if(class(additional.metadata) == "character"){
    additional.metadata <- 
      data.frame(read_excel(additional.metadata))
  }
  
  if(!is.null(additional.metadata)){
    if(!(namevar %in% colnames(additional.metadata)))
      stop('the column name "', namevar, 
           '" is not included in the supplied external metadata')
  }
  
  
  #############################################
  ##
  ##     Create the metadata if 
  ##     class(data)[1] == list()
  ##
  #############################################
  
  metadata <- as.data.frame(list, stringsAsFactors=FALSE)    
  
  metadata[, dbf.file.var] <- no.extension.vec(basename(dbf.files))
  
  metadata[, protocolvar]  <- no.extension.vec(basename(protocol.files))
  
  idcombwide  <- c(namevar, drugvar, identifier)
  idcomb  <- c(namevar, drugvar, identifier)
  
  idcomblong  <- c(namevar, drugvar, identifier, timevar)
  
  
  metadata[, idvar] <- apply(metadata[, idcomb], 1, 
                             function(x) paste(x, collapse = " "))
  if(!("idwide" %in% colnames(metadata)))
    metadata[, "idwide"] <- apply(metadata[, idcombwide], 1, 
                                  function(x) paste(x, collapse = " "))
  if(!("idlong" %in% colnames(metadata)))
    metadata[, "idlong"] <- apply(metadata[, idcomblong], 1, 
                                  function(x) paste(x, collapse = " "))
  
  metadata <- metadata[metadata$idlong != "" & !is.na(metadata$idlong ), ]  
  
  if(!"unit" %in% colnames(metadata)){
    metadata[, "unit"] <- unit
  }else{
    if(any(is.na( metadata[, "unit"])))
      metadata[is.na( metadata[, "unit"]), "unit"] <- unit
  }
  
  #############################################
  ##
  ##    Create the file.info list
  ##    This list keeps track of the files
  ##    It is used for reading the files
  ##    and changing the files
  ##
  #############################################
  
  
  if(!are.paths.full & !is.null(dbf.files))
    dbf.files <- file.path(dbf.path, dbf.files)
  
  if(!are.paths.full & !is.null(protocol.files))
    protocol.files <- file.path(protocol.path, protocol.files)
  
  if(is.null(dbf.files))    
    dbf.files <- dir(dbf.path, pattern = file.extension, full.names=TRUE)
  
  dbf.info <- data.frame(
    file = basename(dbf.files),
    file.path = dirname(dbf.files) ,                          
    file.full = dbf.files, 
    last.change = file.info(dbf.files)$mtime,
    stringsAsFactors=FALSE)
  
  row.names(dbf.info) <- no.extension.vec(basename(dbf.files))
  
  
  if(is.null(protocol.files))    
    protocol.files <- dir(protocol.path, pattern = ".xls|.xlsx", full.names=TRUE)
  
  protocol.files <- unique(protocol.files)
  
  protocol.info <- data.frame(
    file = basename(protocol.files),
    file.path = dirname(protocol.files) ,                          
    file.full = protocol.files, 
    last.change = file.info(protocol.files)$mtime,
    stringsAsFactors=FALSE)
  
  row.names(protocol.info) <- no.extension.vec(basename(protocol.files))
  
  if(!prior){
    file.info <- list(dbf = dbf.info , protocol = protocol.info)
  }else{
    file.info <- old$auxiliary$file.info
    file.info$dbf <- dbf.info
    file.info$protocol = protocol.info
  }
  
  
  ###########################################
  ##
  ##  Identifiers added
  ##
  ###########################################
  
  idcomb  <- c(namevar, drugvar, identifier, timevar)
  
  
  idcombwide  <- c(namevar, drugvar, identifier)
  
  idcomblong  <- c(namevar, drugvar, identifier, timevar)
  
  
  metadata[, idvar] <- apply(metadata[, idcomb], 1, function(x) paste(x, collapse = " "))
  
  metadata[, "idwide"] <- apply(metadata[, idcombwide], 1, function(x) paste(x, collapse = " "))
  metadata[, "idlong"] <- apply(metadata[, idcomblong], 1, function(x) paste(x, collapse = " "))
  
  if(!is.null(additional.metadata)){
    rownames(additional.metadata) <- additional.metadata[,namevar]
    cols <-colnames(additional.metadata) %w/o% namevar
    metadata[, cols] <-  additional.metadata[metadata[, namevar], cols]
  }
  
  metadata[,timevar] <- as.numeric(metadata[,timevar])
  
  metadata.full <- metadata
  
  metadata.correction <- metadata[metadata[, namevar] == correctionname, ]
  
  metadata.new <- metadata[ metadata[, namevar] != correctionname, ]
    
  #############################################
  ##
  ##     Saving the data
  ##
  #############################################
  meta.list <- list(metadata.full       = metadata.full, 
                    metadata.correction = metadata.correction, 
                    metadata.new        = metadata.new,
                    additional          = additional.metadata)
  
  if(prior){
    old$call$record <- callNumbering("createMetaData", record = old$call$record)
    old$meta.list$metadata.full <- meta.list$metadata.full
    old$meta.list$metadata.correction <- meta.list$metadata.correction
    old$meta.list$metadata.new <- meta.list$metadata.new
    old$meta.list$additional <- meta.list$additional
    old$call[["createMetaData"]] <- this.call
    old$auxiliary$file.info <- file.info
    old$auxiliary$shiny.input <- shiny.input
  }else{  
    call <- list()
    call[["createMetaData"]] <- this.call
    call$record <- callNumbering("createMetaData")
    auxiliary <- list(file.info = file.info)
    auxiliary$shiny.input <- shiny.input
    old <- list(meta.list = meta.list, call = call,
                auxiliary = auxiliary)
  }
  
  ###############################################
  ##
  ##   Creating the metadata.long
  ##
  #################################################
  
  if(format == "wide"){
    data <- metadata.full
    if(idvar == "id")
      data$id2 <- data[,idvar]
    if(timevar == "time")
      data$time2 <- data[,timevar]
    
    
    data.long <- 
      reshape(data, direction = "long", 
              varying = list(
                protocolvar = colnames(data)[
                  grepl(paste(protocolvar,".",sep = ""), colnames(data))], 
                timevar = colnames(data)[
                  grepl(paste(timevar,".",sep = ""), colnames(data))],
                date = colnames(data)[
                  grepl("date.", colnames(data))],
                dbf.file.var = colnames(data)[
                  grepl(dbf.file.var, colnames(data))],
                omit = colnames(data)[
                  grepl("omit.", colnames(data))],
                reason = colnames(data)[
                  grepl("reason.", colnames(data))],
                incubation = colnames(data)[
                  grepl(paste(incubationvar, ".", sep = ""), colnames(data))],
                idlong = colnames(data)[
                  grepl("idlong.", colnames(data))] 
              ),
              v.names = c(protocolvar, timevar,"date", 
                          dbf.file.var, "omit", "reason",incubationvar,
                          "idlong"))
    
    
    
    data.long <- data.long[!is.na(data.long[dbf.file.var]),]
    
    if(idvar == "id"){
      data.long$id <- data.long$id2
      data.long <- data.long[,colnames(data.long) != "id2"]
    }else{
      data.long <- data.long[,colnames(data.long) != "id"]
    }
    if(timevar == "time"){
      data.long$time <- data.long$time2
      data.long <- data.long[,colnames(data.long) != "time2"]
    }else{
      data.long <- data.long[, colnames(data.long) != "time"]
    }
    
    
    idcombwide  <- c(namevar, drugvar, identifier)
    idcomb      <- c(namevar, drugvar, identifier)
    
    idcomblong  <- c(namevar, drugvar, identifier, timevar)
    
    # data <- data.long[data.long$idlong != "", ]  
    data.long[, idvar] <- apply(data.long[, idcomblong], 1, 
                                function(x) paste(x, collapse = " "))
    if(!("idwide" %in% colnames(data.long)))
      data.long[, "idwide"] <- apply(data.long[, idcombwide], 1, 
                                     function(x) paste(x, collapse = " "))
    if(!("idlong" %in% colnames(data.long)))
      data.long[, "idlong"] <- apply(data.long[, idcomblong], 1, 
                                     function(x) paste(x, collapse = " "))
    
    data <- data.long[data.long$idlong != "" &
                        !is.na(data.long$idlong ), ]  
    
    
    
    old$meta.list$metadata.long <- data.long
  }else{
    old$meta.list$metadata.long <- metadata.full
  }

  
  class(old) <- c("createMetaData", class(old) %w/o% "createMetaData")
  if(save)
  save(old, file = data.file)
  return(old)
}

##
##
## Create metadata from an data.frame
##
##


createMetaDataDataFrame <- function(
  metadata         = NULL,
  dbf.path         = getwd(),
  protocol.path    = getwd(),
  namevar          = "Cellline",
  drugvar          = "chemo",
  protocolvar      = "R.protocol",
  dbf.file.var     = "dbf.file.name",
  identifier       = "identifier",
  timevar          = "Hour",
  additional.metadata = NULL, # data.frame with extra information for the cell lines
  incubationvar    = "incubation",
  doublingvar       = NULL,
  correctionname   = "Control",
  unit             = "ug/ml",
  idvar            = "sampleid",
  show.warnings    = TRUE,
  update           = TRUE,
  format           = "long",
  file.extension   = "dbf",
  save             = TRUE,
  data.file      = file.path(getwd(), "Absorbance"),
  shiny.input      = NULL)
{
  
  #############################################
  ##
  ##     Function starts 
  ##
  #############################################
  call2 <- match.call()
  this.call <- list()
  
  myfor <- formals(createMetaDataDataFrame)               ## formals with default arguments
  for ( v in names(myfor)){
    if (!(v %in% names(call2))){
      this.call[[v]] <- myfor[[v]]
    }else{
      this.call[[v]] <- call2[[v]]
    }  ## if arg is missing I add it
  }
  this.call$call <- call2
  
  
  ######### old data is loaded if available
  
  data.file   <- paste(data.file, ".RData", sep = "")
  prior <- all(file.exists(data.file), update)
  if(prior)
    load(data.file)
  
  
  if(class(additional.metadata) == "character"){
    additional.metadata <- 
      data.frame(read_excel(additional.metadata))
  }
  
  if(!is.null(additional.metadata)){
    if(!(namevar %in% colnames(additional.metadata)))
      stop('the column name "', namevar, 
           '" is not included in the supplied external metadata')
  }
  
  if(format == "wide"){
    idcomb  <- c(namevar, drugvar, identifier)
  }else{
    idcomb  <- c(namevar, drugvar, identifier, timevar)
  }
  
  idcombwide  <- c(namevar, drugvar, identifier)
  
  idcomblong  <- c(namevar, drugvar, identifier, timevar)
  
  
  metadata[, idvar] <- apply(metadata[, idcomb], 1, function(x) paste(x, collapse = " "))
  
  metadata[, "idwide"] <- apply(metadata[, idcombwide], 1, function(x) paste(x, collapse = " "))
  metadata[, "idlong"] <- apply(metadata[, idcomblong], 1, function(x) paste(x, collapse = " "))
  
  if(!is.null(additional.metadata)){
    rownames(additional.metadata) <- additional.metadata[,namevar]
    cols <-colnames(additional.metadata) %w/o% namevar
    metadata[, cols] <-  additional.metadata[metadata[, namevar], cols]
  }
  
  if(!"unit" %in% colnames(metadata)){
    metadata[, "unit"] <- unit
  }else{
    if(any(is.na( metadata[, "unit"])))
      metadata[is.na( metadata[, "unit"]), "unit"] <- unit
  }
  
  metadata[,timevar] <- as.numeric(metadata[,timevar])
  
  metadata.full <- metadata
  
  metadata.correction <- metadata[metadata[, namevar] == correctionname, ]
  
  metadata.new <- metadata[ metadata[, namevar] != correctionname, ]
  
  #############################################
  ##
  ##    Create the file.info list
  ##    This list keeps track of the files
  ##    It is used for reading the files
  ##    and changing the files
  ##
  #############################################
  
  
  dbf.files <- dir(dbf.path, pattern = file.extension, full.names=TRUE)
  
  dbf.info <- data.frame(
    file = basename(dbf.files),
    file.path = dirname(dbf.files) ,                          
    file.full = dbf.files, 
    last.change = file.info(dbf.files)$mtime,
    stringsAsFactors=FALSE)
  
  row.names(dbf.info) <- no.extension.vec(basename(dbf.files))
  
  protocol.files <- dir(protocol.path, pattern = ".xls|.xlsx", full.names=TRUE)
  
  protocol.files <- unique(protocol.files)
  
  protocol.info <- data.frame(
    file = basename(protocol.files),
    file.path = dirname(protocol.files) ,                          
    file.full = protocol.files, 
    last.change = file.info(protocol.files)$mtime,
    stringsAsFactors=FALSE)
  
  row.names(protocol.info) <- no.extension.vec(basename(protocol.files))
  
  if(!prior){
    file.info <- list(dbf = dbf.info , protocol = protocol.info)
  }else{
    file.info <- old$auxiliary$file.info
    file.info$dbf <- dbf.info
    file.info$protocol = protocol.info
  }
  
  #############################################
  ##
  ##     creating the object old
  ##
  #############################################
  meta.list <- list(metadata.full       = metadata.full, 
                    metadata.correction = metadata.correction, 
                    metadata.new        = metadata.new,
                    additional          = additional.metadata)
  
  if(prior){
    old$call$record <- callNumbering("createMetaData", record = old$call$record)
    old$meta.list$metadata.full <- meta.list$metadata.full
    old$meta.list$metadata.correction <- meta.list$metadata.correction
    old$meta.list$metadata.new <- meta.list$metadata.new
    old$meta.list$additional <- meta.list$additional
    
    old$call[["createMetaData"]] <- this.call
    old$auxiliary$file.info <- file.info
    old$auxiliary$shiny.input <- shiny.input
  }else{
    call <- list()
    call[["createMetaData"]] <- this.call
    call$record  <- callNumbering("createMetaData")
    auxiliary <- list(file.info = file.info)
    auxiliary$shiny.input <- shiny.input
    old <- list(meta.list = meta.list, call = call,
                auxiliary = auxiliary)
  }
  
  
  ###############################################
  ##
  ##   Creating the metadata.long
  ##
  #################################################
  
  if(format == "wide"){
    data <- metadata.full
    if(idvar == "id")
      data$id2 <- data[,idvar]
    if(timevar == "time")
      data$time2 <- data[,timevar]
    
    
    data.long <- 
      reshape(data, direction = "long", 
              varying = list(
                protocolvar = colnames(data)[
                  grepl(paste(protocolvar,".",sep = ""), colnames(data))], 
                timevar = colnames(data)[
                  grepl(paste(timevar,".",sep = ""), colnames(data))],
                date = colnames(data)[
                  grepl("date.", colnames(data))],
                dbf.file.var = colnames(data)[
                  grepl(dbf.file.var, colnames(data))],
                omit = colnames(data)[
                  grepl("omit.", colnames(data))],
                reason = colnames(data)[
                  grepl("reason.", colnames(data))],
                incubation = colnames(data)[
                  grepl(paste(incubationvar, ".", sep = ""), colnames(data))],
                idlong = colnames(data)[
                  grepl("idlong.", colnames(data))] 
              ),
              v.names = c(protocolvar, timevar,"date", 
                          dbf.file.var, "omit", "reason",incubationvar,
                          "idlong"))
    
    
    
    data.long <- data.long[!is.na(data.long[dbf.file.var]),]
    
    if(idvar == "id"){
      data.long$id <- data.long$id2
      data.long <- data.long[,colnames(data.long) != "id2"]
    }else{
      data.long <- data.long[,colnames(data.long) != "id"]
    }
    if(timevar == "time"){
      data.long$time <- data.long$time2
      data.long <- data.long[,colnames(data.long) != "time2"]
    }else{
      data.long <- data.long[, colnames(data.long) != "time"]
    }
    
    
    idcombwide  <- c(namevar, drugvar, identifier)
    idcomb      <- c(namevar, drugvar, identifier)
    
    idcomblong  <- c(namevar, drugvar, identifier, timevar)
    
    # data <- data.long[data.long$idlong != "", ]  
    data.long[, idvar] <- apply(data.long[, idcomblong], 1, 
                                function(x) paste(x, collapse = " "))
    if(!("idwide" %in% colnames(data.long)))
      data.long[, "idwide"] <- apply(data.long[, idcombwide], 1, 
                                     function(x) paste(x, collapse = " "))
    if(!("idlong" %in% colnames(data.long)))
      data.long[, "idlong"] <- apply(data.long[, idcomblong], 1, 
                                     function(x) paste(x, collapse = " "))
    
    data <- data.long[data.long$idlong != "" &
                        !is.na(data.long$idlong ), ]  
    
    
    
    old$meta.list$metadata.long <- data.long
  }else{
    old$meta.list$metadata.long <- metadata.full
  }
  
  class(old) <- c("createMetaData", class(old) %w/o% "createMetaData")
  if(save)
    save(old, file = data.file)
  return(old)
}


A0T1fun <- function(X, timevar = NULL, control, entity, BC = "BC", dataset = "bs.mean",  namevar = "name") {
  if("bootstrap" %in% class(X) |"simult" %in% class(X)){
    bs.names <- ##as.numeric(gsub("BS:","",
      c(BC, colnames(X$data$bs.mean)[grepl("BS:", colnames(X$data$bs.mean))])#))
    data <- X$data$bs.mean
  }
  
  if(!("bootstrap" %in% class(X) |"simult" %in% class(X))){
    bs.names <- c(BC)
    data <- X$data$bc.mean
  }
  
  
  call <- as.list(X$call$readDBFData)
  type <- "new.type"
  for(times in unique(data[, timevar])){
    
    X1 <- data[data[, type] == control & data[, timevar] %in% times,
               c("plate", entity, "outlier", namevar,
                 timevar, bs.names, type)]
    
    rownames(X1) <- X1[, entity]
    
    data[, paste("A0T", times, bs.names, sep = "")] <- 
      X1[as.character(data[, entity]), bs.names, drop = FALSE]
    
  }
  xx <- list()
  for(times in unique(data[,timevar])){
    xx[[paste(times)]] <- data[data[,timevar] == times, ]
    xx[[paste(times)]][, paste("A0T", bs.names, sep = "")] <-
      xx[[paste(times)]][paste("A0T", times, bs.names, sep = "")]
  }
  data2 <- xx[[1]]
  if(length(xx) > 1 )
    for(i in 2:length(xx))
      data2 <- rbind(data2, xx[[i]])
  X$data$GI.mean <- data2[rownames(data),]
  return(X)
}

AC0fun <- function(X, timevar = NULL, control, entity, BC = "BC"){
  if("bootstrap" %in% class(X)){
    bs.names <- ##as.numeric(gsub("BS:","",
      c(BC, colnames(X$data$bs.mean)[grepl("BS:", colnames(X$data$bs.mean))])#))
    data <- X$data$GI.mean
  }
  if("simult" %in% class(X)){
    bs.names <- ##as.numeric(gsub("BS:","",
      c(colnames(X$data$bs.mean)[grepl("BS:", colnames(X$data$bs.mean))])#))
    data <- X$data$GI.mean
  }
  if(!("bootstrap" %in% class(X) |"simult" %in% class(X))){
    bs.names <- c(BC)
    data <- X$data$GI.mean
  }
  call <- as.list(X$call)
  
  X1 <- data[data[, timevar] %in% 0 & data$outlier == 0,
             c("new.type",  "plate", entity,
               timevar, bs.names)]
  
  rownames(X1) <- paste(X1[,entity], X1[, "new.type"])
  
  data[, paste("A0C", bs.names, sep = "")] <-
    X1[paste(data[, entity], data[, "new.type"]),
       bs.names]
  X$data$GI.mean <- data
  return(X)
}




doublingTime <- function(X, timevar = "time", control = "X1",
                         use.supplied.T0 = FALSE, doublingvar = "supT0",
                         entity = "entity", namevar = "name", BC = "BC"){
  
  if("bootstrap" %in% class(X)){
    bs.names <- ##as.numeric(gsub("BS:","",
      c(BC, colnames(X$data$bs.mean)[grepl("BS:", colnames(X$data$bs.mean))])#))
    data <- X$data$GI.mean
  }
  if("simult" %in% class(X)){
    bs.names <- ##as.numeric(gsub("BS:","",
      c(colnames(X$data$bs.mean)[grepl("BS:", colnames(X$data$bs.mean))])#))
    data <- X$data$GI.mean
  }
  if(!("bootstrap" %in% class(X) |"simult" %in% class(X))){
    bs.names <- c(BC)
    data <- X$data$GI.mean
  }
  
  call <- as.list(X$auxiliary$passed.var)
  if(min(data[, timevar]) == max(data[, timevar]))
    use.supplied.T0 <- TRUE
  if(!use.supplied.T0){
    for(bs.iter in bs.names){
      xx1 <- data[data[, timevar] ==  min(data[, timevar])&
                    data[, "new.type"] == control,][,
                                                    c(namevar, entity, bs.iter)]
      
      xx2 <- data[data[, timevar] == max(data[, timevar]) &
                    data[, "new.type"] == control,][,
                                                    c(namevar, entity, bs.iter)]
      
      xx1 <- xx1[order(xx1[, entity]), ]
      xx2 <- xx2[order(xx2[, entity]), ]
      
      dob         <- merge(xx1, xx2, by =c(entity, namevar))
      dob$T0      <- (max(data[, timevar])- min(data[, timevar])) / log2(dob[,paste(bs.iter, ".y", sep = "")]/
                                                                           dob[,paste(bs.iter, ".x", sep = "")])
      formula <- as.formula(paste("T0 ~ ", namevar))
      T0 <- aggregate(formula, data = dob, FUN = median) #skiftet tilformula
      
      rownames(T0) <- T0[,namevar]
      X$data$GI.mean[,paste("T0", bs.iter, sep = "")] <- T0[X$data$GI.mean[,namevar], "T0"]
    }
  }else{
    for(bs.iter in bs.names){
      if(doublingvar %in% colnames(X$meta.list$additional)){
        T0 <- X$meta.list$additional#[, doublingvar] #skiftet tilformula
        
        #rownames(T0) <- T0[,namevar]
        X$data$GI.mean[,paste("T0", bs.iter, sep = "")] <- T0[X$data$GI.mean[,namevar], doublingvar]
      }else{
        X$data$GI.mean[,paste("T0", bs.iter, sep = "")] <- NA
      }
    }
  }
  return(X)
}



DGrowth <- function(X, timevar = "time", BC = "BC", update = TRUE,
                    doublingvar = NULL, DG = TRUE, suplied.doubling.time = FALSE, 
                    n.samples = NULL) {
  
  if("bootstrap" %in% class(X) |"simult" %in% class(X)){
    if(is.null(n.samples)){
      bs.names <- ##as.numeric(gsub("BS:","",
        c(BC, colnames(X$data$bs.mean)[grepl("BS:", colnames(X$data$bs.mean))])#))
    }else{
      bs.names <-  c(BC, paste("BS:", 1:n.samples, sep = ""))
    }
    data <- X$data$GI.mean
  }
  if(!("bootstrap" %in% class(X) |"simult" %in% class(X))){
    bs.names <- c(BC)
    data <- X$data$GI.mean
  }
  
  for(bs.iter in bs.names){
    data[, paste("D", bs.iter, sep = ".")] <-
      ifelse(data[, bs.iter] >=
               data[, paste("A0T0", bs.iter, sep = "")],
             (data[, bs.iter] - data[, paste("A0T0", bs.iter, sep = "")]) /
               (data[, paste("A0T", bs.iter, sep = "")] -
                  data[, paste("A0T0", bs.iter, sep = "")]),
             (data[, bs.iter] - data[, paste("A0T0", bs.iter, sep = "")]) /
               data[, paste("A0T0", bs.iter, sep = "")])
    
    data[, paste("D.upper", bs.iter, sep = ".")] <- 
      (data[, bs.iter] - data[, paste("A0T0", bs.iter, sep = "")]) /
      (data[, paste("A0T", bs.iter, sep = "")] -
         data[, paste("A0T0", bs.iter, sep = "")])
    
    data[, paste("D.lower", bs.iter, sep = ".")] <-  
      (data[, bs.iter] - data[, paste("A0T0", bs.iter, sep = "")]) /
      data[, paste("A0T0", bs.iter, sep = "")]
    
    
    if(DG){
      
      data[, paste("DG.upper", bs.iter, sep = ".")] <-
        (data[, paste("T0", bs.iter, sep = "")] / data[,timevar]) * 
        log2((2^(data[,timevar] /data[, paste("T0", bs.iter, sep = "")]) -1) *
               data[, paste("D.upper", bs.iter, sep = ".")] + 1) * 100
      
      
      data[, paste("DG.lower", bs.iter, sep = ".")] <-
        1/data[,timevar]*log2(data[, paste("D.lower", bs.iter, sep = ".")] + 1)
    }    
    data[, paste("D", bs.iter, sep = ".")]  <- 
      data[, paste("D", bs.iter, sep = ".")]  * 100
    
    data[, paste("D.upper", bs.iter, sep = ".")]  <- 
      data[, paste("D.upper", bs.iter, sep = ".")]  * 100
    data[, paste("D.lower", bs.iter, sep = ".")]  <- 
      data[, paste("D.lower", bs.iter, sep = ".")]  * 100  
  }
  
  X$data$GI.mean <- data
  return(X)
}

#table(bs.mean[,timevar], bs.mean[,namevar], bs.mean[,drugvar]) 

RGrowth <- function(X, timevar  = "time", BC = "BC", update = TRUE,
                    doublingvar = NULL, RG = TRUE, 
                    n.samples = NULL) {
  
  if("bootstrap" %in% class(X) |"simult" %in% class(X)){
    if(is.null(n.samples)){
      bs.names <- ##as.numeric(gsub("BS:","",
        c(BC, colnames(X$data$bs.mean)[grepl("BS:", colnames(X$data$bs.mean))])#))
    }else{
      bs.names <-  c(BC, paste("BS:", 1:n.samples, sep = ""))
    }
    data <- X$data$GI.mean
  }
  if(!("bootstrap" %in% class(X) |"simult" %in% class(X))){
    bs.names <- c(BC)
    data <- X$data$GI.mean
  }
  
  
  
  if(!RG){
    for(bs.iter in bs.names){
      data[, paste("R", bs.iter, sep = ".")] <-
        data[, bs.iter] / data[, paste("A0T", bs.iter, sep = "")]
      
      data[, paste("R", bs.iter, sep = ".")]   <-  data[, paste("R", bs.iter, sep = ".")]  * 100
    }
  }else{   
    for(bs.iter in bs.names){
      if(!is.null(doublingvar))
        data[is.na(data[, paste("T0", bs.iter, sep = "")]), paste("T0", bs.iter, sep = "")] <-
        data[is.na(data[, paste("T0", bs.iter, sep = "")]), doublingvar]
      
      data[, paste("R", bs.iter, sep = ".")] <-
        data[, bs.iter] / data[, paste("A0T", bs.iter, sep = "")]
      
      
      data[, paste("RG.upper", bs.iter, sep = ".")] <-
        (1 + data[, paste("T0", bs.iter, sep = "")] / data[, timevar] *
           log2(data[, paste("R", bs.iter, sep = ".")])) * 100
      
      
      data[, paste("RG.lower", bs.iter, sep = ".")] <-
        (1/ data[, paste("T0", bs.iter, sep = "")] +1 / data[, timevar] *
           log2(data[, paste("R", bs.iter, sep = ".")])) 
      
      
      data[, paste("R", bs.iter, sep = ".")]   <-  data[, paste("R", bs.iter, sep = ".")]  * 100
    }
  }
  
  X$data$GI.mean <- data
  return(X)
}





findGI <- function(x, y, GI = 0.5){
  n <- length(y)
  if(min(y) <= GI){
    res <- x[1]
  }else{
    res <- Inf
  }
  
  if(max(y) >= GI)
    res <- x[length(x)]
  
  
  if(min(y) < GI & max(y) > GI ){
    
    i <- 1
    repeat {
      if(y[i + 1] < GI | i > (n-1)){
        a <- (y[i + 1] - y[i]) /(x[i + 1] - x[i])
        b <- - 1 * a * x[i] + y[i]
        res <- (GI - b) / a
        break()
      }
      i <- i + 1
    }
  }
  res
}





findOneLine <- function(x, upper, lower = NULL , at = 0){
  
  fit <- NULL
  data <- NULL
  if(min(upper, na.rm = TRUE) >= at & max(upper, na.rm = TRUE) >= at | is.null(lower)){
    
    data <- data.frame(x = x, upper = upper)  
    data <- data[order(data$x),]
    data <- data[!apply(data, 1, function(x) any(is.na(x))), ]
    fit <- isoreg(data$x, 100 - data$upper)   
    data$fitted <- 100 - fit$yf 
    data$orig <- data$upper
    data$which <- "upper" 
  }
  
  if(!(min(upper, na.rm = TRUE) >= at & max(upper, na.rm = TRUE) >= at) & !
       is.null(lower)){
    data <- data.frame(x = x, upper = upper, lower = lower)  
    data <- data[order(data$x),]
    data <- data[!apply(data, 1, function(x) any(is.na(x))), ]
    if(min(lower, na.rm = TRUE) < at & max(lower, na.rm = TRUE) < at){   
      fit <- isoreg(data$x, 100 - data$lower)    
      data$fitted <- 100 - fit$yf    
      data$orig  <- data$lower
      data$which <- "lower" 
    }
    if(min(lower, na.rm = TRUE) < at & max(lower, na.rm = TRUE) > at){ 
      fit.lower <- isoreg(data$x, 100 - data$lower)  
      fit.upper <- isoreg(data$x, 100 - data$upper)  
      
      data$fitted <- ifelse(100 - fit.upper$yf >= at, 
                            100 - fit.upper$yf,
                            100 - fit.lower$yf) 
      
      data$orig  <- ifelse(data$upper >= at, data$upper, data$lower)
      
      data$which  <- ifelse(100 - fit.upper$yf >= at, "upper", "lower")
      
      TGI <- findGI(fit.upper$x, 100 - fit.upper$yf, at)
      
      data <- rbind(data, data.frame(x = TGI, upper =at, 
                                     lower = at, fitted = at, orig = at,
                                     which = "at"))
    }
  }
  data <- data[order(data$x), c("x", "fitted", "orig", "which")]
  return(data) 
}





AUCcalc <- function (x, y) {
  m <- length(x)
  xp <- c(x, x[m:1])
  yp <- c(numeric(m), y[m:1])
  n <- 2 * m
  p1 <- sum(xp[1:(n - 1)] * yp[2:n]) + xp[n] * yp[1]
  p2 <- sum(xp[2:n] * yp[1:(n - 1)]) + xp[1] * yp[n]
  return(0.5 * (p1 - p2))
}


findAUCq <- function(x, y, q = 0, lim = NULL){
  
  auc.q <- findGI(x, y, q)
  
  x2 <- c(x[x< auc.q], auc.q)
  
  y2 <- c(y[x< auc.q] - q, 0)
  
  if(!is.null(lim) && !all(range(x) == lim)){
    new <- restictDoses(x=x2, y=y2, lim)
    x2 <- new$x
    y2 <- new$y
  }
  
  AUCcalc(x2, y2)
}

restictDoses <- function(x, y, lim = NULL){
  x.min <- min(lim, na.rm = TRUE)
  x.max <- max(lim, na.rm = TRUE)
  
  if(min(x, na.rm = TRUE) == x.min)
    y.min <- y[x == x.min]
  
  if(min(x, na.rm = TRUE) < x.min)
    if(any(x.min == x)){
      y.min <- y[x == x.min]
    }else{
      x1 <- max(x[x< x.min])
      x2 <- min(x[x > x.min])
      y1 <- y[x == x1]
      y2 <- y[x == x2]
      y.min <- (y2 - y1) / (x2 - x1) * (x.min - x1) + y1
      
    }
  
  if(max(x, na.rm = TRUE) == x.max)
    y.max <- y[x == x.max]
  
  if(max(x, na.rm = TRUE) < x.max)
    y.max <- y[x == max(x, na.rm = TRUE)]
  
  if(max(x, na.rm = TRUE) > x.max)
    if(any(x.max == x)){
      y.max <- y[x == x.max]
    }else{
      x1 <- max(x[x< x.max])
      x2 <- min(x[x > x.max])
      y1 <- y[x == x1]
      y2 <- y[x == x2]
      y.max <- (y2 - y1) / (x2 - x1) * (x.max - x1) + y1
    }
  
  x.new <- c(x.min, x[x > x.min & x < x.max], x.max)
  y.new <- c(y.min, y[x > x.min & x < x.max], y.max) 
  
  return(list(x = x.new, y = y.new))
}


BGfunction <- function(A, B = NULL, type = NULL,
                       plate = NULL,
                       plateset = NULL,
                       data = NULL, ...,
                       temp = NULL,
                       weights = NULL, fitted.a = FALSE,
                       platesets = platesets,
                       contr = c("sum", "helmert", "treatment"),
                       outlier.test = FALSE,
                       outlier.iter = 4,
                       varpower.min = 10e-4,
                       varpower.iter = 50,
                       verbose = verbose,
                       parametrisation = c("none", "exp", "expinv", "square", "abs")){
  
  temp$plate.help <- as.factor(as.character(temp[, plate]))
  temp[, plate] <- as.factor(as.numeric(as.factor(as.character(temp[, plate]))))
  
  n.rep <- length(levels(temp[, plate]))
  contr.orig <- contr
  fit.back2 <- NULL
  Btype <- temp[, type][temp[, B]][1]
  if(Btype %in% temp[, type][temp[, B]== 0]) Btype <- "sdfkj445"
  ## Start values
  temp$NBcut <- ifelse(temp$NB < 0.01, 0.01, temp$NB)
  
  st.ab <- switch(parametrisation,
                  none = aggregate(as.formula(paste("NB ~ ", type)) ,
                                   data = temp, FUN = mean),
                  abs = aggregate(as.formula(paste("(NBcut) ~ ", type)) ,
                                  data = temp, FUN = mean),
                  square = aggregate(as.formula(paste("sqrt(NBcut) ~ ", type)) ,
                                     data = temp, FUN = mean),
                  exp =aggregate(as.formula(paste("log(NBcut) ~ ", type)) ,
                                 data = temp, FUN = mean),
                  expinv =aggregate(as.formula(paste("1/log(NBcut) ~ ", type)) ,
                                    data = temp, FUN = mean))
  
  st.ab <-st.ab[st.ab[,type] != Btype,]
  st.a <- st.ab[, 2]
  
  st.b <- aggregate(as.formula(paste(A, "~", plate)),
                    data = temp[temp[, B] == 1,], FUN = mean)[,2]
  st.b[is.na(st.b)] <- 0.4
  ## contrasts
  temp$content2 <- factor(temp[, type])
  levels(temp$content2)[levels(temp$content2) == Btype] <-
    levels(temp$content2)[levels(temp$content2) != Btype][1]
  
  if(n.rep > 1){
    contr <- switch(contr.orig[1],
                    sum       = contr.sum(levels(temp[,plate])),
                    treatment = contr.treatment(levels(temp[,plate])),
                    helmert   = contr.helmert(levels(temp[,plate])))
    
    if(contr.orig[1] == "helmert"|contr.orig[1] == "sum")
      colnames(contr) <- rownames(contr)[-1]
    colnames(contr) <- paste("plate", colnames(contr), sep = "")
    mat <- ((contr[temp[,plate], ]))
    if(!any(class(mat) == "matrix")){
      mat <- as.matrix(mat)
      colnames(mat) <- paste("plate", rownames(contr)[-1], sep = "")
    }
    temp <- cbind(temp, mat)
    
    ## Dummy variables
    
    dummys <- paste("d", levels(temp[,plate])[-1], "*",
                    colnames(contr),
                    sep = "", collapse = "+")
  }
  ## Setup for gnls
  temp$B2 <- ifelse(temp[, B], 0, 1)
  
  if(n.rep > 1){
    formula <- switch(parametrisation,
                      none = as.formula(paste(A, " ~ B2 * a * exp(",
                                              dummys, ") + b", sep ="")),
                      abs = as.formula(paste(A, " ~ B2 * abs(a) *exp(",
                                             dummys, ") + b", sep ="")),
                      square = as.formula(paste(A, " ~ B2*(a)^2 *exp(",
                                                dummys, ") + b", sep ="")),
                      exp = as.formula(paste(A, " ~ B2*exp(a) *exp(",
                                             dummys, ") + b", sep ="")),
                      expinv = as.formula(paste(A, " ~ B2*exp(1/a) *exp(",
                                                dummys, ") + b", sep ="")))
    
    params <- list(as.formula(a ~ content2 -1),
                   as.formula(paste("b ~", plate)))
    
    d.par <- paste("d", levels(temp[, plate])[-1] ,
                   " ~ 1", sep = "")
    
    for(i in 1:length(d.par))
      params[[length(params) + 1]] <- as.formula(d.par[i])
  }else{
    formula <- switch(parametrisation,
                      none = as.formula(paste(A, " ~ B2 * a + b", sep ="")),
                      abs = as.formula(paste(A, " ~ B2 * abs(a) + b", sep ="")),
                      square = as.formula(paste(A, " ~ B2*(a)^2 + b", sep ="")),
                      exp = as.formula(paste(A, " ~ B2*exp(a) + b", sep ="")),
                      expinv = as.formula(paste(A, " ~ B2*exp(1/a) + b", sep ="")))
    params <- list(as.formula(a ~ content2 -1),
                   as.formula(paste("b ~", 1)))
  }
  
  
  start <- list(a = st.a,
                b = st.b)
  
  if(n.rep > 1){
    temp.d <- temp[temp[, type] != Btype,
                   eval(parse(text = paste('c("NB", "NBcut", "absorbance", "plate","', type, '",',
                                           paste('"', 'plate', levels(temp[, plate])[-1], '"',
                                                 sep = '', collapse = ','), ')', sep = "" )))]
    
    d.y <- aggregate(NB ~ plate, FUN = mean, data=temp.d)[,2]
    # aggregate(NB ~  new.type, FUN = mean, data=temp.d)
    
    coef <- coef(lm((d.y/mean(d.y)) ~ contr))
    
    for(i in paste("d", levels(temp[,plate])[-1], sep = ""))
      start[[i]] <- 0
    
    if(length(coef) > 2){
      for(i in levels(temp[,plate])[-1])
        start[[paste("d", i, sep= "")]] <- coef[paste("contrplate",i, sep = "")]
    }else{
      for(i in levels(temp[,plate])[-1])
        start[[paste("d", i, sep= "")]] <- coef[2]
    }
  }
  fit.back <- fit.back2
  fit.back2 <- "error"
  ## fit using gnls from the nlme package
  # temp <- temp[sample(1:nrow(temp)),]
  try(fit.back2 <-
        nlme::gnls(formula,
                   params = params,
                   data = temp,
                   start =start,
                   control =gnlsControl(
                     maxIter = 1000,
                     returnObject = TRUE)), silent = TRUE)
  
  fit <- NULL
  diff2 <- Inf
  temp$outlier <- 0
  if(class(fit.back2)[1] != "character") {
    temp$fitted <- fitted(fit.back2)
    fit <- "error"
    if(!is.null(weights)) {
      try(fit <-
            update(fit.back2, weights = varPower(form =  "fitted")))
      if(class(fit)[1] != "character"){
        fit.back2 <- fit
        power <- fit.back2$modelStruct$varStruct[1]
        
        temp$outlier <- 0
        
        diff2 <- Inf
        #varpower.min <- 10e-4
        
        iter <- 0
        temp$fitted <- NA
        
        count.stop <- 0
        
        while(diff2[length(diff2)] > varpower.min & iter <= varpower.iter & count.stop < 2){
          iter <- iter + 1
          temp$fitted <- fitted(fit.back2)# + 0.15
          if(fitted.a){
            if(n.rep > 1){
              for(k in 1:length(levels(temp$plate))){
                if(k == 1){
                  temp$fitted[temp$plate == levels(temp$plate)[1]] <-
                    temp$fitted[temp$plate == levels(temp$plate)[1]] -
                    (coef(fit.back2)["b.(Intercept)"]/2)
                }else{
                  temp$fitted[temp$plate == levels(temp$plate)[k]] <-
                    temp$fitted[temp$plate == levels(temp$plate)[k]] -
                    (coef(fit.back2)[paste("b.plate",
                                           levels(temp$plate)[k], sep = "")]/2)
                }
              }
            }else{
              temp$fitted[temp$plate == levels(temp$plate)[1]] <-
                temp$fitted[temp$plate == levels(temp$plate)[1]] -
                (coef(fit.back2)["b"]/2)
            }
          }
          
          try(fit.back2 <-
                update(fit.back2, weights = varPower(form =  "fitted")))
          
          power[iter + 1] <- fit.back2$modelStruct$varStruct[1]
          diff2[iter] <- abs(power[iter+1]-power[iter])
 
          if(any(duplicated(power))){
            number <- which.min(abs(power[2:length(power)])) - 1
            power <- power[2:length(power)][number]
            count.stop <- count.stop + 1
            if(verbose & count.stop ==2)            
              cat("the optimisation routine is running in circles\n")
          }
        }
        #temp$outlier <- 0
      }
    }
    
    if(outlier.test != FALSE) {
      outlier <- 0
      out.iter <-  0
      temp$outlier <- 3
      
      while(!all(outlier == temp$outlier) & out.iter < outlier.iter){
        if(all(temp$outlier ==3)) temp$outlier <- 0
        out.iter <- out.iter + 1
        
        outlier <- temp$outlier
        power <- fit.back2$modelStruct$varStruct[1]
        res <- resid(fit.back2)
        res <- res/attributes(res)$std
        temp[temp$outlier == 0, ]$outlier <- ifelse(abs(res) > outlier.test, 1, 0)
        
        #if(!is.null(weights)){
        #  #res <- scale(resid(fit.back2)/fitted(fit.back2)^power, center = FALSE)
        #  res <- res/attributes(res)$std
        #  temp[temp$outlier == 0, ]$outlier <- ifelse(abs(res) > outlier.test, 1, 0)
        #}else{
        #  res <- res/attributes(res)$std
        #  temp[temp$outlier == 0, ]$outlier <- ifelse(abs(res) > outlier.test, 1, 0)
        #}
        for(plate.iter in levels(temp[, plate]))
          if(sum(temp[outlier == 0, ]$outlier[temp[outlier == 0, type] ==
                                                Btype & temp[outlier == 0, plate] == plate.iter]== 0) < 2)
            temp[outlier == 0, ]$outlier[temp[outlier == 0, type] ==
                                           Btype & temp[outlier == 0, plate] == plate.iter] <- 0
        
        for (plate.iter in levels(temp[, plate]))
          for(l in unique(temp[, type]))
            if(sum(temp[outlier == 0, ]$outlier[temp[outlier == 0, type] ==
                                                  l] == 0) < ceiling(length(
                                                    temp$outlier[temp[, type] ==
                                                                   l] )/2))
              temp[outlier == 0, ]$outlier[temp[outlier == 0, type] ==
                                             l & temp[outlier == 0, plate] == plate.iter] <- 0
        
        fit.back3 <- NULL
        if(!is.null(weights)){
          try(fit.back3 <-
                nlme::gnls(formula,
                           weights = varPower(form = as.formula(paste(" ~", "fitted"))),
                           #...,
                           params = params,
                           data = temp[temp$outlier == 0,],
                           start =start,
                           control =gnlsControl(
                             maxIter = 1000,
                             returnObject = TRUE)), silent = TRUE)
        }else{
          try(fit.back3 <-
                nlme::gnls(formula,
                           params = params,
                           data = temp[temp$outlier == 0,],
                           start =start,
                           control =gnlsControl(
                             maxIter = 1000,
                             returnObject = TRUE)), silent = TRUE)
        }
        if(!is.null(fit.back3)){
          fit.back2 <- fit.back3
          if(!is.null(weights)){
            power <- fit.back2$modelStruct$varStruct[1]
            diff2 <- Inf
            iter <- 0
            temp$fitted <- NA
            
            count.stop <- 0
            
            
            while(diff2[length(diff2)] > varpower.min & iter <= varpower.iter & count.stop < 2) {
              iter <- iter+1
              temp[temp$outlier == 0,]$fitted <- fitted(fit.back2)# + 0.15
              if(fitted.a){
                for(plate.iter in 1:length(levels(temp$plate))){
                  if(n.rep > 1) {
                    if(plate.iter == 1){
                      temp$fitted[temp$plate == levels(temp$plate)[1]] <-
                        temp$fitted[temp$plate == levels(temp$plate)[1]] -
                        (coef(fit.back2)["b.(Intercept)"] / 2)
                    }else{
                      temp$fitted[temp$plate == levels(temp$plate)[plate.iter]] <-
                        temp$fitted[temp$plate == levels(temp$plate)[plate.iter]] -
                        (coef(fit.back2)[paste("b.plate",
                                               levels(temp$plate)[plate.iter], sep = "")] / 2)
                    }
                  }else{
                    temp$fitted[temp$plate == levels(temp$plate)[1]] <-
                      temp$fitted[temp$plate == levels(temp$plate)[1]] -
                      (coef(fit.back2)["b"] / 2)
                  }
                }
              }
              try(fit.back2 <-
                    update(fit.back2, weights = varPower(form =  "fitted")))
              
              power[iter + 1] <- fit.back2$modelStruct$varStruct[1]
              diff2[iter] <- abs(power[iter+1]-power[iter])
   
              if(any(duplicated(power))){
                number <- which.min(abs(power[2:length(power)])) - 1
                power <- power[2:length(power)][number]
                count.stop <- count.stop + 1
                if(verbose & count.stop ==2)            
                  cat("the optimisation routine is running in circles\n")
              }
            }
          }
        }else{
          diff2 <- Inf
          iter <- 0
          temp$fitted <- NA
          
          count.stop <- 0
          
          temp$outlier <- outlier
          if(!is.null(weights)){
            power[iter + 1] <- fit.back2$modelStruct$varStruct[1]
            diff2[iter] <- abs(power[iter+1]-power[iter])
          }
        }
        if(verbose)
          cat(out.iter,  ": number of outliers = ", sum(temp$outlier), "\n")
      }
    }
  }
  
  col <- as.numeric(temp[temp$outlier == 0, plate])
  pch <- ifelse(temp[temp$outlier == 0, type] == Btype, 2, 1)
  if(verbose)
    if(class(power) == "numeric" & class(diff2) == "numeric"){
      cat(platesets, "power = ", round(power, 3) , "\n")
      cat(platesets, "power = ", round(diff2, 3) , "\n")
    }
  
  if(n.rep > 2){
    COEF <- coef(fit.back2)
    bs   <- COEF["b.(Intercept)"] + c(0, COEF[grepl("b.plate", names(COEF))])
    names(bs) <- levels(temp$plate)
    d <- exp(contr[temp$plate, ] %*% COEF[paste("d", levels(temp$plate)[-1], sep ="")])
    
    temp$A <- (temp$absorbance - bs[temp$plate]) / d
  }
  if(n.rep == 2){
    COEF <- coef(fit.back2)
    bs   <- COEF["b.(Intercept)"] + c(0, COEF[grepl("b.plate", names(COEF))])
    names(bs) <- levels(temp$plate)
    d <- exp(contr[temp$plate, ] * COEF[paste("d", levels(temp$plate)[-1], sep ="")])
    
    temp$A <- (temp$absorbance - bs[temp$plate]) / d
  }
  if(n.rep == 1)
    temp$A <- temp$absorbance - coef(fit.back2)["b"]
  
  
  temp[, plate] <- temp$plate.help
  return(list(temp = temp, fit = fit.back2, fit2 = fit,
              col = col, pch = pch, start = start, Btype = Btype, contr = contr))
}



growth.function <- function(data, cut = 0.025,
                            parametrisation = "reci",
                            type            = "new.type",
                            timevar         = "time",
                            dosevar         = "Concentration",
                            incubationvar   = "incubation",
                            additivevar     = "Additive",
                            controlval      = "Control",
                            A               = "BC",
                            n.start.values  = 10,
                            fit.reci.first  = TRUE,
                            ...){
  
  fit.list <- NULL
  fit <- NULL
  fit2 <- NULL
  if(parametrisation == "reci")
    fit.reci.first = FALSE
  
  temp <- data.frame(Bcut = ifelse(data[, A] < cut, cut, data[, A]),
                     control = ifelse(data[, additivevar] == controlval, 0, 1),
                     type = data[, type],
                     protocol.type = data[, "type"],
                     t = data[,timevar] + data[,incubationvar] / 2,
                     Concentration = data[, dosevar])
  
  tc.mat <- data.frame(matrix(NaN, ncol = 6, nrow = length(unique(temp$type))))
  
  rownames(tc.mat) <- unique(temp$type)
  
  colnames(tc.mat) <- c(dosevar, "type", "st.N0", "st.T0", "st.TC", "st.TCp")
  
  conc.mat <- aggregate(Concentration ~ 
                          type + protocol.type, FUN = mean, data = temp)
  tc.mat[as.character(conc.mat[, "type"]), dosevar] <- conc.mat[,"Concentration"]
  tc.mat[as.character(conc.mat[, "type"]), "type"]  <- conc.mat[, "protocol.type"]
  for(content in unique(temp$type)) {
    
    fit.Bcut <- lm(log(temp[temp$type == content, "Bcut"]) ~
                     temp[temp$type == content, "t"])
    
    tc.mat[content, "st.TC"] <- 1 / coef(fit.Bcut)[2]
    
  }
  
  tc.mat[, "st.T0"] <- tc.mat[paste(0, controlval), "st.TC"]
  
  tc.mat[, "st.TCp"] <- 1/(1/tc.mat[, "st.TC"]-1/tc.mat[, "st.T0"])
  
  tc.mat[tc.mat[, "st.TCp"] == Inf, "st.TCp"] <- -1000
  
  tc.mat[, "st.N0"] <- mean(temp[temp$t == min(temp$t), "Bcut"])
  
  temp$type <- as.character(temp$type)
  temp$type[temp$type == paste(0, controlval)] <-
    c(unique(temp$type) %w/o% paste(0, controlval))[1]
  temp$type <- as.factor(temp$type)
  
  ## remove the ones which have been deemed outliers
  #tc.mat <- as.data.frame(tc.mat)
  temp <- temp[temp$type %in% rownames(tc.mat[!is.na(tc.mat[, "st.TCp"]),]),]
  temp$type <- as.factor(as.character(temp$type)) 
  
  fit <- list()
  #class(fit) <- "error"
  #fit[[1]] <- "error"
  #fit$logLik <- 0
  if(fit.reci.first){
    parametrisation2 <- "reci"
    
    formula <- switch(parametrisation2,
                      none =   formula(Bcut ~ N0 * 2^((1/T0 - control/Tcp) * t)),
                      reci =   formula(Bcut ~ N0 * 2^((T0 - control*Tcp) * t)),
                      abs  =   formula(Bcut ~ N0 * 2^((T0 - control*abs(Tcp)) * t)),
                      square = formula(Bcut ~ N0 * 2^((T0 - control*(Tcp)^2) * t)),
                      exp  =   formula(Bcut ~ N0 * 2^((T0 - control*exp(Tcp)) * t)))  
    
    params <- list(formula(Tcp ~ type -1 ),
                   formula(T0  ~ 1),
                   formula(N0  ~ 1))
    
    stp <- abs(ifelse(- tc.mat[levels(temp$type), "st.TCp"] < 0,
                      - 200000, tc.mat[levels(temp$type), "st.TCp"]))
    
    start.TCp <- switch(parametrisation2,
                        none = - tc.mat[levels(temp$type), "st.TCp"],
                        reci = - 1/tc.mat[levels(temp$type), "st.TCp"],
                        abs  =  1 / stp ,
                        square =  sqrt(1/(stp)),
                        exp = log(1/stp))
    
    start.T0 <- switch(parametrisation2,
                       none = tc.mat[levels(temp$type), "st.T0"][1],
                       reci = 1/tc.mat[levels(temp$type), "st.T0"][1],
                       abs  = abs(1/tc.mat[levels(temp$type), "st.T0"])[1],
                       square = abs(1/tc.mat[levels(temp$type), "st.T0"])[1],
                       exp = abs(1/tc.mat[levels(temp$type), "st.T0"])[1])
    
    start <- list(Tcp = start.TCp,
                  T0  = start.T0,
                  N0  = mean(temp[temp$t == min(temp$t), "Bcut"]))
    
    fit2 <- NULL
    
    try(fit2 <-
          nlme::gnls(formula,# ...,
                     params = params,
                     data = temp,
                     start =start,
                     control =gnlsControl(
                       maxIter = 1000,
                       returnObject = TRUE)), silent = TRUE)
    
    if(!is.null(fit2) & parametrisation != "exp"){
      stp <- coef(fit2)[1:(length(coef(fit2))-2)]
      
      start.TCp <- switch(parametrisation,
                          none = - 1/stp,
                          reci = - stp,
                          abs  =  abs(stp) ,
                          square = sqrt(abs(stp)),
                          exp =    log(abs(stp)),
                          expinv = log(1/abs(stp)))
      start <- list(Tcp = start.TCp,
                    T0  = coef(fit2)["T0"],
                    N0  = coef(fit2)["N0"])
      
      formula <- switch(parametrisation,
                        none =   formula(Bcut ~ N0 * 2^((1/T0 - control/Tcp) * t)),
                        reci =   formula(Bcut ~ N0 * 2^((T0 - control*Tcp) * t)),
                        abs  =   formula(Bcut ~ N0 * 2^((T0 - control*abs(Tcp)) * t)),
                        square = formula(Bcut ~ N0 * 2^((T0 - control*(Tcp)^2) * t)),
                        exp =    formula(Bcut ~ N0 * 2^((T0 - control*exp(Tcp)) * t)),
                        expinv = formula(Bcut ~ N0 * 2^((T0 - control/exp(Tcp)) * t)))  
      
      try(fit <-
            nlme::gnls(formula,# ...,
                       params = params,
                       data = temp,
                       start =start,
                       control =gnlsControl(
                         maxIter = 1000,
                         returnObject = TRUE)), silent = TRUE)
      
      fit.list <- list()
      
      #       start.TCp <- coef(fit)[1:(length(coef(fit))-2)]
      #       fit.list[[1]] <- fit 
      #       number <- 1
      #       for(i in 2:n.start.values){
      #         if(parametrisation == "abs")
      #         start <- list(Tcp = start.TCp + rnorm(length(start.TCp), 0, start.TCp/30),
      #                       T0  = coef(fit2)["T0"] + rnorm(1, 0, coef(fit2)["T0"]/30),
      #                       N0  = coef(fit2)["N0"] + rnorm(1, 0, coef(fit2)["N0"]/30))
      #          if(parametrisation == "square")
      #             start <- list(Tcp = start.TCp + rnorm(length(start.TCp), 0, start.TCp/30),
      #                       T0  = coef(fit2)["T0"] + rnorm(1, 0, coef(fit2)["T0"]/30),
      #                       N0  = coef(fit2)["N0"] + rnorm(1, 0, coef(fit2)["N0"]/30))
      #         fit.l <- NULL
      #         try(fit.l <-
      #               nlme::gnls(formula,# ...,
      #                          params = params,
      #                          data = temp,
      #                          start =start,
      #                          control =gnlsControl(
      #                            maxIter = 1000,
      #                            returnObject = TRUE)), silent = TRUE)
      #          
      #         
      #         
      #         if(!is.null(fit.l)){ 
      #           number <- number + 1
      #           fit.list[[number]] <- fit.l
      #       }}
      #       
      #       fit <- fit.list[[which.max((unlist(lapply(fit.list, function(x) x$logLik))))]]
    }
    
    
    if(!is.null(fit2) & parametrisation == "exp"){
      parametrisation3 <- "abs" 
      stp <- coef(fit2)[1:(length(coef(fit2))-2)]
      
      start.TCp <- switch(parametrisation3,
                          none = - 1/stp,
                          reci = - stp,
                          abs  =  abs(stp) ,
                          square =  sqrt(abs(stp)),
                          exp = log(abs(stp)))
      start <- list(Tcp = start.TCp,
                    T0  = coef(fit2)["T0"],
                    N0  = coef(fit2)["N0"])
      
      formula <- 
        switch(parametrisation3,
               none =   formula(Bcut ~ N0 * 2^((1/T0 - control/Tcp) * t)),
               reci =   formula(Bcut ~ N0 * 2^((T0 - control*Tcp) * t)),
               abs  =   formula(Bcut ~ N0 * 2^((T0 - control*abs(Tcp)) * t)),
               square = formula(Bcut ~ N0 * 2^((T0 - control*(Tcp)^2) * t)),
               exp =    formula(Bcut ~ N0 * 2^((T0 - control*exp(Tcp)) * t)))  
      
      try(fit2 <-
            nlme::gnls(formula,# ...,
                       params = params,
                       data = temp,
                       start =start,
                       control =gnlsControl(
                         maxIter = 1000,
                         returnObject = TRUE)), silent = TRUE)
      
      stp <- coef(fit2)[1:(length(coef(fit2))-2)]
      
      start.TCp <- switch(parametrisation3,
                          none = log(abs(1/stp)),
                          reci = log( abs(stp)),
                          abs  = log( abs(stp)) ,
                          square =  log(stp^2),
                          exp = stp)
      
      start <- list(Tcp = start.TCp,
                    T0  = coef(fit2)["T0"],
                    N0  = coef(fit2)["N0"])
      
      formula <- switch(parametrisation,
                        none =   formula(Bcut ~ N0 * 2^((1/T0 - control/Tcp) * t)),
                        reci =   formula(Bcut ~ N0 * 2^((T0 - control*Tcp) * t)),
                        abs  =   formula(Bcut ~ N0 * 2^((T0 - control*abs(Tcp)) * t)),
                        square = formula(Bcut ~ N0 * 2^((T0 - control*(Tcp)^2) * t)),
                        exp =    formula(Bcut ~ N0 * 2^((T0 - control*exp(Tcp)) * t)))  
      
      try(fit <-
            nlme::gnls(formula,# ...,
                       params = params,
                       data = temp,
                       start =start,
                       control =gnlsControl(
                         maxIter = 1000,
                         returnObject = TRUE)), silent = TRUE)
    }
    
  }
  #if(class(fit)[1] == "error" | !fit.reci.first){
  if(is.null(fit) | !fit.reci.first){
    formula <- switch(parametrisation,
                      none =   formula(Bcut ~ N0 * 2^((1/T0 - control/Tcp) * t)),
                      reci =   formula(Bcut ~ N0 * 2^((T0 - control*Tcp) * t)),
                      abs  =   formula(Bcut ~ N0 * 2^((T0 - control*abs(Tcp)) * t)),
                      square = formula(Bcut ~ N0 * 2^((T0 - control*(Tcp)^2) * t)),
                      exp =    formula(Bcut ~ N0 * 2^((T0 - control*exp(Tcp)) * t)))  
    
    params <- list(formula(Tcp ~ type -1 ),
                   formula(T0  ~ 1),
                   formula(N0  ~ 1))
    
    stp <- abs(ifelse(- tc.mat[levels(temp$type), "st.TCp"] < 0,
                      - 200000, tc.mat[levels(temp$type), "st.TCp"]))
    
    start.TCp <- switch(parametrisation,
                        none   = - tc.mat[levels(temp$type), "st.TCp"],
                        reci   = - 1/tc.mat[levels(temp$type), "st.TCp"],
                        abs    =  1 / stp ,
                        square =  sqrt(1/(stp)),
                        exp    = log(1/(stp)))
    
    start.T0 <- switch(parametrisation,
                       none = tc.mat[levels(temp$type), "st.T0"][1],
                       reci = 1/tc.mat[levels(temp$type), "st.T0"][1],
                       abs  = abs(1/tc.mat[levels(temp$type), "st.T0"])[1],
                       square = abs(1/tc.mat[levels(temp$type), "st.T0"])[1],
                       exp = abs(1/tc.mat[levels(temp$type), "st.T0"])[1])
    
    start <- list(Tcp = start.TCp,
                  T0  = start.T0,
                  N0  = mean(temp[temp$t == min(temp$t), "Bcut"]))
    
    try(fit <-
          nlme::gnls(formula,# ...,
                     params = params,
                     data = temp,
                     start =start,
                     #verbose =TRUE,
                     control =gnlsControl(
                       maxIter = 1000,
                       returnObject = TRUE)), silent = TRUE)
  }
  
  
  
  #if(class(fit)[1] == "error"){
  if(is.null(fit)){
    
    start.TCp <- switch(parametrisation,
                        none = - tc.mat[levels(temp$type), "st.TCp"],
                        reci = - 1/tc.mat[levels(temp$type), "st.TCp"],
                        abs  =  abs(1/tc.mat[levels(temp$type), "st.TCp"]) ,
                        square = sqrt(1/abs(tc.mat[levels(temp$type), "st.TCp"])),
                        exp = log(1/abs(tc.mat[levels(temp$type), "st.TCp"])))
    start <- list(Tcp = start.TCp,
                  T0  = start.T0,
                  N0  = mean(temp[temp$t == min(temp$t), "Bcut"]))
    try(fit <-
          nlme::gnls(formula,# ...,
                     params = params,
                     data = temp,
                     start =start,
                     #verbose =TRUE,
                     control =gnlsControl(
                       maxIter = 1000,
                       returnObject = TRUE)), silent = TRUE)
    
  }
  #id.fits[[id]] <- fit
  
  #if(class(fit)[1] != "error"){
  if(!is.null(fit)){
    COEF <- coef(fit)
    
    tc.mat[, "N0"] <- COEF["N0"]
    
    tc.mat[, "T0"] <-
      switch(parametrisation,
             none = COEF["T0"],
             reci = 1/COEF["T0"],
             abs  = 1/COEF["T0"],
             square = 1/COEF["T0"],
             exp = 1/COEF["T0"])
    
    
    Tcp.coef <- COEF[grepl("Tcp.type", names(coef(fit)))]
    names(Tcp.coef) <- gsub("Tcp.type", "", names(Tcp.coef) )
    tc.mat[names(Tcp.coef), "TCp"] <-
      switch(parametrisation,
             none =  - Tcp.coef,
             reci =  - 1/Tcp.coef,
             abs  =  - 1/abs(Tcp.coef),
             square = - 1/(Tcp.coef)^2,
             exp = -1/exp(Tcp.coef))
    
    
    tc.mat[,"TC"] <- 1/(1/tc.mat[,"T0"] + 1/tc.mat[,"TCp"])
    tc.mat[,"G"] <- (tc.mat[,"T0"] / tc.mat[,"TC"] ) *100
    
  }else{
	fit <- list()
	class(fit) <- "error"
	fit[[1]] <- "error"
	fit$logLik <- 0
  }
  return(list(fit = fit, summary = tc.mat, fit.list = fit.list))
}  



