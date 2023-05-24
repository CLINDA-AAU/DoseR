##
##
## Convert between units
##
##

concConvert <- function(x, mol.mass = 1, from = "ug/ml", to = "mol/l",
                        logfun.from = "",
                        logfun.to   = ""){
  
  logfun.from <- ifelse(logfun.from == "", "nolog", logfun.from)
  logfun.to   <- ifelse(logfun.to == "",   "nolog", logfun.to)
  
  if(length(mol.mass) == 1)
    mol.mass <- rep(mol.mass, length(x))
  
  if(length(to) == 1)
    to <- rep(to, length(x))
  
  if(length( logfun.from) == 1)
    logfun.from <- rep( logfun.from, length(x))
  
  if(length(logfun.to) == 1)
    logfun.to <- rep(logfun.to, length(x))
  
  
  multiples <- data.frame(
    name = c("", "deca", "hecto", "kilo", "mega", "giga", 
             "tera", "peta", "exa", "zetta", "yotta"),
    
    prefix = c("", "da", "h", "k",  "M",  "G",	"T",	"P",	"E",	"Z",	"Y"),
    factor = 10^c(0, 1, 2, 3, 6, 9, 12, 15, 18, 21, 24))
  
  fractions <- data.frame(
    name = c("deci",  "centi",	"milli",	"micro",	"nano",	
             "pico",	"femto",	"atto",	"zepto",	"yocto"),
    prefix = c("d",  "c",	"m",	"u",	"n",	"p",	"f",	"a",	"z",	"y"),
    factor = 10^(-c(1, 2, 3, 6, 9, 12, 15, 18, 21, 24)))
  
  prefixes <- rbind(multiples, fractions )
  
  rownames(prefixes) <- prefixes$prefix
  
  
  x <- switch(logfun.from[1],
              nolog    = x,
              log10 = 10^x,
              log2  =  2^x,
              log   = exp(x))
  
  
  
  from.1 <- #sapply(strsplit(from, "/"), function(x) x[1])
    unlist(lapply(strsplit(from, "/"), function(x) x[[1]]))
  
  from.1 <- ifelse(grepl("mol", from.1), gsub("mol", "", from.1),
                   substr(x=from.1,1, nchar(from.1)-1))
  
  from.2 <- unlist(lapply(strsplit(from, "/"), function(x) x[[2]]))
  #sapply(strsplit(from, "/"), function(x) x[2])
  from.2 <- substr(x=from.2,1, nchar(from.2)-1)
  
  
  to.1 <- #sapply(strsplit(to, "/"), function(x) x[1])
    unlist(lapply(strsplit(to, "/"), function(x) x[[1]]))
  
  to.1 <- ifelse(grepl("mol", to.1), gsub("mol", "", to.1), 
                 substr(x=to.1, 1, nchar(to.1)-1))
  
  to.2 <-  unlist(lapply(strsplit(to, "/"), function(x) x[[2]]))
  to.2 <- substr(x=to.2,1, nchar(to.2)-1)
  
  mol.mass.2 <- ifelse(grepl("mol", to),  1 / mol.mass, mol.mass)
  
  mol.mass.2 <- ifelse(grepl("mol", from), mol.mass, mol.mass.2)
  
  mol.mass.2 <- ifelse(grepl("mol", from) & grepl("mol", to) | 
                         !grepl("mol", from) & !grepl("mol", to), 1, mol.mass.2 )
  
  
  #which(to.1[1] == row.names(prefixes))
  
  from.wh.1 <- sapply(from.1, function(x) which(x == row.names(prefixes)))
  from.wh.2 <- sapply(from.2, function(x) which(x == row.names(prefixes)))
  
  to.wh.1 <- sapply(to.1, function(x) which(x == row.names(prefixes)))
  to.wh.2 <- sapply(to.2, function(x) which(x == row.names(prefixes)))
  
  
  conc <- x * mol.mass.2 / prefixes[from.wh.2, 3]  * prefixes[from.wh.1, 3] *
    prefixes[to.wh.2, 3]  / prefixes[to.wh.1, 3]
  if(!(logfun.to[1] == "nolog"))
    conc <- do.call(logfun.to[1], list(x = conc))
  conc
}




##
##
## Create metadata
##
##

# dbf.path       = file.path(SeedingConc.ext.dir, "DBF files")
# protocol.path  = file.path(SeedingConc.ext.dir, "Protocols")
# correctionname = "Correction"
# date.format    = "%d%m%y"
# data.file      = "GeneratedData/SeedingConc"

createMetaData <- function( 
  data             = NULL,
  data.file           = file.path(getwd(), "Absorbance"),
  save             = TRUE,
  namevar          = "Cellline",
  drugvar          = "chemo",
  protocolvar      = "R.protocol",
  identifier       = "identifier",
  timevar          = "Hour",
  correctionname   = "Control",
  incubationvar    = "incubation",
  doublingvar       = NULL,
  format           = "long",
  dbf.path         = getwd(),
  protocol.path    = getwd(),
  dbf.files        = NULL,
  file.extension = ".dbf",
  protocol.files   = NULL,
  are.paths.full    = TRUE,
  colnames         = c("namevar", "drugvar", "protocolvar", 
                       "identifier", "timevar"),
  sep              = c(";"),
  namesep          = " ",
  identifiersep    = "_",
  date.format      = NULL,
  unit             = "ug/ml",
  additional.metadata = NULL,
  show.warnings    = TRUE,
  update           = TRUE,
  idcomb           = NULL,#c(namevar, drugvar, identifier) # vars used to generate id variable
  idvar            = "sampleid", 
  namecols         = NULL,
  dbf.file.var     = "dbf.file.name",
  shiny.input      = NULL)
{
  
  
  #############################################
  ##
  ##     Function starts 
  ##
  #############################################
  
  call2 <- match.call()
  
  this.call <- list()
  
  myfor <- formals(createMetaData)               ## formals with default arguments
  for (v in names(myfor)){
    if (!(v %in% names(call2))){
      this.call[[v]] <- myfor[[v]]
    }else{
      this.call[[v]] <- call2[[v]]
    }  ## if arg is missing I add it
  }
  
  this.call$call <- call2  
  
  ##############################################
  ##
  ##   Read from directory
  ##
  ##############################################
  if(is.null(data)) {  
    A.data <- 
      do.call(createMetaDataDir, list(
        dbf.files        = eval(this.call$dbf.files),
        file.extension   = eval(this.call$file.extension),
        protocol.files   = eval(this.call$protocol.files),
        dbf.path         = eval(this.call$dbf.path),
        protocol.path    = eval(this.call$protocol.path),
        sep              = eval(this.call$sep),
        colnames         = eval(this.call$colnames),
        namevar          = eval(this.call$namevar),
        drugvar          = eval(this.call$drugvar),
        protocolvar      = eval(this.call$protocolvar),
        identifier       = eval(this.call$identifier),
        timevar          = eval(this.call$timevar),
        incubationvar    = eval(this.call$incubationvar),
        doublingvar       = eval(this.call$doublingvar),
        namesep          = eval(this.call$namesep),
        identifiersep    = eval(this.call$identifiersep),
        namecols         = eval(this.call$namecols),
        date.format      = eval(this.call$date.format),
        correctionname   = eval(this.call$correctionname),
        unit             = eval(this.call$unit),
        additional.metadata = eval(this.call$additional.metadata),
        idcomb           = eval(this.call$idcomb),
        idvar            = eval(this.call$idvar), 
        dbf.file.var     = eval(this.call$dbf.file.var),
        show.warnings    = eval(this.call$show.warnings),
        update           = eval(this.call$update),
        format           = eval(this.call$format),
        data.file        = eval(this.call$data.file),
        save             = eval(this.call$save),
        shiny.input      = eval(this.call$shiny.input)))
  }
  
  
  ##############################################
  ##
  ##   Read from excel file
  ##
  ##############################################
  
  if(class(data) == "character"){
    A.data <- do.call(createMetaDataXls, list(
      excel.file       = eval(this.call$data),
      dbf.path         = eval(this.call$dbf.path),
      protocol.path    = eval(this.call$ protocol.path),
      namevar          = eval(this.call$namevar),
      drugvar          = eval(this.call$drugvar),
      protocolvar      = eval(this.call$protocolvar),
      dbf.file.var     = eval(this.call$dbf.file.var),
      identifier       = eval(this.call$identifier),
      timevar          = eval(this.call$timevar),
      additional.metadata = eval(this.call$additional.metadata),
      incubationvar    = eval(this.call$incubationvar),
      doublingvar       = eval(this.call$doublingvar),
      correctionname   = eval(this.call$correctionname),
      unit             = eval(this.call$unit),
      idvar            = eval(this.call$idvar),
      show.warnings    = eval(this.call$show.warnings),
      update           = eval(this.call$update),
      format           = eval(this.call$format),
      save             = eval(this.call$save),
      data.file        = eval(this.call$data.file),
      shiny.input      = eval(this.call$shiny.input)))
    
  }
  
  
  
  
  
  
  ##############################################
  ##
  ##   Read from data.fram data.frame
  ##
  ##############################################
  
  if(class(data) == "data.frame"){
    A.data <- do.call(createMetaDataDataFrame, list(
      metadata         = eval(this.call$data),
      dbf.path         = eval(this.call$dbf.path),
      protocol.path    = eval(this.call$protocol.path),
      namevar          = eval(this.call$namevar),
      drugvar          = eval(this.call$drugvar),
      protocolvar      = eval(this.call$protocolvar),
      dbf.file.var     = eval(this.call$dbf.file.var),
      identifier       = eval(this.call$identifier),
      timevar          = eval(this.call$timevar),
      additional.metadata = eval(this.call$additional.metadata),
      incubationvar    = eval(this.call$incubationvar),
      doublingvar       = eval(this.call$doublingvar),
      correctionname   = eval(this.call$correctionname),
      unit             = eval(this.call$unit),
      idvar            = eval(this.call$idvar),
      show.warnings    = eval(this.call$show.warnings),
      update           = eval(this.call$update),
      format           = eval(this.call$format),
      save             = eval(this.call$save),
      data.file      = eval(this.call$data.file ),
      shiny.input      = eval(this.call$shiny.input)))
  }
  
  ##############################################
  ##
  ##   Read from data.fram list
  ##
  ##############################################
  
  if(class(data) == "list"){
    this.call$format <- "long" 
    A.data <- do.call(createMetaDataList, list(
      list             = eval(this.call$data),
      dbf.files        = eval(this.call$dbf.files),
      protocol.files   = eval(this.call$protocol.files),
      are.paths.full   = eval(this.call$are.paths.full), 
      dbf.path         = eval(this.call$dbf.path),
      protocol.path    = eval(this.call$protocol.path),
      namevar          = eval(this.call$namevar),
      drugvar          = eval(this.call$drugvar),
      protocolvar      = eval(this.call$protocolvar),
      dbf.file.var     = eval(this.call$dbf.file.var),
      identifier       = eval(this.call$identifier),
      timevar          = eval(this.call$timevar),
      additional.metadata = eval(this.call$additional.metadata), 
      incubationvar    = eval(this.call$incubationvar),
      doublingvar      = eval(this.call$doublingvar),
      correctionname   = eval(this.call$correctionname),
      unit             = eval(this.call$unit),
      idvar            = eval(this.call$idvar),
      format           = eval(this.call$format),
      show.warnings    = eval(this.call$show.warnings),
      update           = eval(this.call$update),
      save             = eval(this.call$save),
      data.file        = eval(this.call$data.file),
      shiny.input      = eval(this.call$shiny.input)))
  }
  
  ##############################################
  ##
  ##   End of function
  ##
  ##############################################
  
  
  call2 <- match.call()
  
  this.call <- list()
  
  myfor <- formals(createMetaData)               ## formals with default arguments
  for (v in names(myfor)){
    if (!(v %in% names(call2))){
      this.call[[v]] <- myfor[[v]]
    }else{
      this.call[[v]] <- call2[[v]]
    }  ## if arg is missing I add it
  }
  this.call$shiny.input <- NULL
  
  
  vec <- vector()
  for(i in 1:length(this.call))
    vec[i] <- ifelse(is.character(this.call[[i]]) & length(this.call[[i]]) == 1, '"', "")
  
  this.shiny.call <- 
    paste('createMetaData', '(', paste(names(this.call), ' = ', vec, 
                                       this.call, vec, collapse = ', ', sep =""), ')', 
          collapse = '', sep = '')
  
  
  
  this.call$call <- call2
  #############
  ## add extension
  w.ext <- substr(basename(data.file), 
                  start=nchar(basename(data.file))- 5, 
                  nchar(basename(data.file))) == ".RData"
  if(!w.ext)
    data.file <- paste(data.file, ".RData", sep = "")
  this.call$data.file <- data.file
  
  
  A.data$call[["createMetaData"]] <- this.call
  
  A.data$call$record <- callNumbering("createMetaData", A.data$call$record)
  class(A.data)  <- c("createMetaData", class(A.data) %w/o% "createMetaData")
  
  if(is.null(A.data$auxiliary$passed.var)){
    A.data$auxiliary$passed.var <- this.call
  }else{    
    for(input.iter in (names(this.call) %w/o% "call"))
      A.data$auxiliary$passed.var[[input.iter]] <-
      this.call[[input.iter]] 
  }
  
  if(is.null(A.data$auxiliary$shiny.calls)){
    shiny.calls <- list(createMetaData=this.shiny.call)
    A.data$auxiliary$shiny.calls <- shiny.calls
  }else{
    A.data$auxiliary$shiny.calls$createMetaData <- this.shiny.call
  }
  
  
  
  old <- A.data
  if(save)
    save(old, file = data.file)
  return(A.data)
}



##
##
## Function for reading dbf data into R
##
##

readDBFData <- function(
  A.data,  # Absorbance, data.frame, list  
  update         = TRUE, 
  perl           = "perl",
  well.id        = "WELLNUM", # well identifier in the dbf file
  absorbance.id  = "M1",      # absorbance identifier in the dbf file 
  discard.lines  = c(1, 2),   # Discarded lines in the dbf file
  dosevar        = "Concentration",
  additivevar    = "Additive",
  controlval     = "Control",
  backgroundval  = "Background", 
  mistakeval     = "X",
  remove.rows    = c("A","H"),# indicate rows that should not be used
  remove.cols    = c(1, 12),  # indicate cols that should not be used
  progressbar    = "text",
  verbose        = FALSE,      # print output
  save           = TRUE,
  shiny.input    = NULL,
  session        = NULL
) {  
  slow <- FALSE
  data <- A.data
  ###############################
  ## create call with formals included
  call2 <- match.call()
  
  this.call <- createCall(call = call2, fun = "readDBFData")#match.call())
  
  
  #############################################
  ##
  ##     Create the call based on the A.data
  ##
  #############################################
  
  A.data.metadata <- NULL
  
  if("createMetaData" %in% class(data)){
    file.info <- data$auxiliary$file.info
    A.data.metadata <- data
    
    CMD.call <- data$call$createMetaData
    
    dbf.path       <- this.call$dbf.path <- CMD.call$dbf.path
    idvar          <- this.call$idvar <-    CMD.call$idvar
    namevar        <- this.call$namevar <-  CMD.call$namevar
    drugvar        <- this.call$drugvar <-  CMD.call$drugvar
    timevar        <- this.call$timevar <-  CMD.call$timevar
    correctionname <- this.call$correctionname <- CMD.call$correctionname
    incubationvar  <- this.call$incubationvar<-  CMD.call$incubationvar
    data           <- A.data.metadata$meta.list$metadata.long
    identifier     <- this.call$identifier <-  CMD.call$identifier 
    
    format <- "long"#this.call$format <- CMD.call$format
    dbf.file.var <- this.call$dbf.file.var <-  eval(CMD.call$dbf.file.var)
    
    protocol.path <-  this.call$protocol.path <-eval(CMD.call$protocol.path)
    
    protocolvar <- this.call$protocolvar <- eval(CMD.call$protocolvar)
    
    data.file <- this.call$data.file <- CMD.call$data.file
    
  } 
  
  
  w.ext <- substr(basename(data.file), 
                  start=nchar(basename(data.file))- 5, 
                  nchar(basename(data.file))) == ".RData"
  if(!w.ext)
    data.file <- paste(data.file, ".RData", sep = "")
  this.call$data.file <- data.file
  
  
  #######################################
  #incubation is added to the established dataset
  
  if(!incubationvar %in% data)
    data[,incubationvar] <- 2
#   if(!is.null(shiny.input)){
#     progressbar = "none"
#   }
  
  if(!progressbar == "none"){
    on.exit(close(pb))
    pb <- progressBar(title = "Reading dbf files", min = 0,
                      max = nrow(data), width = 300, window = progressbar == "window" )
  }
  
#  pb <- DoseR:::progressBar(title = "Reading dbf files", min = 0,
#                    max = 200, width = 300, window = TRUE)
  
  #############################################
  ##
  ##     Create a record that keeps track of changes
  ##     in the various functions using the A.data
  ##
  #############################################
  prior <- FALSE
  priortest <- all(file.exists(data.file), update)
  
  if(file.exists(data.file)){
    load(data.file)
  }
  if(priortest){
    if("raw.data" %in% names(old$data))
      prior <- TRUE
  }
  
  if(prior)
    if("raw.data.all" %in% names(old$data))
      if(old$call$record["readDBFData", "last.visited"] < 
           old$call$record["drugColorCorrection", "last.visited"])
        old$data$raw.data <- old$data$raw.data.all 
  
  ##############################################
  ##
  ## Check that the old call resembles the new call
  ##
  ######################################
  
  if(prior){
    this.call <- this.call
    old.call  <- old$call$readDBFData
    if(length(intersect(names(old.call) , names(this.call))) != 
         length(names(old.call)) |
         length(intersect(names(old.call) , 
                          names(this.call))) != length(names(this.call) ))
      prior <- FALSE
    vec <- ""
    if(prior){
      call.names <- names(old.call) %w/o% 
        c("call", "format", "update", "data", "verbose", 
          "progressbar", "shiny.input", "shiny.call", "session")
      
      vec <- list()
      for(call.iter in call.names){
        vec[[call.iter]] <- this.call[[call.iter]] == old.call[[call.iter]]
      }
      if(!any(unlist(vec)))
        prior <- FALSE
    }
    if(!prior){
      
      if(createMetaDataDir %in% names(old$call)){
        xls.file <- A.data$auxiliary$passed.var$data.file
        excel.2 <- paste( no.extension(xls.file), " last changed ",
                          file.info(xls.file)$mtime, 
                          ".xls", sep = "")
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
      }else{
        data.file2 <- paste(no.extension(data.file), " last changed ",
                            file.info(xls.file)$mtime, 
                            ".RData", sep = "")
        file.copy(from=data.file, to=data.file2)
        
        warning("The input variables have changed but the data is saved to an allready exisiting project without setting update to FALSE.\n\n", 
                "The input update is changed to FALSE.\n\n",
                "The old R.data file is backed up to:\n",
                data.file2)
      }
    }
  }
  
  #########################################
  ##
  ##  Create the record 
  ##
  ########################################
  
  early <- Sys.time()
  ## Read the data.file into R
  if(prior){
    record <- old$auxiliary$record
    # if(!is.null(correctionname))
    #record <- record[ record[, namevar] != correctionname, ]
    info   <- old$auxiliary$info
  }else{
    record <- data[, c(idvar, namevar, drugvar)]
    record$readDBFData <- Sys.time()
  }
  
  ## change potential factors to characters
  for(cols in c(protocolvar, dbf.file.var)){
    data[, cols] <- as.character(data[, cols])
  }
  
  #############################################
  ##
  ##     Check that all dbf and protocol files exists
  ##
  #############################################
  
  file.info$protocol$exists <- 
    file.exists(file.info$protocol$file.full)
  
  if(!all(file.info$protocol$exists)){
    print(file.info$protocol[!file.info$protocol$exists, 
                             "file.full",  drop = FALSE])
    if(progressbar != "none" & max > 0)
    close(pb)
    stop("protocols listed above not found")
  }
  
  protocol.error <- !data[,protocolvar] %in% row.names(file.info$protocol)
  
  if(any(protocol.error)){
    cat("\n")
    print(unique(data[protocol.error, protocolvar] ))
    if(progressbar != "none")
    close(pb)
    stop("protocols listed above not found")
  }
  
  
  file.info$dbf$exists <- 
    file.exists(file.info$dbf$file.full)
  
  if(!all(file.info$dbf$exists)){
    print(file.info$dbf[!file.info$dbf$exists, 
                        "file.full",  drop = FALSE])
    if(progressbar != "none")
    close(pb)
    stop("dbf files listed above not found")
  }
  
  
  #############################################
  ##
  ##     Remove protocols not used by the 
  ##     current metasheet
  ##
  #############################################
  
  protocol.used <- unique(data[, protocolvar])
  
  file.info$protocol <-
    file.info$protocol[intersect(protocol.used,row.names(file.info$protocol)), ]
  
  #############################################
  ##
  ##     Read the protocols
  ##
  #############################################
  
  ## Read the protocols into the list protocols
  
  protocol.vec <- rownames(file.info$protocol)
  if(prior){
    protocols <- old$data$protocols
    extra <- unique(protocol.vec) %w/o% protocols$record$protocol
    missing <- protocols$record$protocol %w/o% unique(protocol.vec)
    protocols$record <- protocols$record[!(protocols$record$protocol %in% missing),]
    if(length(extra) > 0){
      change2 <-  file.info$protocol[extra, "last.change"]
      
      xx <- data.frame(protocol = as.character(extra), last = NA, 
                       now = change2, change = TRUE)
      # rownames(xx) <- xx$protocol
      protocols$record <- rbind(protocols$record, xx)   
      #protocols$record[is.na(protocols$record$last),"last"] <- 0
    }
  }else{
    protocols <- list()
    protocols$record <- data.frame(protocol = unique(protocol.vec), last = 0)
  }
  
  protocols$record$protocol <- as.character(protocols$record$protocol)
  
  protocols$record$now <- file.info$protocol[protocols$record$protocol, "last.change"]
  if(prior){
    protocols$record$change <- protocols$record$now != protocols$record$last
  }else{
    protocols$record$change <- TRUE
  }
  protocols$record$change[is.na(protocols$record$change)] <- TRUE
  prots <- protocols$record$protocol[protocols$record$change]
  if(length(prots) > 0 ){
    for(prot in prots) {
      filesetup <-  file.info$protocol[prot, "file.full"]
      
      is64 <- grepl("64", sessionInfo()[1]$R.version$arch)
      sos  <- grepl("darwin", sessionInfo()[1]$R.version$os)
      
      if(sos){
        protocols[[prot]][["setup"]] <- data.frame(read_excel(filesetup, sheet = 1))
        protocols[[prot]][["conc"]]  <- data.frame(read_excel(filesetup, sheet = 2))
      }else{
        if(is64){
          protocols[[prot]][["setup"]] <- 
            data.frame(read_excel(filesetup, sheet = 1))
          protocols[[prot]][["conc"]]  <- 
            data.frame(read_excel(filesetup, sheet = 2))
        }else{
          protocols[[prot]][["setup"]] <- 
            data.frame(read_excel(filesetup, col_names = TRUE, sheet = 1))
          protocols[[prot]][["conc"]]  <- 
            data.frame(read_excel(filesetup, col_names = TRUE, sheet = 2))
        }
      }
      
      protocols[[prot]][["conc"]] <- 
        protocols[[prot]][["conc"]][
          ,!(colnames(protocols[[prot]][["conc"]]) %in%
               c("X", paste("X", 1: ncol(protocols[[prot]][["conc"]]), sep = ".")))]
      
      row.names(protocols[[prot]][["conc"]]) <- 
        as.character(protocols[[prot]][["conc"]][, 1])
      protocols[[prot]][["conc"]] <- protocols[[prot]][["conc"]][, -1]
      
      protocols[[prot]][["setup"]] <- 
        protocols[[prot]][["setup"]][
          ,!(colnames(protocols[[prot]][["setup"]])[-1] %in%
               c("X", paste("X", 1: ncol(protocols[[prot]][["setup"]]), sep = ".")))]
      
      row.names(protocols[[prot]][["setup"]]) <- 
        as.character(protocols[[prot]][["setup"]][, 1])
      
      protocols[[prot]][["setup"]] <- protocols[[prot]][["setup"]][, -1]
    }  
  }
  
  protocols$record$last <- protocols$record$now
  
  #############################################
  ##
  ##     Create the raw.data object,
  ##     either from new or from old$data$raw.data
  ##
  #############################################
  
  if(prior){
    raw.data <- old$data$raw.data
    names <- colnames(data)  %w/o% c(timevar, protocolvar, dbf.file.var) 
  }else{
    names <- colnames(data)  %w/o% c(timevar, protocolvar, dbf.file.var) 
    #names <- c(names)
    raw.data <- (matrix(NaN, ncol = length(names) ))
    colnames(raw.data) <- names
    raw.data <- as.data.frame(raw.data)  
  }
  
  ## remove the platesets which, has been removed from metadata
  if(prior){
    rem <- old$meta.list$readDBFData[, idvar] %w/o% data[, idvar]
    if(length(rem) > 0){
      raw.data <- raw.data[!(raw.data[,idvar] %in% rem) ,]
      record <- record[!(record[, idvar] %in%  rem),]
      
      
      rem.level <- unique(old$data$raw.data[old$data$raw.data[, idvar]  %in% rem, "level"])
      old$auxiliary$bootstrap.levels <- 
        old$auxiliary$bootstrap.levels[!(old$auxiliary$bootstrap.levels$level %in% rem.level),]
    }
  }
  
  
  if(prior){
    old.meta <- old$meta.list$readDBFData    
    row.names(old.meta) <- old.meta[,idvar]
    
    new.meta <- data
    row.names(new.meta) <- new.meta[,idvar]
    int.id <- intersect(row.names(new.meta), row.names(old.meta))
    wh.drug <- old.meta[int.id, drugvar] != new.meta[int.id, drugvar]
    wh.name <- apply(old.meta[int.id, namevar, drop = FALSE] != 
                       new.meta[int.id, namevar, drop = FALSE], 1, any)
    
    
    
    h1 <- old.meta[int.id, timevar, drop = FALSE]
    h1[is.na(h1)] <- ""
    
    h2 <- new.meta[int.id, timevar, drop = FALSE]
    h2[is.na(h2)] <- ""
    
    wh.time <- apply(h1 != h2, 1, any)
    
    
    h1 <- old.meta[int.id, dbf.file.var, drop = FALSE] 
    h1[is.na(h1)] <- ""
    
    h2 <- new.meta[int.id, dbf.file.var, drop = FALSE]
    h2[is.na(h2)] <- ""
    
    wh.dbf <- apply(h1 != h2, 1, any)
    
    changed.data <- apply(cbind(wh.time, wh.dbf), 1, any)
    
    xx <- rbind(new.meta[wh.drug, c(namevar, drugvar)],
                old.meta[wh.drug, c(namevar, drugvar)])
    
    changed.data[new.meta[new.meta[,drugvar] %in% xx[,drugvar] &
                            new.meta[,namevar] %in% xx[,namevar]  , idvar]] <- TRUE       
    
    if("T0.list" %in% names(old)){    
      rem.drug <- names(old$data$T0.list) %w/o% unique(new.meta[, drugvar])
      for(i in rem.drug)  
        T0.list[[i]] <- NULL
      for(i in names(old$data$T0.list)){
        rem.times <- names(old$data$T0.list[[i]]) %w/o% 
          unique(new.meta[new.meta[, drugvar] == i, timevar])
        rem.times %w/o% "all"
        for(j in rem.times)
          T0.list[[i]][[j]] <- NULL
      }
    }
    xx <- data[changed.data, c(namevar, drugvar)]
    ids <- data[data[, namevar] %in% xx[, namevar] &  
                  data[, drugvar] %in% xx[, drugvar], idvar]
    changed.data[ids] <- TRUE
    colnames(data)
    raw.data <- raw.data[!(raw.data[, idvar] %in% 
                             names(changed.data)[changed.data]),]  
  }
  
  
  raw.data <- raw.data[, colnames(raw.data) %w/o% 
                         c(backgroundval, controlval, "plateset", "plate", "level", "new.type" )] 
  
  names <- colnames(data) # %w/o% c(timevar, protocolvar, dbf.file.var) 
  data.list <- list()
  #n.plates <- 0
  record2 <- vector()
  for(i in 1:nrow(data)) { 
    #dbf.file.var2 <- (1:length(dbf.file.var))[!(data[i, dbf.file.var] %in% c("", NA, NaN))]
    change <- TRUE# rep(TRUE, length(dbf.file.var2))
    
    prot <- data[i, protocolvar]
    setup <- protocols[[prot]][["setup"]]
    conc  <- protocols[[prot]][["conc"]]
    
    ## Have the data allready been treated ?
    
    if(prior) {
      change <- FALSE    
      if(!(data[i, idvar] %in% old$meta.list$readDBFData[,idvar])) ## dette er old$id er korrekt
        change <- TRUE
      
      if(!change){
        if(changed.data[data[i, idvar]])
          change <- TRUE 
        
        if(!change & !(data[i, protocolvar] %in% rownames(old$auxiliary$file.info$backup$protocol)))
          change<- TRUE
        if(!change & data[i, protocolvar] !=
             old$meta.list$readDBFData[old$meta.list$readDBFData[, idvar] == data[i, idvar], protocolvar])
          change <- TRUE
        if(data[i, dbf.file.var] %in% rownames(old$auxiliary$file.info$backup$dbf))
          if(!change & file.info$dbf[data[i, dbf.file.var], "last.change"] != 
               old$auxiliary$file.info$backup$dbf[data[i, dbf.file.var], "last.change"])
            change <- TRUE
        if(!change)
          change <-
          any(!identical(protocols[[prot]][["conc"]] ,
                         old$data$protocols[[prot]][["conc"]]),
              !identical(protocols[[prot]][["setup"]],
                         old$data$protocols[[prot]][["setup"]]))
      }
    }
    
    
    
    ## when new data is added or there is changes in the protocols or dbf files
    ## the data is read in again
    if(change) {
      record2 <- c(record2, data[i, idvar])
      ## old data with same ID is removed and replaced with the new
      raw.data <- raw.data[raw.data[, idvar] != data[i, idvar],  ]
      
      #n.plates <- n.plates + 1
      setup <- protocols[[data[i, protocolvar]]][["setup"]]
      conc  <- protocols[[data[i, protocolvar]]][["conc"]]
      
      a <- paste(rownames(setup), rep( as.numeric(as.numeric(gsub("[^[:digit:]]","",colnames(setup)))),
                                       each = nrow(setup)), sep = "")
      
      ## each column is converted into a character
      for(k in 1:ncol(setup))
        setup[, k] <- as.character(setup[, k])
      
      ## data frame with the rownames and content is created
      b <- as.character(unlist(setup))
      concsetup <- data.frame(a, b)
      rownames(concsetup) <- concsetup$a
      
      ## The plate is read into R
      if(verbose)
        cat(data[i, dbf.file.var], "\n")
      
      cell0        <- read.dbf(file.info$dbf[data[i, dbf.file.var], "file.full"],
                               as.is=TRUE)
      missingcolumns <- c(well.id, absorbance.id) %w/o% colnames(cell0)
      
      if(length(missingcolumns )){
        well.id.mis <- !well.id %in% colnames(cell0) 
        head(cell0)
        if(well.id.mis)
          cat("The given well.id (", well.id  ,") is not found in colnames of the dbf file.\n", sep = "")
        absorbance.id.mis <- !absorbance.id %in% colnames(cell0) 
        if(absorbance.id.mis)
          cat("The given absorbance.id (", absorbance.id  ,") is not found in colnames of the dbf file.\n", sep = "")      
        cat("The header of the dbf file is written above. Check which colnames are appropriate for well.id and absorbance.id\n")
        cat("current dbf file:" , data[i,dbf.file.var] ,"\n")
        if(progressbar != "none")
        close(pb)
        stop()
      }
    
      if(length(discard.lines) > 0)
        cell0        <- cell0[ - discard.lines, ]
      
      cell0$colnum <- as.numeric(as.numeric(gsub("[^[:digit:]]","",cell0[, well.id])))
      cell0$rownum <-  gsub("[[:digit:]]","",cell0[, well.id])
      
      cell0$wellnum <- paste(cell0$rownum, cell0$colnum, sep = "")
      rownames(cell0) <- cell0$wellnum
      
      cell0$type <- as.character(concsetup[rownames(cell0),]$b)
      
      cell0 <- cbind(cell0, (conc[as.character(cell0$type), ]))
      
      
      cell0 <- cell0[!(cell0$colnum %in% remove.cols), ]
      cell0 <- cell0[!(cell0$rownum %in% remove.rows), ]
      
      
      orders <- c(colnames(conc), "type")
      for(order in rev(orders))
        cell0 <- cell0[order(cell0[, order]), ]
      
      cell0$doserep <- 1
      iter <- 1
      dups <- duplicated(paste(cell0$type, cell0$doserep))
      
      while(any(dups)){
        iter <- iter +1
        cell0$doserep <- ifelse(!dups, cell0$doserep, iter)
        dups <- duplicated(paste(cell0$type, cell0$doserep))
      }
      
      cell0$absorbance <- cell0[, absorbance.id]
      
      ##"plateset", "plate"
      for(k in c(timevar, "change", "type",
                 colnames(conc), "absorbance", "wellnum", "rownum", "colnum",
                 dbf.file.var, protocolvar)){
        if(!(k %in% colnames(raw.data)) )
          raw.data[, k] <- NaN
      }
      suppressWarnings(
        xx <- cbind(data[i, names],
                    timevar = data[i, timevar],
                    change = any(change),
                    cell0[, c("type", colnames(conc), "absorbance", "wellnum",
                              "rownum", "colnum")],
                    dbf.file.var = data[i, dbf.file.var],
                    protocolvar = data[i, protocolvar],
                    absorbance.nc = NaN)
      )
      if(slow){
        raw.data <- rbind(raw.data, xx[, colnames(raw.data)])
      }else{
        ek <- colnames(raw.data) %w/o% colnames(xx)
        xx[, ek] <- NA
        data.list[[length(data.list) + 1]] <- xx[, colnames(raw.data)]
      }
      
      #options(warn=-1)
    }
    if(!progressbar == "none")
      setProgressBar(pb, i, label=paste( round(i/nrow(data)*100, 0),
                                         "% done"),
                     window = progressbar == "window")
  #  if(!is.null(shiny.input))
  #  setProgress(value = i,
  #              detail=paste( round(i/nrow(data)*100, 0),
  #                            "% done"))
    
    
  }
  
  
  new <- record2 %w/o% record[,idvar]
  if(length(new) > 0){
    for(news in new){
      for(i in 4:length(colnames(record)))
        record[news, i] <- as.POSIXct("2011-04-04 14:18:58", tz="GB")
      record[news,  c(idvar, namevar, drugvar)] <- data[data[,idvar] == news, c(idvar, namevar, drugvar)] 
    }   
  }
  
  
  record$readDBFData[record[,idvar] %in% record2] <- Sys.time()
  row.names(record) <- record[, idvar]
  
  
  raw.data <- rbind(raw.data, rbind.fill(data.list))
  if(raw.data[1, idvar ] == NaN) raw.data <- raw.data[-1, ]
  raw.data[,backgroundval] <- ifelse(raw.data[, additivevar] == backgroundval, 1, 0)
  raw.data$control    <- ifelse(raw.data[, additivevar] == controlval, 1, 0)
  
  # remove those marked by a capital X as errors
  raw.data <- raw.data[raw.data$type != mistakeval,]
  
  if(is.null(A.data.metadata)){
    raw.data$plateset <- raw.data[, idvar]
    raw.data$plate <- paste(raw.data[, idvar],  raw.data[, timevar])
    raw.data$level <- paste(raw.data[, drugvar], raw.data[, namevar], raw.data[, timevar])
  }else{
    raw.data$plateset <- paste(raw.data[, drugvar], raw.data[, namevar], raw.data[, identifier])
    raw.data$plate <- paste(raw.data[, drugvar], raw.data[, namevar], raw.data[, identifier],
                            raw.data[, timevar])
    raw.data$level <- paste(raw.data[, drugvar], raw.data[, namevar], raw.data[, timevar])
  }
  ## problem that the colnames need to state "Concentration" and "Additive"
  raw.data$new.type <- paste(raw.data[, dosevar], raw.data[, additivevar])
  
  
  if(!prior){
    info <- list()
    for(drug.iter in unique(data[, drugvar]))
      info[[drug.iter]] <- drugInfo(drug.iter)
    
    mol.data <- data.frame(drug = "drug", mol.mass = "mol.mass", stringsAsFactors=FALSE)
    for(drug.iter in names(info))
      if(any(info[[drug.iter]] != "Drug not found"))
        mol.data[drug.iter, ] <- 
      c(drug.iter, info[[drug.iter]]$"Chemical data"["Molar mass", ])
    
    mol.data <- mol.data[-1,]
    
  }else{
    mol.data <- old$auxiliary$mol.data
    remove <- names(info) %w/o% unique(data[, drugvar])
    
    for(i in remove){
      info[[i]] <- NULL
    }
    
    for(drug.iter in (unique(data[, drugvar]) %w/o% names(info))){
      info[[drug.iter]] <- drugInfo(drug.iter)
      if(any(info[[drug.iter]] != "Drug not found"))
        mol.data[drug.iter, ] <- 
        c(drug.iter, info[[drug.iter]]$"Chemical data"["Molar mass
                                                       
                                                       ", ])
    }
    
    mol.data <- mol.data[unique(data[, drugvar]), ]   
  }
  mol.data[,2] <- gsub("\\s[g]", "", mol.data[,2])    
  mol.data[,2] <- as.numeric(gsub("/mol",   "", mol.data[,2]))
    
  if(!file.exists(data.file)){
    data.b <- data
    call <- list()
    call[["readDBFData"]] <- this.call    
    call$record <- callNumbering("readDBFData")
    meta.list <- list()
    
    file.info$backup <- file.info
    meta.list[["readDBFData"]] <- data
    auxiliary <- list(file.info = file.info, record = record, info = info,  
                      mol.data =  mol.data, shiny.input = shiny.input)
    
    data <- list(raw.data = raw.data, protocols = protocols)
    old <- list(data = data, meta.list = meta.list, call = call)
    data <- data.b
  }
  
  if(file.exists(data.file)){
    file.info$backup <- file.info
    old$call$record <- callNumbering("readDBFData", record = old$call$record)
    old$call[["readDBFData"]] <- this.call  
    
    old.meta <- old$meta.list$readDBFData
    row.names(old.meta) <- old.meta[, idvar]
    
    old$data$raw.data     <- raw.data
    old$meta.list$readDBFData     <- data
    old$data$protocols    <- protocols
    # old$protocol.vec <- protocol.vec
    old$auxiliary$record       <- record
    old$auxiliary$info         <- info
    old$auxiliary$file.info    <- file.info
    old$auxiliary$mol.data     <- mol.data
    old$auxiliary$shiny.input <- shiny.input
  }
  
  
  if(prior) {
    
    new.meta <- old$meta.list$readDBFData
    row.names(new.meta) <- new.meta[,idvar, drop = T]
    
    com.rows <- intersect(old.meta[, idvar], new.meta[, idvar])
    
    old.meta <- old.meta[com.rows, ]
    new.meta <- new.meta[com.rows, ]
    
    cols <- intersect(colnames(old.meta), colnames(new.meta))
    cols <- unique(c(cols, colnames(new.meta))) %w/o% c(idvar, "replicate", "reason", "omit")
    for(col in cols){
      
      if(col %in% (colnames(old.meta) %w/o% "replicate")){
        wh.row <- old.meta[, col] != new.meta[, col]
        wh.row[is.na(old.meta[, col]) | is.na(new.meta[, col]) ] <- TRUE
        #wh.row[is.na(old.meta[, col]) & is.na(new.meta[, col]) ] <- FALSE
      }else{
        wh.row <- rep(TRUE, nrow(new.meta))
      }
      if(any(wh.row)){
        
        for(old.iter in (names(old) %w/o% c("meta.list", "call", "drug.color.correct"))){
          data.frames <- names(old[[old.iter]]) %w/o% 
            c("metadata", "e.data", "protocols", 
              "protocol.vec", "id",  "fits", "T0.mat",
              "bootstrap.levels", "T0.list", "info", "mol.data",
              "DR.model", "shiny.input")
          for(frame in data.frames){
            
            if(frame %in% "record" & old.iter == "auxiliary"){
              if(col %in% colnames(old[[old.iter]][[frame]])){
                wh.row.new <- old[[old.iter]][[frame]][, idvar ] %in% new.meta[wh.row , idvar]
                id.v <- old[[old.iter]][[frame]][, idvar ][wh.row.new]       
                old[[old.iter]][[frame]][wh.row.new, col] <- new.meta[id.v, col]
              }
            }else{         
              if(class(old[[old.iter]][[frame]])[1] == "data.frame"){
                wh.row.new <- old[[old.iter]][[frame]][, idvar ] %in% new.meta[wh.row , idvar]
                id.v <- old[[old.iter]][[frame]][, idvar ][wh.row.new]       
                old[[old.iter]][[frame]][wh.row.new, col] <- new.meta[id.v, col]
              }
            }
          }
        }      
      }    
    }
    
    rem.cols <- colnames(old.meta) %w/o% colnames(new.meta)
    if(!is.null(rem.cols)){
      for(old.iter in c("data") ){
        data.frames <- names(old[[old.iter]]) %w/o% 
          c("metadata", "e.data", "protocols", 
            "protocol.vec", "id", "record", "fits",
            "bootstrap.levels", "T0.list", "info", "mol.data",
            "DR.model", "shiny.input", "iso.fits")
        for(frame in data.frames){
          for(rem.col in rem.cols) 
            old[[old.iter]][[frame]] <- 
            old[[old.iter]][[frame]][, (colnames(old[[old.iter]][[frame]]) %w/o% rem.col)]
        }    
      }
    }
  }
  
  
  if(is.null(old$auxiliary$passed.var)){
    old$auxiliary$passed.var <- this.call
  }else{    
    for(input.iter in names(this.call)[-1])
      old$auxiliary$passed.var[[input.iter]] <-
      this.call[[input.iter]] 
  }
  
  class(old) <- c("Absorbance", union(class(old),
                                      class(A.data.metadata)) %w/o% "Absorbance")
  
  if(save)
    save(old, file = data.file)
  
  #close(pb)
  

  return(old)
  if(progressbar != "none")
  invisible()               
}


##
##
##     Function for fitting the nonlinear function dA + b
##
##


bgModel <- function(A.data, #..., 
                    update          = TRUE,
                    parametrisation = "unrestricted", 
                    outlier.test    = 3,
                    outlier.iter    = 2,
                    weights         = "fitted", 
                    fitted.a        = FALSE,
                    varpower.min    = 10e-4,
                    varpower.iter   = 50,
                    contr           = c("sum", "helmert", "treatment"),
                    progressbar     = "text",
                    verbose         = FALSE,
                    save            = TRUE,
                    shiny.input     = NULL,
                    session         = NULL){
  

  
  ###########################################
  ##
  ## The Call is saved
  ##
  ###########################################
  
  if(!parametrisation %in% c("unrestricted", "restricted"))
  stop("Parametrisation need to be either restricted and unrestricted")
  if(parametrisation == "unrestricted")
    parametrisation <- "none"
  if(parametrisation == "restricted")
    parametrisation <- c("square", "abs", "exp", "expinv")
  
  
  ###############################
  ## create call with formals included
  this.call <- createCall(call = match.call(), "bgModel")
  
  
  
  ###########################################
  ##
  ## update the call according to read data
  ##
  ###########################################
  
  
  prior <- FALSE
  
  drug.color.correction <- "correction.data" %in% class(data)
  
  if("Absorbance" %in% class(A.data)){
    call         <- A.data$auxiliary$passed.var
    A            <- "absorbance"
    B            <- eval(call$backgroundval)
    backgroundval<- eval(call$backgroundval)
    idvar        <- eval(call$idvar)
    drugvar      <- eval(call$drugvar)
    namevar      <- eval(call$namevar)
    timevar      <- eval(call$timevar)
    dbf.file.var <- eval(call$dbf.file.var)
    additivevar  <- eval(call$additivevar)
    dosevar      <- eval(call$dosevar)
    plateset     <- "plateset"
    plate        <- "plate"
    type         <- "new.type"
    record       <- A.data$auxiliary$record
    data.file    <- eval(call$data.file) 
    
    prior <- all("bc.data" %in% names(A.data$data), update)
    if(prior){
      data   <- A.data$data$bc.data
    }else{
      data   <- A.data$data$raw.data
    }
  }
  
  data <- data[!is.na(data[namevar]),]
  ## Changes compared to the previous version
  if(prior) { 
    ## Changes compared to old
    if("drugColorCorrection" %in% colnames(record) ){
      new.platesets <-  record[record$bgModel < record$readDBFData |
                                 record$bgModel < record$drugColorCorrection, idvar]
    }else{
      new.platesets <-  record[record$bgModel < record$readDBFData, idvar]
    }
    new.levels    <- unique(A.data$data$raw.data[
      A.data$data$raw.data[, idvar] %in% new.platesets, "level"])
    
    ## removed compared to old
    rem.platesets <- unique(A.data$data$bc.data[, idvar]) %w/o% 
      unique(A.data$data$raw.data[, idvar])
    rem.levels <- unique(A.data$data$bc.data[
      A.data$data$bc.data[, idvar] %in% rem.platesets, "level"])
    
    level.relevant <- unique(A.data$data$bc.data[A.data$data$bc.data$level %in% rem.levels, 
                                                 "plateset"]) %w/o% rem.platesets
    
    
    rem.dbf <- unique(A.data$data$bc.data[, dbf.file.var]) %w/o% 
      unique(A.data$data$raw.data[, dbf.file.var])
    
    ## if an replicate is removed a new background correction takes place.
    # if(length(rem.platesets) > 0){
    
    #data <- data[!(data$level %in% cange.levels), ]
    
    
    ## if an replicate is removed a new background correction takes place.
    data.2      <- A.data$data$bc.data[!(A.data$data$bc.data[, idvar] %in% rem.platesets), ]
    rem.levels2 <- intersect(rem.levels, data.2$level)
    complete.rem.levels <- setdiff(rem.levels, rem.levels2)
    new.levels  <- c(new.levels, rem.levels2)
    
    if(length(new.platesets) > 0){
      data.new <- A.data$data$raw.data[A.data$data$raw.data[,idvar] %in% new.platesets, ]    
      int.col <- intersect(colnames(data.2), colnames(data.new))
      data.2[(nrow(data.2)+1):(nrow(data.2)+nrow(data.new)), int.col] <- data.new[, int.col]
    }
    
    data <- data.2
    
    data <- data[!(data[, dbf.file.var] %in% rem.dbf),]
    
    change.platesets <- c(new.platesets, rem.platesets )    
    change.levels <- c(new.levels, rem.levels)
    
    fits <- A.data$fits$bgModel
    for(remove in change.levels){  
      grepl(remove, names(A.data$fits$bgModel))
      fits <- fits[!grepl(remove, names(fits))]
    }
    
    if( length(complete.rem.levels) > 0){
      e.data <- A.data$auxiliary$e.data[!(A.data$auxiliary$e.data$level %in% complete.rem.levels),]    
    }else{
      e.data <- A.data$auxiliary$e.data
    }
    
    new.levels2 <- new.levels %w/o% e.data$level
    if(length(new.levels2) > 0){
      e.data2 <- data.frame(level = new.levels2, error = NaN)
      row.names(e.data2) <- e.data2$level
      e.data <- rbind(e.data, e.data2)
    }
    
    new.platesets <- new.levels
    
  } else{
    ## lists for errors and fits
    record$bgModel <- Sys.time()
    e.data <- data.frame(level = unique(data$level), error = NaN)
    row.names(e.data) <- e.data$level
    data$para  <- NA
    fits        <- list()
    new.platesets <- unique(data$level)
  }
  
  contr.orig <- contr  
  
  ## background correct according to plate in naive
  if(!is.null(B)){
    mean.b <- aggregate(formula(paste(A, "~" ,plate)), FUN = mean,
                        data = data[data[, backgroundval] == 1,])
    rownames(mean.b) <- mean.b[, plate]
    data[, "NB"] <- data[,A] - mean.b[data[,plate], A]
  }else{
    data$NB <- data[,A]
  }
  
  ## Truncate NB at zero
  data$NB[data$NB < 0] <- 0
  
  lev.table <- table(paste(data$level[!duplicated(paste(data$plateset, data[,timevar]))] ))
  
  prob.list <- list()
  
  ids <-  unique(A.data$data$raw.data[A.data$data$raw.data$level %in% new.platesets, idvar])
  if(length(ids) > 0 )
    record[ids, "bgModel"] <- Sys.time()
  
  if(!progressbar == "none" & length(new.platesets) > 0)
    pb <- progressBar(title = "Normalizing absorbance data ", min = 0,
                      max = length(new.platesets), width = 300, 
                      window = progressbar == "window")
  
  count <- 0
  for(platesets in new.platesets) {
    
    
    temp <- data[data$level == platesets,]
    repli <- lev.table[platesets]
    
    fit.back2 <- "error"
    para.iter <- 0
    while(class(fit.back2)[1] == "character" & para.iter <= length(parametrisation)){
      para.iter <- para.iter + 1
      #B2 <- backgroundval
      #print(varpower.min)
      suppressWarnings(
      fit <- BGfunction(A = A, B = backgroundval, type = type,
                        plate = plate,
                        plateset = plateset,
                        platesets = platesets,
                        data = data, #...,
                        temp = temp,
                        weights = weights, fitted.a = fitted.a,
                        contr = contr.orig[1],
                        outlier.test  = outlier.test,
                        outlier.iter  = outlier.iter,
                        varpower.min  = varpower.min,
                        varpower.iter = varpower.iter,
                        parametrisation = parametrisation[para.iter],
                        verbose = verbose)
      )
      #if(fit$fit2[1] == "error")
      fit.back2 <- fit$fit
    }
    
    unique(temp$plate)
    
    #prob <- data.frame(type =temp[temp$plate ==  unique(temp$plate)[1], "type" ])
    #for(pl in as.character(unique(temp$plate))){
    #  prob[, paste(pl, "A")] <- temp[temp$plate ==  pl, "absorbance"]
    #  prob[, paste(pl, "NB")] <- temp[temp$plate ==  pl, "NB"]
    #}
    #prob.list[[platesets]][["plates"]] <- prob
    #prob.list[[platesets]][["plates.table"]] <- table(temp$plate, temp$id)
    temp <- fit$temp
    
    #aggregate(A~new.type, FUN = mean, data = temp)
    
    if(class(fit.back2)[1] == "character"){
      e.data[platesets, "error"] <-  "error"
      
      rows <- rownames(temp[temp$outlier == 0,])
      data[rownames(temp), "outlier"] <- temp$outlier
      data[rows, "BC"] <- NaN
      data[rows, c("value")] <- NaN
      data[rows, c("Std.Error")] <- NaN
    }else{   
      e.data[platesets, "error"] <-  "succes"
      fits[[platesets]] <- fit.back2
      fits[[paste(platesets, "col", sep = ":")]] <- fit$col
      fits[[paste(platesets, "pch", sep = ":")]] <- fit$pch
      fits[[paste(platesets, "start", sep = ":")]] <- fit$start
      a.coef <- coef(fit.back2)
      a.coef <- a.coef[grepl("a.content2", names(a.coef))]
      names(a.coef) <- gsub("a.content2", "", names(a.coef))
      rows <- rownames(temp[temp$outlier == 0,])
      data[rownames(temp), "outlier"] <- temp$outlier
      
      data[rows, "para"] <- parametrisation[para.iter]
      data[rownames(temp), "BC2"]  <- as.numeric(temp$A)
      
      data[rows, "BC"]   <- switch(parametrisation[para.iter],
                                   none   = a.coef[data[rows, type]],
                                   abs    = abs(a.coef[data[rows, type]]),
                                   square = (a.coef[data[rows, type]])^2,
                                   exp    = exp(a.coef[data[rows, type]]),
                                   expinv = exp(1/a.coef[data[rows, type]]))
      
      
      a.conf <- summary(fit.back2)$tTable
      a.conf <- a.conf[grepl("a.content2", rownames(a.conf)),]
      rownames(a.conf) <- gsub("a.content2", "", rownames(a.conf))
      
      data[rows, c("value")] <-
        a.conf[,1][data[rows, type]]
      
      data[rows, c("Std.Error")] <-
        a.conf[,2][data[rows, type]]
      
      if(parametrisation[para.iter] == "none"){
        lo <- data[rows, c("value")] - 1.96 * data[rows, c("Std.Error")]
        up <- data[rows, c("value")] + 1.96 * data[rows, c("Std.Error")]
      }else{
        up <- abs(data[rows, c("value")]) + 1.96 * data[rows, c("Std.Error")]
        lo <- abs(data[rows, c("value")]) - 1.96 * data[rows, c("Std.Error")]
      }       
      
      lo <- switch(parametrisation[para.iter],
                   none   = lo,
                   abs    = lo,
                   square = lo^2,
                   exp    = exp(lo),
                   expinv = exp(1 / lo))
      
      # if(parametrisation[para.iter] != "none"){
      data[rows, c("lo")] <- ifelse(lo <= 0, 0, lo)
      
      up <- switch(parametrisation[para.iter],
                   none = up,
                   abs = up,
                   square = up^2,
                   exp = exp(up),
                   expinv = exp(1/up))
      data[rows, c("up")] <- up
    }
    if(!progressbar == "none" & length(new.platesets) > 0){
      count <- count + 1
      setProgressBar(pb, count, label=paste(round(count/length(new.platesets)*100, 1),
                                            "% done"),
                     window = progressbar == "window")
    }
  }
  
  
  
  A.data$auxiliary$e.data <- e.data
  if(sum(A.data$auxiliary$e.data$error == "error") > 0 )
    print(A.data$auxiliary$e.data[A.data$auxiliary$e.data$error == "error",])
  
  mean.data <- aggregate(cbind(absorbance, NB, BC2) ~ new.type + level,
                         FUN = mean, data = data[ data$outlier != 1 ,])
  
  mean.data <- (merge(data, as.data.frame(mean.data),
                      by = c("new.type", "level"), all = TRUE))
  
  mean.data <- mean.data[, !(colnames(mean.data) == "NB.x")]
  #mean.data <- mean.data[, !(colnames(mean.data) == "BC2.x")]
  mean.data <- mean.data[, !(colnames(mean.data) == "absorbance.x")]
  mean.data <- mean.data[, !(colnames(mean.data) == "BC.col.name.x")]
  mean.data <- mean.data[, !(colnames(mean.data) == "BC.col.x")]
  colnames(mean.data)[colnames(mean.data) == "NB.y"] <- "NB"
  colnames(mean.data)[colnames(mean.data) == "BC2.y"] <- "BC2"
  colnames(mean.data)[colnames(mean.data) == "absorbance.y"] <- "absorbance"
  colnames(mean.data)[colnames(mean.data) == "BC.col.y"] <- "BC.col"
  colnames(mean.data)[colnames(mean.data) == "BC.col.name.y"] <- "BC.col.name"
  
  mean.data <- mean.data[mean.data$outlier != 1, ]
  mean.data <-
    mean.data[!duplicated(paste(mean.data[, "level"],
                                mean.data[, type])), ]
  
  mean.data <- mean.data[mean.data[, backgroundval] != 1, ]
  
  #  mean.pl.data <- aggregate(cbind(absorbance, NB,
  #                                 BC2) ~ new.type + level + plate,
  #                            FUN = mean, data = data)
  
  #  mean.pl.data <- (merge(data, as.data.frame(mean.pl.data),
  #                         by = c("new.type", "level", "plate"), all = TRUE))
  
  
  #  mean.pl.data <- mean.pl.data[, !(colnames(mean.pl.data) == "NB.x")]
  #  mean.pl.data <- mean.pl.data[, !(colnames(mean.pl.data) == "BC2.x")]
  #  mean.pl.data <- mean.pl.data[, !(colnames(mean.pl.data) == "absorbance.x")]
  
  
  #  colnames(mean.pl.data)[colnames(mean.pl.data) == "NB.y"] <- "NB"
  #  colnames(mean.pl.data)[colnames(mean.pl.data) == "BC2.y"] <- "BC2"
  #  colnames(mean.pl.data)[colnames(mean.pl.data) == "absorbance.y"] <- "absorbance"
  
  #  mean.pl.data <-
  #    mean.pl.data[!duplicated(paste(mean.pl.data[, "level"],
  #                                   mean.pl.data[, type],
  #                                   mean.pl.data[, "plate"]   )), ]
  
  # mean.pl.data <- mean.pl.data[mean.pl.data[, backgroundval] != 1, ]
  
  A.data$data$bc.data <- data
  A.data$data$bc.mean <- mean.data
  #  A.data$data$bc.pl.mean <- mean.pl.data
  
  A.data$fits$bgModel <- fits
  if(!prior)
    A.data$auxiliary$btype = fit$Btype
  
  A.data$auxiliary$record <- record
  #A.data$absorb.list <- prob.list 
  
  A.data$auxiliary$shiny.input <- shiny.input
  
  #########################################
  
  A.data$call$bgModel <- this.call
  
  A.data$call$record <- callNumbering("bgModel", record = A.data$call$record)
  
  old <- A.data
  class(old) <- c("bgModel", class(old) %w/o% "bgModel")
  
  if(save & !drug.color.correction)
    save(old, file = data.file)
  if(!progressbar == "none" & length(new.platesets) > 0)
  close(pb)
  return(old)
}





##
##
## Function for correction of drug colour
##
##
drugColorCorrection <- 
  function(A.data = A.data,
           weights = "fitted", 
           fitted.a = FALSE, # used to be FALSE
           contr = c("sum", "helmert", "treatment"),
           outlier.test = 3,
           outlier.iter = 2,
           varpower.min = 10e-4,
           varpower.iter = 50,
           parametrisation = "unrestricted",
           update = TRUE,
           save  = TRUE,
           progressbar    = "text",
           shiny.input = NULL,
           session     = NULL)
  {
    
    ##############################################
    ##
    ##  The call is saved
    ##
    ###############################################
    
    ###############################
    ## create call with formals included
    
    this.call <- createCall(call = match.call(), "drugColorCorrection")
    
    ##############################################
    ##
    ##  The call is updated
    ##
    ###############################################
    
    call <- A.data$auxiliary$passed.var
    B              <- call$backgroundval
    idvar          <- call$idvar
    namevar        <- call$namevar
    drugvar        <- call$drugvar
    dosevar        <- call$dosevar 
    timevar        <- call$timevar
    data.file      <- call$data.file
    controlval     <- call$controlval
    additivevar    <- call$additivevar
    correctionname <- call$correctionname
    
    record  <- A.data$auxiliary$record
    
    ################################################
    ##
    ## make sure the right extension is used
    ##
    ################################################
    w.ext <- substr(basename(data.file), 
                    start=nchar(basename(data.file))- 5, 
                    nchar(basename(data.file))) == ".RData"
    
    if(!w.ext)
      data.file <- paste(data.file, ".RData", sep = "")
    this.call$data.file <- data.file
    
    #################################################
    prior <- FALSE
    if("drug.color.correct" %in% names(A.data))
      prior <- TRUE
    
    
    #if(prior){
    #  raw.data.all <- A.data$data$raw.data.all    
    #}else{
    raw.data.all <- A.data$data$raw.data
    if(prior)
      if(A.data$call$record["readDBFData", "last.visited"] < 
           A.data$call$record["drugColorCorrection", "last.visited"])
        raw.data.all <- A.data$data$raw.data.all
    
    #}
    ##########################
    ##
    ##  The correction data is established 
    ##
    ###########################
    if(any(raw.data.all[, namevar] == correctionname)){
      if(prior){
        correction.data <-  A.data$drug.color.correct$correction.data
        correction.data$auxiliary$record <- 
          record[record[, namevar] == correctionname, ]
        correction.data$auxiliary$passed.var <-
          A.data$auxiliary$passed.var
        colnames(correction.data$auxiliary$record)[
          colnames(correction.data$auxiliary$record) == "drugColorCorrection"] <- "bgModel"
        
        correction.data$data$raw.data <- 
          raw.data.all[raw.data.all[, namevar] == correctionname, ]
        
        correction.data$call <- A.data$call
        
        class(correction.data) <- class(A.data)
        class(correction.data) <- 
          c("correction.data", class(correction.data) %w/o% "correction.data")
      }else{
        correction.data <- list()
        
        correction.data$auxiliary$record <- 
          record[record[, namevar] == correctionname, ]
        
        correction.data$data$raw.data <- 
          raw.data.all[raw.data.all[, namevar] == correctionname, ]
        
        correction.data$call <- A.data$call
        
        correction.data$auxiliary$passed.var <-
          A.data$auxiliary$passed.var
        
        class(correction.data) <- class(A.data)
        class(correction.data) <- 
          c("correction.data", class(correction.data) %w/o% "correction.data")
      }
      ##########################
      ##
      ##  The correction data is normalised 
      ##
      ###########################
      
      correction.data <- do.call(
        bgModel, list(A.data = quote(correction.data), 
                      weights = eval(weights), 
                      fitted.a = eval(fitted.a), # used to be FALSE
                      contr = eval(contr),
                      outlier.test = eval(outlier.test),
                      outlier.iter = eval(outlier.iter),
                      varpower.min = eval(varpower.min),
                      varpower.iter = eval(varpower.iter),
                      parametrisation = eval(parametrisation),
                      progressbar     = eval(progressbar),
                      update = eval(update)))
      
      
      #       
      #      correction.data <- eval(substitute(
      #         bgModel(A.data = A.data.var,
      #                 weights = weights.var,
      #                 fitted.a = fitted.a.var, # used to be FALSE
      #                 contr = contr.var,
      #                 outlier.test = outlier.test.var,
      #                 outlier.iter = outlier.iter.var,
      #                 varpower.min = varpower.min.var,
      #                 varpower.iter = varpower.iter.var,
      #                 parametrisation = parametrisation.var,
      #                 update = update.var),
      #         list(A.data.var = correction.data,
      #                 weights.var = weights,
      #                 fitted.a.var = fitted.a, # used to be FALSE
      #                 contr.var = contr,
      #                 outlier.test.var = outlier.test,
      #                 outlier.iter.var = outlier.iter,
      #                 varpower.min.var = varpower.min,
      #                 varpower.iter.var = varpower.iter,
      #                 parametrisation.var = parametrisation,
      #                 update.var = update)))
      
      
      #  eval(substitute(glm(formula, poisson(), data, weights=1-w),
      #            list(w=as.name(zname))))
      
      ##########################
      ##
      ##  Amount of colour correction is calculated 
      ##
      ###########################
      
      background.list <- list()
      
      data <- correction.data$data$bc.data
      
      for(drugColor in unique(data[, drugvar])){
        data.color <- data[data[, timevar] == 0 & data[, B] != 1 &
                             data[, drugvar] == drugColor, ] 
        data.color$up <- data.color$BC + 1.96 * data.color$Std.Error
        data.color$lo <- data.color$BC - 1.96 * data.color$Std.Error
        
        formula <- as.formula(paste("cbind(BC, lo, up ) ~ ", dosevar))
        col.C <-
          aggregate(formula, FUN = mean,
                    data = data.color)
        
        col.C$BC2 <- col.C$BC -
          mean(data.color[data.color[, additivevar] == controlval , 
                          "BC"], na.rm = TRUE)
        
        col.C$Conc <- (col.C[, dosevar])
         col.C$Conc[col.C$Conc == 0] <-
           min(col.C$Conc[col.C$Conc != 0])/2
         col.C$Conc <- concConvert(col.C$Conc, 
                                   mol.mass = A.data$auxiliary$mol.data[drugColor, 2], 
                                   logfun.to = "log10")  
        
        background.list[[drugColor]] <- col.C
      }
      
      ##########################
      ##
      ##  The color correction data is used on raw data 
      ##
      ###########################
      
      raw.data <- 
        raw.data.all[raw.data.all[, namevar] != correctionname, ]
      data.cl <- correction.data$data$bc.data
      for(drugColor in unique(data.cl[, drugvar])){
        col.C <- background.list[[drugColor]]
        if("absorbance.nc" %in% colnames(raw.data))
          raw.data[raw.data[, drugvar] == drugColor, "absorbance"][
            !is.na(raw.data[raw.data[, drugvar] == drugColor, "absorbance.nc"])] <- 
          raw.data[raw.data[, drugvar] == drugColor, "absorbance.nc"][
            !is.na(raw.data[raw.data[, drugvar] == drugColor, "absorbance.nc"])]
        
        raw.data[raw.data[, drugvar] == drugColor, "absorbance.nc"] <- 
          raw.data[raw.data[, drugvar] == drugColor, "absorbance"]
        
        conc <- raw.data[raw.data[,drugvar] == drugColor, dosevar]
        rownames(col.C) <- col.C[,dosevar]
        
        raw.data[raw.data[, drugvar] == drugColor,   "absorbance"] <- 
          raw.data[raw.data[, drugvar] == drugColor, "absorbance"] -
          col.C[as.character(conc), "BC2"]
      }
      
      
      drug.color.correct <- list(background.list  = background.list, 
                                 correction.data  = correction.data)
    }else{
      drug.color.correct <- list(background.list  = NULL, 
                                 correction.data  = NULL)
    }
    
    
    
    # edit the record call 
    if(any(raw.data.all[, namevar] == correctionname)){
      if(! "drugColorCorrection" %in% colnames(record))
        record[, "drugColorCorrection"] <- 
        as.POSIXct("2011-04-04 14:18:58", tz="GB")
      
      record[is.na(record[, "drugColorCorrection"]), "drugColorCorrection"] <- 
        as.POSIXct("2011-04-04 14:18:58", tz="GB")
      
      record2 <- correction.data$auxiliary$record
      record2 <- record2[!duplicated(record2[, drugvar]),]
      rownames(record2) <- record2[,drugvar]
      for(drugColor in record2[, drugvar])
        record[record[,drugvar] == drugColor, "drugColorCorrection"] <- 
        record2[drugColor, "bgModel"]
    }else{
      if(! "drugColorCorrection" %in% colnames(record))
        record[, "drugColorCorrection"] <- 
        as.POSIXct("2011-04-04 14:18:58", tz="GB")
      
      record[is.na(record[, "drugColorCorrection"]), "drugColorCorrection"] <- 
        as.POSIXct("2011-04-04 14:18:58", tz="GB")
    }
    
    rem.drug <- names(A.data$drug.color.correct$background.list) %w/o%
      names(drug.color.correct$background.list)
    
    record[record[,drugvar] %in% rem.drug, "drugColorCorrection"] <- Sys.time()
    A.data$auxiliary$record <- record 
    A.data$drug.color.correct <- drug.color.correct
    if(any(raw.data.all[, namevar] == correctionname)){
      A.data$data$raw.data     <- raw.data
      A.data$data$raw.data.all <- raw.data.all  
    }else{
      A.data$data$raw.data <- raw.data.all
      A.data$data$raw.data.all <- raw.data.all  
    }
    ##########################################
    ##
    ##  updating the call
    A.data$call$record <- 
      callNumbering("drugColorCorrection", record = A.data$call$record)
    A.data$call$drugColorCorrection <- this.call
    
    A.data$auxiliary$shiny.input <- shiny.input
    
    ## giving a class
    class(A.data) <- c("drugColorCorrection", class(A.data) 
                       %w/o% "drugColorCorrection")
    
    ## saving it as old
    old <- A.data
    if(save)
      save(old,   file = data.file)
    return(old)
  }

##
##
## Function for creating bootstrapped datasets
##
##

bootstrap <- function(A.data, update = TRUE, n.samples = 50, max.iter = 100,
                      type = c("parametric", "residual", "nonparametric"),
                      progressbar = "text",
                      verbose = FALSE,
                      save = TRUE,
                      shiny.input = NULL,
                      session = NULL){
  
  type <- type[1]
  
  ###############################
  ## create call with formals included
  
  this.call <- createCall(call = match.call(), "bootstrap")
  
  ##########################################
  ##
  ##  The call is updated 
  ##
  ##########################################
  
  call <- A.data$auxiliary$passed.var
  
  data.file   <- eval(call$data.file)
  namevar     <- eval(call$namevar)
  controlval  <- eval(call$controlval)
  drugvar     <- eval(call$drugvar)
  timevar     <- eval(call$timevar)
  additivevar <- eval(call$additivevar)
  dosevar     <- eval(call$dosevar)
  
  BC <- "BC"
  ##########################################
  ##
  ##  The call to bgModel is used to fit the data
  ##
  ##########################################
  data <- A.data  # 
  call <- data$call$bgModel$call
  
  call$update        <- FALSE
  call$outlier.test  <- FALSE
  call$save          <- FALSE
  call$verbose       <- verbose
  call$progressbar       <- "none"
  
  for(form in names(formals(bgModel)))   
  if(is.null(call[[form]]) & !is.null(formals(bgModel)[[form]]))
    call[[form]] <- formals(bgModel)[[form]]

  
  record <- A.data$auxiliary$record
  idvar  <- eval(data$call$readDBFData$idvar)
  
  #call$data.file     <- file.path(data$auxiliary$passed.var$dbf.path, 
  #                                "boostrapSample.bg.RData")
  
  #############################################
  errors   <- 0
  
  prior <- all("bs.raw.data" %in% names(data$data), update)
  
  if(!prior){
    ## bs.data contains all the data including replicates
    bs.raw.data <- data$data$bc.data
    
    ## bs.data2 do not include outliers
    ## bs.raw.data2 <- data$data$bc.data[data$data$bc.data$outlier == 0, ]
    
    bootstrap.levels <- data.frame(level = unique(A.data$data$raw.data[, "level"]), n.samples = 0)
    bootstrap.levels$level <- as.character(bootstrap.levels$level)
    rownames(bootstrap.levels) <- bootstrap.levels$level
    
    record$bootstrap <- Sys.time()
    
    new.platesets <- bootstrap.levels$level
    #names(data$fits$bgModel[rep(1:4, length(data$fits$bgModel)/4) == 1])
  }
  
  if(prior){
    new.platesets <-  record[record$bootstrap < record$bgModel, idvar]
    new.levels    <- unique(A.data$data$raw.data[
      A.data$data$raw.data[, idvar] %in% new.platesets, "level"])
    
    ## removed compared to old
    rem.platesets <- unique(A.data$data$bs.raw.data[, idvar]) %w/o% 
      unique(A.data$data$raw.data[, idvar])
    rem.levels    <- unique(A.data$data$bs.raw.data[
      A.data$data$bs.raw.data[, idvar] %in% rem.platesets, "level"])
    
    level.relevant <- unique(A.data$data$bs.raw.data[A.data$data$bs.raw.data$level %in% rem.levels, 
                                                     "plateset"]) %w/o% rem.platesets
    
    
    ## if an replicate is removed a new background correction takes place.
    data.2      <- A.data$data$bs.raw.data[!(A.data$data$bs.raw.data[, idvar] %in% rem.platesets), ]
    rem.levels2 <- intersect(rem.levels, data.2$level)
    complete.rem.levels <- setdiff(rem.levels, rem.levels2)
    new.levels  <- c(new.levels, rem.levels2)
    
    
    #  if(length(new.platesets) > 0){
    #    data.new <- A.data$data$bc.data[A.data$data$bc.data[,idvar] %in% new.platesets, ]    
    #    int.col <- intersect(colnames(data.2), colnames(data.new))
    #    data.2[(nrow(data.2)+1):(nrow(data.2)+nrow(data.new)), int.col] <- data.new[, int.col]
    #  }
    bs.raw.data <- data.2
    
    rows <- rownames(A.data$data$bc.data) %w/o% rownames(bs.raw.data)
    if(length(rows) > 0)
      bs.raw.data[rows, colnames(A.data$data$bc.data)] <- 
      A.data$data$bc.data[rows, ]
    
    change.platesets <- c(new.platesets, rem.platesets )    
    change.levels <- c(new.levels, rem.levels)
    
    
    if(length(complete.rem.levels) > 0){
      bootstrap.levels <- 
        A.data$auxiliary$bootstrap.levels[!(A.data$auxiliary$bootstrap.levels$level 
                                            %in% complete.rem.levels),]    
    }else{
      bootstrap.levels <- A.data$auxiliary$bootstrap.levels
    }
    
    bootstrap.levels2 <- new.levels %w/o% bootstrap.levels$level
    if(length(new.levels) > 0){
      bootstrap.levels2 <- data.frame(level = new.levels, n.samples = 0)
      row.names(bootstrap.levels2) <- bootstrap.levels2$level
      
      bootstrap.levels <- rbind(bootstrap.levels, bootstrap.levels2)
    }
    
    new.platesets <- new.levels
    if(length(new.platesets) > 0)
      bootstrap.levels[new.levels,  "n.samples"] <- 0
    
    bootstrap.levels <- bootstrap.levels[rownames(bootstrap.levels) == bootstrap.levels[, "level"],]
  }
  
  # rownames(bs.raw.data) <- paste(bs.raw.data$idvar, bs.raw.data$level, 
  #                                 bs.raw.data$new.type, bs.raw.data$wellnum)
  
  bs.iter <- min(bootstrap.levels$n.samples)
  bs.iter2 <- bs.iter
  
  
  max <- sum(n.samples - bootstrap.levels$n.samples)
  if(!progressbar == "none" & max > 0)
    pb <- progressBar(title = paste("Bootstrapping", n.samples, "samples"), min = 0,
                      max = max, 
                      width = 300, 
                      window = progressbar == "window")
  
  count <- 0
  while(bs.iter < n.samples & bs.iter2 < max.iter) {
    
    bs.iter2 <-  bs.iter2 +1
    bs.iter <- min(bootstrap.levels$n.samples) + 1
    a <- Sys.time()
    new.platesets <- bootstrap.levels[bootstrap.levels$n.samples == bs.iter -1, "level"]
    bs.raw.data2 <- data$data$bc.data[data$data$bc.data$outlier == 0 &
                                        data$data$bc.data$level %in% new.platesets , ]
    
    for(i in new.platesets) {
      
      power <- data$fits$bgModel[[i]]$modelStruct$varStruct[1]
      
      names <- names(resid(data$fits$bgModel[[i]]))
      if(is.null(power)){
        power <- rep(1, length(names))
      }
      
      if(type[1] == "residual"){
        res <- resid(data$fits$bgModel[[i]])
        
        sc.res <- sample(sample( res /
                                   fitted(data$fits$bgModel[[i]])^power))
        
        bs.raw.data2[names ,"absorbance"] <-
          fitted(data$fits$bgModel[[i]]) + sc.res * fitted(data$fits$bgModel[[i]])^power
      }
      
      if(type[1] == "parametric" | type[1] == "nonparametric" ){
        sigma <- data$fits$bgModel[[i]]$sigma
        
        bs.raw.data2[names ,"absorbance"] <-
          fitted(data$fits$bgModel[[i]]) +
          fitted(data$fits$bgModel[[i]])^power * rnorm(length(names), sd = sigma)
      }
      count <- count + 1
      if(!progressbar == "none" & max > 0)
        setProgressBar(pb, count, label=paste( round(count/max*100, 2),
                                               "% done"),
                       window = progressbar == "window")
    }
    # not in for loop
   
  #  if(type[1] == "nonparametric"){
  #    bs.raw.data2[,"absorbance"] <-
  #      BootstrapPlateFun(X = bs.data, Y = bs.raw.data2)
  #  }
    
    #rownames(bs.raw.data2) <- paste(bs.raw.data2[,idvar], bs.raw.data2$level, 
    #                             bs.raw.data2$new.type, bs.raw.data2$wellnum)
    bs <- "error"
    A.data$data$raw.data <- bs.raw.data2[bs.raw.data2$outlier == 0 , ]
    
    
    #What ever the name the of the absorbance data used  it is replaced here
    eval(parse(text = paste(call$A.data, "<- A.data")))
    call$shiny.input <- NULL
  
  
  ##############
  ##

  
#   update<- FALSE
#   
#   parametrisation<-"unrestricted"
#   
#   outlier.test<-FALSE
#   
#   outlier.iter<-2
#   
#   weights<- "fitted"
#   
#   fitted.a<- FALSE
#   
#   progressbar<- "none"
#   
#   save<- FALSE
#   
#   verbose<- FALSE
#   
#   varpower.min<-0.001
#   
#   varpower.iter<- 50
#   
#   contr <-c("sum", "helmert", "treatment")
  
  ##
  ####################
  
  
  
  
  
  
  
  
    (bs <- do.call("bgModel", as.list(call)))
    
    xxx <- 
      bs$data$bc.data[bs$data$bc.data[,additivevar] == controlval, 
                      c("level", idvar, namevar, additivevar, timevar, BC, drugvar)]
    
    xxx$level1 <- paste(xxx[, namevar], xxx[,drugvar])
    xxx <- xxx[!duplicated(paste(xxx[, drugvar], xxx[, namevar], xxx[, timevar])), ]
    id.check <- vector()
    for(id.xxx in unique(xxx[, "level1"])){ # id istedet for level
      xx2 <- xxx[id.xxx  == xxx[, "level1"], ]
      if(length(unique(xx2[, timevar])) > 1){
        id.check <- c(id.check, xx2[which.max(xx2[, timevar]), BC] - 
                        xx2[which.min(xx2[, timevar]), BC] > 0)
      }else{
        id.check <- c(id.check, TRUE)
      }
    }
    checked <- xxx[ ,"level"][xxx[ ,idvar] %in% unique(xxx[ ,idvar]) [id.check]]
    
    if(class(bs)[1]=="character"){
      errors <- errors + 1
      if(verbose)
        cat("number of errors =", errors)
    }else{      
      bs.raw.data[rownames(bs$data$bc.data), 
                  paste("BS:", bs.iter, sep = "")] <- bs$data$bc.data[, BC]
      
      if(verbose){
        cat("bs iter", bs.iter, "\n")
        print(table(bs$data2$para))
        b <- Sys.time()
        cat("time left:", as.list((b-a) * (n.samples-bs.iter))[[1]] ,
            attributes(b-a)$units, "\n")
      }
      #bootstrap.levels[new.platesets, "n.samples"] <-
      #bootstrap.levels[new.platesets, "n.samples"] + 1
      bootstrap.levels[checked, "n.samples"] <-
        bootstrap.levels[checked, "n.samples"] + 1
      
    }
    ids <-  unique(A.data$data$raw.data[A.data$data$raw.data$level %in% new.platesets, idvar])
    if(length(ids) > 0)
      record[ids, "bootstrap"] <- Sys.time()
  }
  
  data$data$bs.raw.data <- bs.raw.data
  
  
  bs.raw.data <- bs.raw.data[bs.raw.data[,  "new.type"] !=  data$auxiliary$btype &
                               bs.raw.data$outlier != 1, ]
  bs.raw.data <- bs.raw.data[!duplicated(paste(bs.raw.data$level,
                                               bs.raw.data$new.type)), ]
  
  
  ## bs.mean contains mean data i.e. the value corresponding
  ## to A in the model d*A + B
  
  data$data$bs.mean <- bs.raw.data
  data$auxiliary$bootstrap.levels <- bootstrap.levels
  
  data$auxiliary$record <- record
  
  data$call$bootstrap <- this.call
  
  data$call$record <- callNumbering("bootstrap", data$call$record)
  
  old <- data  
  class(old) <- c("bootstrap", class(old) %w/o% "bootstrap")
  
  
  w.ext <- substr(basename(data.file), 
                  start=nchar(basename(data.file))- 5, 
                  nchar(basename(data.file))) == ".RData"
  if(!w.ext)
    data.file <- paste(data.file, ".RData", sep = "")
  
  old$auxiliary$shiny.input <- shiny.input
  if(save)
  save(old, file = data.file)
  if(!progressbar == "none" & max > 0)
  close(pb)
  return(old)
}

##
##
## Function for estimating the G model
##
##

growthModel <- function(A.data, 
                        parametrisation =  "restricted",
                        cut = 0.025, 
                        update = TRUE,
                        gnlscontrol = NULL,
                        save = TRUE,
                        progressbar = "text",
                        verbose = TRUE){
  
  # if(!("bgModel" %in% class(data)))
  #   stop("You need to perform a background correctin using bgModel")
  ## update the call according to read data
  
  
  ###############################
  ## create call with formals included
  
  this.call <- createCall(call = match.call(), "growthModel")
  
  
  # if(!("bgModel" %in% class(data)))
  #   stop("You need to perform a background correctin using bgModel")
  ## update the call according to read data
  
  # print("bootstrap")
  A.data$call$growthModel <- this.call
  
  
  #print("bootstrap")
  
  
  call   <- as.list(A.data$auxiliary$passed.var)
  controlval <- call$controlval
  # bs.call <- data$call$bootstrap
  A         <- "absorbance"
  B         <- eval(call$backgroundval)
  idvar     <- call$idvar
  drugvar   <- eval(call$drugvar)
  namevar   <- eval(call$namevar)
  timevar   <- eval(call$timevar)
  dosevar <- eval(call$dosevar)
  backgroundval <- eval(call$backgroundval)
  incubationvar <- eval(call$incubationvar)
  additivevar <- eval(call$additivevar)
  plateset <- "plateset"
  plate  <- "plate"
  type   <- "new.type"
  control <- A.data$data$bc.mean[A.data$data$bc.mean[, additivevar] == controlval, type][1]
  names(data)
  BC <- "BC"
  record <- A.data$auxiliary$record
  if("bootstrap" %in% class(A.data)){
    n.samples <- A.data$call$bootstrap$n.samples 
  }else{
    n.samples <- 0 
  }
  
  
  
  # time.points <- list()
  
  #if(!("all" %in% time.points))
  #  time.points$all <- sort(unique(data$data$bc.mean[,timevar]))
  
  prior <- all("GM.mean" %in% names(A.data$data), update)
  
  if(!prior){
    
    if("bootstrap" %in% class(A.data)){
      GM.mean <- A.data$data$bs.mean
      
      GM.names <- c(BC, paste("BS:", 1:n.samples, sep = ""))
      
      GM.mean[, paste("G", GM.names, sep = ".")] <- NaN  
    }else{
      GM.mean <- A.data$data$bc.mean
      
      GM.names <- c(BC)
      
      GM.mean[, paste("G", GM.names, sep = ".")] <- NaN  
    }    
    
    T0.list <- list()
    
    for(drug.iter in unique(GM.mean[, drugvar])){
      T0.mat <- matrix(NaN, ncol = length(GM.names),
                       nrow = length(unique(GM.mean[GM.mean[, drugvar] == drug.iter, namevar])))
      
      colnames(T0.mat) <- GM.names
      rownames(T0.mat) <- unique(GM.mean[GM.mean[, drugvar] == drug.iter, namevar])
      
      T0.mat <- as.data.frame(T0.mat)
      T0.list[[drug.iter]] <- T0.mat
    }
    
    record$GModel <- Sys.time()
    
  }
  if(prior){
    GM.mean <- A.data$data$GM.mean
    if("bootstrap" %in% class(A.data)){
      GM.names <- c(BC, paste("BS:", 1:n.samples, sep = ""))
    }else{
      GM.names <- c(BC)
    }
    
    
    new.platesets <-  record[record$GModel < record$bootstrap, idvar]
    new.levels    <- unique(A.data$data$raw.data[
      A.data$data$raw.data[, idvar] %in% new.platesets, "level"])
    
    
    ## removed compared to old    
    if("bootstrap" %in% class(A.data))
      rows <- intersect(rownames(GM.mean), rownames(A.data$data$bs.mean))
    
    if(!"bootstrap" %in% class(A.data))
      rows <- intersect(rownames(GM.mean), rownames(A.data$data$bc.mean))
    
    
    GM.mean <- GM.mean[rows, ]   
    
    
    ## add ekstra columns
    cols <- paste("G", GM.names, sep = ".") %w/o% colnames(GM.mean)
    GM.mean[, cols] <- NaN  
    
    cols <- GM.names %w/o% colnames(GM.mean)
    GM.mean[, cols] <- NaN  
    if("bootstrap" %in% class(A.data))
      GM.mean[rows, GM.names] <- A.data$data$bs.mean[rows, GM.names]
    
    if(!"bootstrap" %in% class(A.data))
      GM.mean[rows, GM.names] <- A.data$data$bc.mean[rows, GM.names]
    
    ## insert new
    rows <- rownames(A.data$data$bs.mean) %w/o% rownames(GM.mean)
    if("bootstrap" %in% class(A.data))
      if(length(rows) > 0)
        GM.mean[rows, colnames(A.data$data$bs.mean)] <- 
      A.data$data$bs.mean[rows, ]  
    if(!"bootstrap" %in% class(A.data))
      if(length(rows) > 0)
        GM.mean[rows, colnames(A.data$data$bc.mean)] <- 
      A.data$data$bc.mean[rows, ]
    
    
    ## levels with changes are removed
    GM.mean[GM.mean$level %in% new.levels, 
            paste("G", GM.names, sep = ".")] <- NaN
    
    T0.list <- A.data$data$T0.list
    
    for(drug.iter in unique(GM.mean[, drugvar])){
      if(!(drug.iter %in% names(T0.list))){
        T0.mat <- matrix(NaN, ncol = length(GM.names),
                         nrow = length(unique(GM.mean[
                           GM.mean[, drugvar] == drug.iter, namevar])))
        
        colnames(T0.mat) <- GM.names
        rownames(T0.mat) <- unique(GM.mean[GM.mean[, drugvar] == drug.iter, namevar])
        
        T0.mat <- as.data.frame(T0.mat)
        T0.list[[drug.iter]] <- T0.mat
      }
      T0.mat <- T0.list[[drug.iter]][unique(GM.mean[GM.mean[, drugvar] == drug.iter, namevar]), ]
    }
  }
  
  if(prior){
    xx <- GM.mean[GM.mean[, additivevar] != controlval, c(namevar, drugvar,
                                                          paste("G", GM.names, sep = "."))]
    wh.na.col <- apply(xx[,-c(1:2)], 2, function(x)  any(is.na(x)))
    wh.na.row <- apply(xx[,-c(1:2)], 1, function(x)  any(is.na(x)))
    
    GM.names <- GM.names[wh.na.col]
    if(length(GM.names > 0)){
      drug.iters <- unique(xx[wh.na.row, drugvar])      
      name.iters <- unique(xx[wh.na.row, namevar])
    }
  }else{
    drug.iters <- unique(GM.mean[, drugvar])
    name.iters <- unique(GM.mean[, namevar])
  }
  
  if(prior){
    fit.list <- A.data$fits$growthModel
  }else{
    fit.list <- list()
  }
  
  
  if(!progressbar == "none"){
    max <- 0
    for(bs.iter in GM.names){
      for(drug.iter in drug.iters){ 
        
        h.namevar <- unique(GM.mean[GM.mean[, drugvar] ==  drug.iter, namevar])
        name.iters.2 <-  h.namevar[h.namevar %in% name.iters]
        for(name.iter in name.iters.2){
          max <- max +1
        }}}
    
    pb <- progressBar(title = paste("Calculating the G model"), 
                      min = 0,
                      max = max, 
                      width = 300, 
                      window = progressbar == "window")
  }
  count <- 0
  for(bs.iter in GM.names){
    for(drug.iter in drug.iters){ 
      
      T0.mat <- T0.list[[drug.iter]]
      h.namevar <- unique(GM.mean[GM.mean[, drugvar] ==  drug.iter, namevar])
      name.iters.2 <-  h.namevar[h.namevar %in% name.iters]
      for(name.iter in name.iters.2){
        
        count <- count + 1
        if(!progressbar == "none")
          setProgressBar(pb, count, label=paste( round(count/max*100, 2),
                                                 "% done"),
                         window = progressbar == "window")
        
        ## name.iter = unique(bs.mean[bs.mean[, drugvar] == drug.iter, namevar])[23]
        name.temp.all <- GM.mean[GM.mean[, drugvar] == drug.iter &
                                   GM.mean[, namevar] == name.iter, ]
        # bs.mean <- data$data$bs.mean
        #name.temp <- bs.mean[bs.mean[, drugvar] == drug.iter &
        #                    bs.mean[, namevar] == name.iter, ]
        
        if(all(is.na(name.temp.all[, paste("G", bs.iter, sep = ".")])) &
             length(unique(name.temp.all[, timevar])) > 1){
          
          record[record[,namevar] == name.iter & record[,drugvar] == 
                   drug.iter, "GModel"] <- Sys.time()
          
          
          name.temp.all <- name.temp.all[name.temp.all[, type] != paste(0, backgroundval), ]
          #name.temp.all <- name.temp.all[name.temp.all[, "outlier"] != 1, ]
          
          #if(!is.null(time.points))
          time.points <- sort(unique(name.temp.all[, timevar])) %w/o% 0
          for(time.iter in 0:length(time.points)){
            if(time.iter != 0){
              name.temp <- name.temp.all[name.temp.all[, timevar] %in% 
                                           c(0, time.points[time.iter]), ]
            }else{
              name.temp <- name.temp.all
            }
            
            if(!(time.iter == 1 & length(time.points) == 1)){
              
              logLik <- 0
              if(parametrisation == "restricted"){
                if(verbose){
                print("------------------------------------------------------")
                print(paste(bs.iter, drug.iter, name.iter, time.iter))
                }# n.start <- ifelse(bs.iter == "BC", n.start.values, n.start.values.bs)
                
                suppressWarnings(
                  fit.abs <- growth.function(data            = name.temp, 
                                             cut             = cut,
                                             parametrisation = "abs",
                                             type            = type,
                                             timevar         = timevar,
                                             dosevar         = dosevar,
                                             incubationvar   = incubationvar,
                                             additivevar     = additivevar,
                                             controlval      = controlval,
                                             A               = bs.iter,
                                             #n.start.values  = n.start,
                                             gnlscontrol     = gnlscontrol)
                  
                )
                suppressWarnings(
                  fit.square <- growth.function(data            = name.temp, 
                                                cut             = cut,
                                                parametrisation = "square",
                                                type            = type,
                                                timevar         = timevar,
                                                dosevar         = dosevar,
                                                incubationvar   = incubationvar,
                                                additivevar     = additivevar,
                                                controlval      = controlval,
                                                A               = bs.iter,
                                                #n.start.values  = n.start,
                                                gnlscontrol     = gnlscontrol)
                )
                #dev.off()
                #                  par(mfrow = c(1,2))
                #                 ylim <- range(unlist(lapply(fit.abs$fit.list, function(x) x$logLik)),
                #                               unlist(lapply(fit.square$fit.list, function(x) x$logLik)))
                #                   plot(unlist(lapply(fit.abs$fit.list, function(x) x$logLik)),
                #                        main = "abs", ylim = ylim)
                #                 plot(unlist(lapply(fit.square$fit.list, function(x) x$logLik)), 
                #                      main = "Square", ylim = ylim)
                logLik <- c(fit.abs$fit$logLik, fit.square$fit$logLik )
                names(logLik) <- c("abs", "square")
                if(verbose){
                print(which.max(logLik))
                print(logLik)
                }# print(paste("square: n.fits = ", length(fit.square$fit.list)))
                # print(paste("abs: n.fits = ", length(fit.abs$fit.list)))
                if(max(logLik) > 0 ) {
                  if(names(which.max(logLik)) == "abs"){
                    #if(logLik["square"]  == 0){
                    fit <- fit.abs
                  }else{
                    fit <- fit.square
                  }
                }
                
              }else{
                if(verbose)
                print(paste(bs.iter, drug.iter, name.iter, time.iter))
                fit <-
                  growth.function(data = name.temp, cut = cut,
                                  parametrisation = "reci",
                                  type = type,
                                  timevar = timevar,
                                  dosevar = dosevar,
                                  incubationvar = incubationvar,
                                  additivevar = additivevar,
                                  controlval  = controlval,
                                  A = bs.iter)
                logLik <- fit$fitlogLik
              }            
            }
            if(class(fit$fit)[1] != "error"){
              if(time.iter == 0){
                name.temp2 <- 
                  GM.mean[GM.mean[, drugvar] == drug.iter &
                            GM.mean[, namevar] == name.iter, ]
                GM.mean[GM.mean[, drugvar] == drug.iter &
                          GM.mean[, namevar] == name.iter, 
                        paste("G", bs.iter, sep =".")] <-
                  fit[["summary"]][name.temp2[, type], "G"]
                
                GM.mean[GM.mean[, drugvar] == drug.iter &
                          GM.mean[, namevar] == name.iter, 
                        paste("T0", bs.iter, sep ="")] <-
                  fit[["summary"]][1, "T0"]          
                
                T0.mat[name.iter, bs.iter] <-  fit[["summary"]][1, "T0"]
              }else{
                name.temp2 <- 
                  GM.mean[GM.mean[, drugvar] == drug.iter &
                            GM.mean[, namevar] == name.iter&
                            GM.mean[, timevar] == time.points[time.iter], ]
                
                GM.mean[GM.mean[, drugvar] == drug.iter &
                          GM.mean[, namevar] == name.iter &
                          GM.mean[, timevar] == time.points[time.iter], 
                        paste("Gres", bs.iter, sep =".")] <-
                  fit[["summary"]][name.temp2[, type], "G"]
                
                GM.mean[GM.mean[, drugvar] == drug.iter &
                          GM.mean[, namevar] == name.iter &
                          GM.mean[, timevar] == time.points[time.iter], 
                        paste("T0res", bs.iter, sep =".")] <-
                  fit[["summary"]][1, "T0"]  
              }             
            }
            a.time <- ifelse(time.iter == 0, "all", time.points[time.iter])
            fit.list[[drug.iter]][[name.iter]][[bs.iter]][[paste(a.time)]] <- fit            
          }       
        }
      }   
      T0.list[[drug.iter]] <-  T0.mat
    }     
  }
  
  
  GM.mean[, paste("G.upper", GM.names, sep = ".")] <- 
    GM.mean[, paste("G", GM.names, sep = ".")]
  
  GM.mean[, paste("Gres.upper", GM.names, sep = ".")] <- 
    GM.mean[, paste("Gres", GM.names, sep = ".")]
  
  
  GM.mean[, paste("G.lower", GM.names, sep = ".")] <- 
    GM.mean[, paste("G", GM.names, sep = ".")] / 100 /
    GM.mean[, paste("T0", GM.names, sep = "")]
  
  GM.mean[, paste("Gres.lower", GM.names, sep = ".")] <- 
    GM.mean[, paste("Gres", GM.names, sep = ".")] / 100 /
    GM.mean[, paste("T0res", GM.names, sep = ".")]
  
  GM.mean <- GM.mean[GM.mean[, "new.type"] != control & GM.mean[, timevar] != 0, ]
  
  A.data$data$GM.mean <- GM.mean
  A.data$data$T0.list <- T0.list 
  A.data$auxiliary$record <- record
  A.data$fits$growthModel <- fit.list
  A.data$call$growthModel <- this.call
  
  
  A.data$call$record <- callNumbering("growthModel", A.data$call$record)
  
  #  A.data$GM.time.points <- time.points
  class(A.data) <- c("growthModel", class(A.data) %w/o% "growthModel")  
  
  old <- A.data
  data.file <- eval(A.data$auxiliary$passed.var$data.file)
  if(save)
  save(old, file = data.file)
  if(progressbar != "none" & max > 0)
  close(pb)
  return(old)
}



## Function for calculating the measures based on the old growth models

createGIData <- function(A.data, update = TRUE, cut = 0.025,
                         use.supplied.T0 = FALSE,
                         doublingvar = "supT0",
                         model  = c("R", "D", "RG", "DG"),
                         progressbar = "text",
                         save = TRUE,
                         shiny.input = NULL){
  
  
  
  if(update == FALSE)
    A.data$data$GI.mean <- NULL
  
  ###############################
  ## create call with formals included
  
  this.call <- createCall(call = match.call(), "createGIData")
  
  
  call        <- A.data$auxiliary$passed.var
  #bs.call    <- A.data$call$bootstrap
  A           <- "absorbance"
  B           <- eval(call$backgroundval)
  idvar       <- call$idvar
  drugvar     <- eval(call$drugvar)
  namevar     <- eval(call$namevar)
  timevar     <- eval(call$timevar)
  dosevar     <- eval(call$dosevar)
  additivevar <- eval(call$additivevar)
  controlval  <- eval(call$controlval)
  if(is.null(doublingvar))
    doublingvar  <- call$doublingvar
  plateset    <- "plateset"
  plate       <- "plate"
  type        <- "new.type"
  
  max <- sum(any(grepl("D", model)), any(grepl("R", model))) + 4
  if(progressbar != "none" & max > 0)
  pb <- progressBar(title = paste("Calculating the", 
                                  paste(model, collapse = ", "), 
                                  "models"), min = 0,
                    max = max, 
                    width = 300, 
                    window = progressbar == "window")
  
  
  if("bootstrap" %in% class(A.data)){
    bs.mean <- A.data$data$bs.mean
  }else{
    bs.mean <- A.data$data$bc.mean
  } 
  
  control <- bs.mean[bs.mean[, additivevar] == controlval, type][1]
  
  BC <- "BC"
  
  X <- A.data
  
  
  if("bootstrap" %in% class(A.data)){
    n.samples <- X$call$bootstrap$n.samples
    bs.names  <- c(BC, paste("BS:", 1:n.samples, sep = "")) 
  }else{
    n.samples <- NULL
    bs.names <- BC
  }
  prior <- all("GI.mean" %in% names(A.data), update)
  
  if(prior)
    prior <- all(model %in% eval(X$call$createGIdata$model))
  if(prior)
    prior <- all(bs.names %in% colnames(X$data$GI.mean))
  
  record <- X$auxiliary$record
  
  if(!prior){    
    record$GIModel <- Sys.time()
  }
  
  if(prior) {
    
    GI.mean  <- A.data$data$GI.mean
    
    if("bootstrap" %in% class(A.data)){
      new.platesets  <-  record[record$GIModel < record$bootstrap, idvar]
    }else{
      new.platesets  <-  record[record$GIModel < record$bgModel, idvar]
    }
    new.platesets2 <- unique(A.data$data$bs.mean[,idvar] %w/o% GI.mean[,idvar])
    new.platesets  <- unique(new.platesets, new.platesets2)
    new.levels     <- unique(A.data$data$raw.data[
      A.data$data$raw.data[, idvar] %in% new.platesets, "level"])
    
    
    
    ## removed compared to old    
    rows <- intersect(rownames(GI.mean), rownames(bs.mean))
    rows.nogo <- setdiff(rownames(GI.mean), rownames(bs.mean))
    GI.mean <- GI.mean[rows, ]   
    
    ## insert new
    #    rows <- rownames(A.data$data$bs.mean) %w/o% rownames(GI.mean)
    #    GI.mean[rows, colnames(A.data$data$bs.mean)] <- 
    #            A.data$data$bs.mean[rows, ]
    
    #X$data$bs.mean <- X$data$bs.mean[X$data$bs.mean$level %in% new.levels, ]
    bs.mean <- bs.mean[bs.mean$level %in% new.levels, ]
    if(length(new.platesets) > 0)
      record[new.platesets, "GIModel"] <- Sys.time()
    
    if(nrow( bs.mean)) {
      mat <- bs.mean[, bs.names] < cut
      for(i in bs.names)
        bs.mean[mat[,i], i] <- cut
    }
  }else{
    GI.mean <- bs.mean
    mat <- bs.mean[, bs.names, drop=FALSE] < cut
    mat[is.na(mat)] <- FALSE 
    for(i in bs.names)
      bs.mean[mat[,i,drop = FALSE], i] <- cut   
  }
  
  count <- 0
  suppressWarnings(
    if(nrow(bs.mean)) {
      bs.mean$entity <- paste(bs.mean[,namevar], bs.mean[, drugvar])
      X$data$bs.mean <- bs.mean
      X$data$bc.mean <- bs.mean
      
      
      X <- A0T1fun(X, timevar = timevar, namevar = namevar, control = control, 
                   entity = "entity", BC = BC)
      
      count <- count + 1
      if(!progressbar == "none")
        setProgressBar(pb, count, label=paste( round(count/max*100, 2),
                                               "% done"),
                       window = progressbar == "window")
      
      X <- AC0fun(X,  timevar = timevar, control = control, 
                  entity = "entity", BC = BC)
      
      
      count <- count + 1
      if(!progressbar == "none")
        setProgressBar(pb, count, label=paste( round(count/max*100, 2),
                                               "% done"),
                       window = progressbar == "window")
      
      X <- doublingTime(X,  timevar = timevar, namevar = namevar, 
                        use.supplied.T0 = use.supplied.T0, doublingvar = doublingvar,
                        control = control, entity = "entity",  BC = BC)
      
      
      count <- count + 1
      if(!progressbar == "none")
        setProgressBar(pb, count, label=paste( round(count/max*100, 2),
                                               "% done"),
                       window = progressbar == "window")
      
      X$data$GI.mean <- X$data$GI.mean[X$data$GI.mean[, "new.type"] 
                                       != control & X$data$GI.mean[, timevar] != 0, ]
      
      count <- count + 1
      if(!progressbar == "none")
        setProgressBar(pb, count, label=paste( round(count/max*100, 2),
                                               "% done"),
                       window = progressbar == "window")
      
      
      if("R" %in% model | "RG" %in% model ){
        X <- RGrowth(X = X, BC = BC, update = update, n.samples = n.samples,
                     timevar = timevar, doublingvar = doublingvar, 
                     RG = eval("RG" %in% model))
        
        count <- count + 1
        if(!progressbar == "none")
          setProgressBar(pb, count, label=paste( round(count/max*100, 2),
                                                 "% done"),
                         window = progressbar == "window")
        
      }
      if("D" %in% model | "DG" %in% model){
        X <- DGrowth(X = X, BC = BC, update = update, 
                     timevar = timevar, 
                     DG = eval("DG" %in% model))
        
        count <- count + 1
        if(!progressbar == "none")
          setProgressBar(pb, count, label=paste( round(count/max*100, 2),
                                                 "% done"),
                         window = progressbar == "window")
      }
      
      if(prior)
        A.data$data$GI.mean <- rbind(GI.mean, X$data$GI.mean)
      if(!prior)
        A.data$data$GI.mean <- X$data$GI.mean
      A.data$auxiliary$record  <- record
      
    }else{
      A.data$data$GI.mean <- GI.mean
    }
    
  )
  
  A.data$call$createGIData <- this.call
  
  A.data$call$record <- callNumbering("createGIData", A.data$call$record)
  
  old <- A.data
  class(old) <- c("createGIData", class(old) %w/o% "createGIData")
  
  data.file <- eval(old$call$readDBFData$data.file)
  if(save)
  save(old, file = data.file)
  if(progressbar != "none")
  close(pb)
  return(old)
}




##
##
##     Combine data objects created by R and D model 
##     with that created using the G model
##
##

combineGIdata <- function(A.data, save = TRUE){
  
  ###############################
  ## create call with formals included
  
  this.call <- createCall(call = match.call(), "combineGIdata")
  #################################################################
  
  call <- A.data$auxiliary$passed.var
  
  namevar <- eval(call$namevar)
  drugvar <- eval(call$drugvar)
  timevar <- eval(call$timevar)
  
  dosevar  <- eval(call$dosevar)
  additivevar <-  eval(call$additivevar)
  
  ##################################################################
  
  if(!is.null(A.data$data$GI.mean) & !is.null(A.data$data$GM.mean)){
  xx <- intersect(names(A.data$data$GI.mean), names(A.data$data$GM.mean)) %w/o% 
    c(namevar, drugvar, timevar, dosevar)
  xx <- colnames(A.data$data$GM.mean) %w/o% xx
  
  A.data$data$DR.data <- merge(A.data$data$GI.mean, A.data$data$GM.mean[, xx],
                               by = c(namevar, drugvar, timevar, dosevar))
  
  }else{
    if(!is.null(A.data$data$GI.mean))
      A.data$data$DR.data <- A.data$data$GI.mean
    if(!is.null(A.data$data$GM.mean))
      A.data$data$DR.data <- A.data$data$GM.mean
  }
  A.data$call$combineGIdata <- this.call
  A.data$call$record <- callNumbering("combineGIdata", A.data$call$record)
  
  ## the pass along list is updated
  
  for(input.iter in names(as.list(this.call$call))[-1]){
    A.data$auxiliary$passed.var[[input.iter]] <-
      this.call[[input.iter]] 
  }
  
  old <- A.data
  
  class(old) <- c("DRdata", class(old) %w/o% "DRdata")
  data.file <- eval(old$call$readDBFData$data.file)
  if(save)
    save(old, file = data.file)
  return(old)  
}




isoreg.DRdata <-
  function(A.data, 
           logfun.from = "nolog",
           logfun.to  = "log10",
           verbose = FALSE,
           save = TRUE,
           progressbar = "text"){
    
    
    ###############################
    ## create call with formals included
    
    this.call <- createCall(call = match.call(), "isoreg.DRdata")
    
    #############################################
    ##
    ##   The variables are passed along from the previous 
    ##   function calls
    ##
    #############################################
    
    ## The variables are passed along from the previous 
    ## function calls
    call        <- A.data$auxiliary$passed.var
    namevar     <- eval(call$namevar)
    dosevar     <- eval(call$dosevar)
    drugvar     <- eval(call$drugvar)
    timevar     <- eval(call$timevar)
    controlval  <- eval(call$controlval)
    additivevar <- eval(call$additivevar)
    
    BC <- "BC"    
    ################################################
    ##
    ## The summary lists are created
    ##
    ###############################################
    
    ## The number of bootstrap samples
    if("bootstrap" %in% class(A.data)){
      n.samples <- eval(A.data$call$bootstrap$n.samples)
      bs.names <- 
        c("BC", paste("BS:", 1:n.samples, sep = ""))
      bs.mean <- A.data$data$bs.mean
    }else{
      bs.names <- "BC"
      bs.mean <- A.data$data$bc.mean
    }
    ## Which models shoul be calculated
    
    models <- c()
    if(paste("D.upper", BC, sep = ".") %in% colnames(A.data$data$DR.data))
      models <- c(models, "D")
    if(paste("DG.upper", BC, sep = ".") %in% colnames(A.data$data$DR.data))
      models <- c(models, "DG")
    if(paste("R", BC, sep = ".") %in% colnames(A.data$data$DR.data))
      models <- c(models, "R")
    if(paste("RG.upper", BC, sep = ".") %in% colnames(A.data$data$DR.data))
      models <- c(models, "RG")
    if(paste("G.upper", BC, sep = ".") %in% colnames(A.data$data$DR.data))
      models <- c(models, "G")
    if(paste("Gres.upper", BC, sep = ".") %in% colnames(A.data$data$DR.data))
      models <- c(models, "Gres")
    
    
    # For each drug available in split.data
    drugs <- unique(A.data$data$DR.data[, drugvar])      
    iso.fit <- list()
    range <- list()
    
    
    
    
    ##########################
    ## Calculate the number of doseresponse curves that needs to be fit
    max <- length(bs.names) * length(models) *
      length(unique(apply(A.data$data$DR.data[
        A.data$data$DR.data[, timevar] != 0, c(drugvar,timevar, namevar)], 1, 
                          function(x) paste(x, collapse = " "))))
    
    if("G" %in% models){
      max <- max - length(bs.names) * 
        length(unique(apply(A.data$data$DR.data[
          A.data$data$DR.data[, timevar] != 0, 
          c(drugvar,timevar, namevar)], 1, 
                            function(x) paste(x, collapse = " "))))
      
      
      max <- max +  length(bs.names) * length(unique(
        apply(A.data$data$DR.data[, c(drugvar, namevar)], 1, 
              function(x) paste(x, collapse = " "))))
    }
    if(progressbar != "none")
    pb <- progressBar(title = paste("Estimating isotonic regression"), 
                      min = 0,
                      max = max, 
                      width = 300, 
                      window = progressbar == "window")
   
    count <- 0
    for(drug.iter in drugs){
      drug.rang <- data.frame(min = NA, max = NA)
      data.drug.all <- 
        A.data$data$DR.data[A.data$data$DR.data[, drugvar] == drug.iter, ]
      
      data.drug.all[,dosevar] <- 
        switch(logfun.from[1],
               nolog    = data.drug.all[,dosevar],
               log10 = 10^data.drug.all[,dosevar],
               log2  =  2^data.drug.all[,dosevar],
               log   = exp(data.drug.all[,dosevar]))
      
      if(!(logfun.to[1] == "nolog")){
        data.drug.all[,dosevar] <- 
          do.call(logfun.to[1], list(x =  data.drug.all[,dosevar]))
      }
      for(m in models){
        
        times <- sort(unique(data.drug.all[, timevar]))
        
        if(m == c("G"))
          times <- "all"
        
        data.list <- list()
        
        for(time.iter in times){
          if(m == "G"){
            data.time <- data.drug.all[data.drug.all[, timevar] ==
                                         data.drug.all[, timevar][1]  , ]
          }else{
            data.time <-  data.drug.all[data.drug.all[, timevar] == time.iter, ]
          }
          
          names <- unique(data.time[,  namevar])
          
          
          for(name.iter in names){
            if(m != "G"){
              data.name <- A.data$meta.list$readDBFData[
                A.data$meta.list$readDBFData[, namevar] == name.iter &
                  A.data$meta.list$readDBFData[, drugvar] == drug.iter &
                  A.data$meta.list$readDBFData[, timevar] == time.iter ,][
                    1,, drop = FALSE]
            }else{
              data.name <- A.data$meta.list$readDBFData[
                A.data$meta.list$readDBFData[, namevar] == name.iter &
                  A.data$meta.list$readDBFData[, drugvar] == drug.iter,][
                    1,, drop = FALSE]
            }
            m.drug <- list()
            for(bs.iter in bs.names){
              count <- count + 1
              if(!progressbar == "none")
                setProgressBar(pb, count, label=paste( round(count/max*100, 2),
                                                       "% done"),
                               window = progressbar == "window")
              
              
              x <-  data.time[, dosevar][data.time[, namevar] == name.iter]
              
              drug.rang[nrow(drug.rang) + 1, ] <- range(x)
              
              if(m == "R"){
                y <-  data.time[data.time[,namevar] == name.iter,
                                paste(m, bs.iter, sep = "."),]
                data.line <- findOneLine(x, upper= y)
              }else{                   
                y.l <-  data.time[data.time[, namevar] == name.iter,
                                  paste(m, "lower", bs.iter, sep = ".")]
                y.h <-  data.time[data.time[, namevar] == name.iter,
                                  paste(m, "upper", bs.iter, sep = ".")]
                if(!all(is.na(y.l)) & !all(is.na(y.h))){
                  data.line <- findOneLine(x, upper= y.h, lower = y.l)
                }else{
                  data.line <- data.frame(x = x, fitted = NA, orig = NA, which = NA)
                    data.line <- data.line[order(data.line$x), c("x", "fitted", "orig", "which")]
                }
              }
              
              colnames(data.line) <- c(dosevar, paste(m, "iso", sep = "."), 
                                       m, "origin")
              
              m.drug[[length(m.drug) + 1]] <- 
                cbind(data.line, data = bs.iter)
              
            } # bs.iter
            suppressWarnings(
              data.list[[length(data.list) + 1]] <- 
                cbind(data.name, rbind.fill(m.drug))
            )
            
          } # name .iter
        } # time.iter
        iso.fit[[drug.iter]][[m]] <- rbind.fill(data.list)
        
        
        
        
      } # model.iter
      range[[drug.iter]] <- drug.rang[-1,]
    } # drug.iter
    A.data$auxiliary$drug.range <- 
      lapply(range, function(x) c(max(x[,1]), min(x[,2])))
    A.data$call$isoreg.DRdata <- this.call
    A.data$data$iso.fits <- iso.fit
    A.data$call$record <- callNumbering("isoreg.DRdata", A.data$call$record)
    class(A.data) <- c("isoreg.DRdata", class(A.data) %w/o% "isoreg.DRdata" )
    if(progressbar != "none")
    close(pb)
    return(A.data)
  }




doseResponseModel <- function(
  A.data,
  models           = c("G", "R", "D", "DG", "RG"),  
  dose.scale       = "mol/l",
  dose.logfun.from = "nolog",
  dose.logfun.to   = "log10",
  parametrisation  = "restricted",  
  AUC.q            = 0,
  t                = 48,
  doublingvar      = NULL,
  cut = 0.025, 
  update = TRUE,
  verbose = FALSE,
  progressbar    = "text",
  save = TRUE,
  shiny.input = NULL,
  session     = NULL
  ){
  
  if(!any(models != "G"))
    models <- unique(c(models, "R"))
  if("G" %in% models){
    if(verbose| progressbar == "text")
      cat("Calculating the G-model \n\n")
    
    
    A.data <- do.call(growthModel,
                      list(A.data = quote(A.data),
                           parametrisation = eval(parametrisation),
                           cut = eval(cut),
                           update = eval(update),
                           verbose = eval(verbose),
                           progressbar = progressbar),
                      )
    
  }
  
  if(any(models != "G")){
    if(verbose | progressbar == "text" )
      cat("Calculating the models:", models %w/o% "G", "\n\n")
    
    A.data <- do.call(createGIData,
                      list(A.data = quote(A.data),
                           cut = eval(cut),
                           doublingvar = eval(doublingvar),
                           model = eval(models %w/o% "G"),
                           update = eval(update),
                           progressbar = progressbar
                      ))   
  }
  
  A.data <- combineGIdata(A.data)
  
  if(verbose| progressbar == "text" )
    cat("Isotonic regression:\n\n")
  
  
  A.data <- do.call(isoreg.DRdata,
                    list(A.data      =  quote(A.data),
                         logfun.from =  dose.logfun.from, 
                         logfun.to   =  dose.logfun.to,
                    progressbar = progressbar))
  
  
  
  if(verbose| progressbar == "text" )
    cat("calculating summary statistics GI50, TGI, LC50, and AUC.\n\n")
  
  A.data <- do.call(summary.DRdata,
                    list(A.data = quote(A.data),
                         t = eval(t), 
                         AUC.q  = eval(AUC.q),
                         type.fit= "iso",
                         dose.scale = eval(dose.scale), 
                         dose.logfun= eval(dose.logfun.to ),
                         progressbar = progressbar))
  
  
  ############################################
  
  A.data$call[["doseResponseModel"]] <- createCall(match.call(), "doseResponseModel")
  
  A.data$call$record <- callNumbering("doseResponseModel", A.data$call$record)
  
  class(A.data)  <- c("doseResponseModel", class(A.data) %w/o% "doseResponseModel")
  
  A.data$auxiliary$shiny.input <- shiny.input
  
  old <- A.data
  if(save)
  save(old, file = A.data$auxiliary$passed.var$data.file)
  return(A.data)  
}


changeMolarMass <- function(A.data, drug, mol.mass, save = TRUE, shiny.input = NULL){
  
  A.data$auxiliary$mol.data[drug, "mol.mass"] <- mol.mass
  
  if(!is.null(shiny.input))
    A.data$auxiliary$shiny.input <- shiny.input
  
  old <- A.data
  if(save)
    save(old, file = A.data$auxiliary$passed.var$data.file)
  return(A.data)  
}



summary.DRdata <-
  function(A.data, 
           t = 48,
           AUC.q = 0, #"GI50", #
           type.fit = "iso",
           dose.scale = "mol/l",
           dose.logfun = "log10", 
           verbose = TRUE,
           progressbar = "text"){
    
    ###############################
    ## create call with formals included
    
    this.call <- createCall(call = match.call(), "summary.DRdata")
    
    
    
    #############################################
    ##
    ##   The variables are passed along from the previous 
    ##   function calls
    ##
    #############################################
    
    ## The variables are passed along from the previous 
    ## function calls
    call       <- A.data$auxiliary$passed.var
    namevar    <- eval(call$namevar)
    dosevar    <- eval(call$dosevar)
    drugvar    <- eval(call$drugvar)
    timevar    <- eval(call$timevar)
    controlval <- eval(call$controlval)
    additivevar <- eval(call$additivevar)
    
    BC <- "BC"#A.data$call$bootstrap$BC
    
    ################################################
    ##
    ## The summary lists are created
    ##
    ###############################################
    
    summary <- list()
    
    GI50.l <- list()
    TGI.l  <- list()
    LC50.l <- list()
    AUC.l  <- list()
    
    ## The GI50, TGI, and LC50 have different values 
    GIov <- matrix(c(50, 50, 75, 50, 50, 50,
                     0,  0, 50,  0, 0, 0,
                     -50, -50, 25, -50,  -1/t, -1/t),
                   ncol = 6, nrow = 3, byrow = TRUE)
    rownames(GIov) <- c("GI50", "TGI", "LC50")
    colnames(GIov) <- c("D", "DG", "R", "RG", "G", "Gres")
    
    
    ## The number of bootstrap samples
    if("bootstrap" %in% class(A.data)){
      n.samples <- eval(A.data$call$bootstrap$n.samples)
      bs.names <- 
        c("BC", paste("BS:", 1:n.samples, sep = ""))
      bs.mean <- A.data$data$bs.mean
    }else{
      bs.names <- "BC"
      bs.mean <- A.data$data$bc.mean
    }
    ## Which models shoul be calculated
    
    models <- c()
    if(paste("D.upper", BC, sep = ".") %in% colnames(A.data$data$DR.data))
      models <- c(models, "D")
    if(paste("DG.upper", BC, sep = ".") %in% colnames(A.data$data$DR.data))
      models <- c(models, "DG")
    if(paste("R", BC, sep = ".") %in% colnames(A.data$data$DR.data))
      models <- c(models, "R")
    if(paste("RG.upper", BC, sep = ".") %in% colnames(A.data$data$DR.data))
      models <- c(models, "RG")
    if(paste("G.upper", BC, sep = ".") %in% colnames(A.data$data$DR.data))
      models <- c(models, "G")
    if(paste("Gres.upper", BC, sep = ".") %in% colnames(A.data$data$DR.data))
      models <- c(models, "Gres")
    
    
    
    summary <- list()
    # For each drug available in split.data
    drugs <- unique(A.data$data$DR.data[,drugvar])      
    
    ##########################
    ## Calculate the number of doseresponse curves that needs to be fit
    max <- length(bs.names) * length(models) *
      length(unique(apply(A.data$data$DR.data[
        A.data$data$DR.data[, timevar] != 0, c(drugvar,timevar, namevar)], 1, 
                          function(x) paste(x, collapse = " "))))
    
    if("G" %in% models){
      max <- max - length(bs.names) * 
        length(unique(apply(A.data$data$DR.data[
          A.data$data$DR.data[, timevar] != 0, 
          c(drugvar,timevar, namevar)], 1, 
                            function(x) paste(x, collapse = " "))))
      
      
      max <- max + length(bs.names) * length(unique(
        apply(A.data$data$DR.data[, c(drugvar, namevar)], 1, 
              function(x) paste(x, collapse = " "))))
    }
    if(progressbar != "none" & max > 0)
    pb <- progressBar(title = paste("Estimating summary statistics"), 
                      min = 0,
                      max = max, 
                      width = 300, 
                      window = progressbar == "window")
    count <- 0
    
    for(drug.iter in drugs){
      data.drug.all <- A.data$data$DR.data[
        A.data$data$DR.data[, drugvar] == drug.iter, ]
      for(m in models){
        
        
        
        ################################################
        ##  The concentration is recalculated to match the call
        ##
        ###############################################
        range <- A.data$auxiliary$drug.range[[drug.iter]]
        
        iso.fit <- A.data$data$iso.fits[[drug.iter]][[m]]
        
        #iso.fit <- iso.fit[iso.fit[, dosevar] >= range[1] &
        #                     iso.fit[, dosevar] <= range[2]  ,]
        
        mol <- A.data$auxiliary$mol.data[drug.iter, "mol.mass"]
        
        iso.fit[, dosevar] <- 
          concConvert(x = iso.fit[, dosevar], mol.mass= mol, 
                      from = iso.fit$unit, 
                      to   = dose.scale,
                      logfun.from = eval(A.data$call$isoreg.DRdata$logfun.to),
                      logfun.to = dose.logfun)
        
        ##### time points are found
        times <- sort(unique(data.drug.all[, timevar]))
        if(m == "G")
          times <- "all"
        for(time.iter in times){
          if(m != "G"){
            names <- unique(iso.fit[iso.fit[,  timevar] == time.iter,  namevar])
          }else{
            names <- unique(iso.fit[,  namevar])
          }
          
          template <- matrix(NaN, ncol = length(names),
                             nrow = length(bs.names))
          
          colnames(template) <- names
          rownames(template) <- bs.names
          
          GI50 <- template
          TGI  <- template
          LC50 <- template
          AUC  <- template
          
          for(name.iter in names){
            for(bs.iter in bs.names){
              
              count <- count + 1
              if(progressbar != "none" & max > 0)
                setProgressBar(pb, count, label=paste( round(count/max*100, 2),
                                                       "% done"),
                               window = progressbar == "window")
              
              
              if(m != "G"){
                wh <- iso.fit[, namevar] == name.iter &
                  iso.fit[, timevar] == time.iter &
                  iso.fit[, "data"] == bs.iter
              }else{
                wh <- iso.fit[, namevar] == name.iter &
                  iso.fit[, "data"] == bs.iter
              }
              y <- iso.fit[wh, paste(m, type.fit = "iso", sep = ".")]
              x <- iso.fit[wh, dosevar] 
              
              AUC.q2 <- AUC.q
              if(AUC.q == "GI50")
                AUC.q2 <- GIov["GI50", m]
              if(AUC.q == "TGI")
                AUC.q2 <- GIov["TGI", m]
              
              if(!all(is.na(y))){
                GI50[bs.iter, name.iter] <- 
                  findGI(x, y, GIov["GI50", m]) # try
                TGI[bs.iter, name.iter] <- 
                  findGI(x, y, GIov["TGI",  m])  # try
                LC50[bs.iter, name.iter] <- 
                  findGI(x, y, GIov["LC50", m]) # try
                AUC[bs.iter, name.iter] <- findAUCq(x, y, q = AUC.q2, lim = range)  
              }                    
            }
          }
          if(m == "G"){
            summary[[drug.iter]][[m]]$GI50 <- GI50
            summary[[drug.iter]][[m]]$TGI  <- TGI
            summary[[drug.iter]][[m]]$LC50 <- LC50
            summary[[drug.iter]][[m]]$AUC  <- AUC
          }else{
            summary[[drug.iter]][[m]][["GI50"]][[paste(time.iter)]] <- GI50
            summary[[drug.iter]][[m]][["TGI"]][[paste(time.iter)]]  <- TGI
            summary[[drug.iter]][[m]][["LC50"]][[paste(time.iter)]] <- LC50
            summary[[drug.iter]][[m]][["AUC"]][[paste(time.iter)]]  <- AUC
          }
        } # time.iter
      } # model.iter
    } # drug.iter
    
    
    A.data$summary <- summary
    A.data$call$summary.DRdata <- this.call
    A.data$call$record <- callNumbering("summary.DRdata", A.data$call$record)
    class(A.data) <- c("summary.DRdata", class(A.data) %w/o% "summary.DRdata" )
    if(progressbar != "none" & max > 0)
    close(pb)
    return(A.data)
  }





getSummary <- function(A.data, drug = 1,
                       model = "G",
                       type = "AUC.iso",
                       time = 1){
  if(model == "G"){
    return(A.data$summary[[drug]][[model]][[type]])
  }else{
    return(A.data$summary[[drug]][[model]][[paste(time)]][[type]])
  }
  
}




CI <- function(A.data, ...) UseMethod("CI")

CI.summary.DRdata <- function(A.data,..., model = "G",
                              types =  c("GI50", "TGI", "LC50", "AUC"),
                              conf = c(2.5, 97.5),
                              time = "48",
                              dose.scale = "mol/l",
                              dose.logfun = "log10"){
  
  
  sum.list <- list()
  if(max(conf) > 1)
    conf <- conf/100
  
  for(drug in names(A.data$summary)){
    
    if(!is.null(A.data$data$T0.list)){
      T0.ci <- apply(A.data$data$T0.list[[drug]][, -1, drop = FALSE],
                     1, FUN = function(x)
                       quantile(x, c(conf[1], conf[2]), na.rm = TRUE ))
      
      T0 <- A.data$data$T0.list[[drug]][, 1, drop = FALSE]
      colnames(T0) <- "T0"
      T0[, paste("T0", conf*100, sep = ".")] <- t(T0.ci)
    }else{
      
      T0 <- data.frame()
    }
    for(type in types){
      if(model == "G"){
        GI <- A.data$summary[[drug]][[model]][[type]]
      }else{
        GI <- A.data$summary[[drug]][[model]][[type]][[time]]
      }
      if(!grepl("AUC", type))       
        GI <- concConvert(GI, mol.mass = A.data$auxiliary$mol.data[drug, "mol.mass"], 
                          logfun.from=A.data$call$summary.DRdata$dose.logfun   ,         
                          from = A.data$call$summary.DRdata$dose.scale, 
                          to = dose.scale,
                          logfun.to = dose.logfun)
      
      
      GI.ci <- t(apply(GI[-1, ,drop = FALSE],
                       2, FUN = function(x)
                         quantile(x, c(conf[1], conf[2]), na.rm = TRUE )))
      
     
      if(nrow(T0) == 0){
        T0 <- data.frame(E = rep(NA, ncol(GI)))
        rownames(T0) <- colnames(GI)
        colnames(T0) <- type
      }
      
      int <- intersect(rownames(T0), colnames(GI))
    
      T0[, type] <- NaN
      
      T0[int, type] <- GI[1, int]
      
      
      T0[, paste(type, conf*100, sep = "")] <- NaN
      int <- intersect(rownames(T0), rownames(GI.ci))
      T0[int, paste(type, conf*100, sep = "")] <- GI.ci[int, ]      
    }
    sum.list[[drug]] <- T0   
  }
  return(sum.list)
}


pdfCI <- 
  function(A.data, 
           model       = "G",
           types       = c("T0", "GI50", "TGI","LCt", "AUCq"),
           conf        = c(2.5, 97.5),
           dec         = 2,
           size        = "small",
           type.order  = "GI50",
           dose.scale  = "mol/l", 
           dose.logfun = "log10",
           drugs        = 1,
           splitvar    = NULL,
           file       = "",
           pdfheight   = "297mm",
           pdfwidth    = "210mm",
           keep.tex    = TRUE,
           compile.tex = TRUE,
           png         = FALSE,
           clean       = TRUE, ...) {
    
    if(is.numeric(drugs))
      drugs <- names(A.data$data$iso.fits)[drugs]
    
    tempfile <- base::tempfile(fileext = ".pdf")
    
    q <- A.data$call$doseResponseModel$AUC.q
    t <- A.data$call$doseResponseModel$t
    
    round2 <- function(x, dec = 2) {
      x.orig <- x
      x.vec <- vector()
      for(x in x.orig){
        x   <- round(x, dec)
        x2  <- strsplit(as.character(x), "\\.")
        nc  <- ifelse(is.na((x2[[1]][2])), 0, nchar(x2[[1]][2]))
        dif <- dec - nc
        n0  <-  paste(rep("0", dif), collapse  = "")
        if(is.na((x2[[1]][2]))){
          x.vec <- c(x.vec, paste(x2[[1]][1], ".", n0, collapse = "", sep = ""))
        }else{
          x.vec <- c(x.vec, paste(x2[[1]][1], ".",x2[[1]][2], n0, collapse = "", sep = ""))
        }
      }
      x.vec
    }
    
    sum.list <- CI(A.data, model= model, conf = conf,
                   dose.scale = dose.scale, 
                   dose.logfun = dose.logfun)
    
    types.av <- colnames(sum.list[[1]])
    types2 <- intersect(types, types.av)
    
    
    if("LCt" %in% types & "LC50" %in% types.av )
      types2 <- c(types2, "LCt")
    if("AUCq" %in% types & "AUC" %in% types.av )
      types2 <- c(types2, "AUCq")
    
    types <- types2
      
    names <- vector()
    
    for(drug in drugs)
      names <- c(names, rownames(sum.list[[drug]]))
    names <- unique(names)
    
    data <- data.frame(cell.line = names)
    row.names(data) <- data$cell.line
    for(drug in drugs){
      x1 <- sum.list[[drug]]
      
      xx <- data.frame(cell.line = names)
      row.names(xx) <- xx$cell.line
      
      if("T0" %in% types)
        xx[, "T0"] <- paste(round(x1$T0), " (", round(x1$T0.2.5), ";", 
                            round(x1$T0.97.5), ")", sep = "")
      
      if("GI50" %in% types)
        xx[, "GI50"] <- paste(round2(x1$GI50, dec), " (", round2(x1$GI502.5, dec), ";", 
                              round2(x1$GI5097.5, dec), ")", sep = "")
      
      if("TGI" %in% types)
        xx[, "TGI"] <-paste(round2(x1$TGI, dec), " (", round2(x1$TGI2.5, dec), ";", 
                            round2(x1$TGI97.5, dec), ")", sep = "")
      
      
      if("LCt" %in% types)
        xx[, "LCt"] <-paste(round2(x1$LC50, dec), " (", round2(x1$LC502.5, dec), ";", 
                             round2(x1$LC5097.5, dec), ")", sep = "")
      
      if("AUCq" %in% types)
        xx[, "AUCq"] <-paste(round2(x1$AUC, dec), " (", round2(x1$AUC2.5, dec), ";", 
                             round2(x1$AUC97.5, dec), ")", sep = "")
      
      xx <- xx[, -1]
      
      xx[xx == "NaN.00 (NA.00;NA.00)"] <- "-"
      xx["NA.00 (NA.00;NA.00)" == xx] <- "-"
#       xx <- data.frame(
#         T0 = paste(round(x1$T0), " (", round(x1$T0.2.5), ";", 
#                    round(x1$T0.97.5), ")", sep = ""),
#         GI50  = paste(round2(x1$GI50, dec), " (", round2(x1$GI502.5, dec), ";", 
#                       round2(x1$GI5097.5, dec), ")", sep = ""),
#         TGI = paste(round2(x1$TGI, dec), " (", round2(x1$TGI2.5, dec), ";", 
#                     round2(x1$TGI97.5, dec), ")", sep = ""),
#         LCt = paste(round2(x1$LC50, dec), " (", round2(x1$LC502.5, dec), ";", 
#                     round2(x1$LC5097.5, dec), ")", sep = ""),
#         AUCq = paste(round(x1$AUC), " (", round(x1$AUC2.5), ";", 
#                      round(x1$AUC97.5), ")", sep = "")
#         
#       )
#       
#       
#       row.names(xx) <- rownames(x1)
#       xx <- xx[, types]
#       
      xx <- xx[rownames(data), ]
      data <- cbind(data, xx)
      
    }
    
    if(!is.null(splitvar)){
      data[, splitvar] = A.data$meta.list$additional[row.names(data), splitvar]
      data <- data[order(data[,splitvar],data[,type.order] ,row.names(data)), ]
    }
    
    #require(tools)
    
    if(file != ""){
      sink(gsub(".pdf",".tex",tempfile))
      cat("\\documentclass{report}\n")
      cat("\\usepackage{booktabs}\n")
      cat("\\usepackage{lscape}\n")
      cat(paste("\\usepackage[centering,",
                "paperwidth=",   pdfwidth,
                ",paperheight=", pdfheight,
                ",noheadfoot,",
                "margin=0in]{geometry}\n", sep = ""))
      
      cat("\\begin{document}\\pagestyle{empty}\n")
    }
    
    
    if(length(types) > 1 & !is.null(splitvar) & length(drugs) == 1){
      
      colheads <-  
        c("Hours", rep(c("$GI_{50}$", "$TGI$" ,
                         paste('$LC_{',t,'}$',sep =''),
                         paste('$AUC_{',q,'}$',sep ='')), 1))
      
      names(colheads) <- c("T0", rep(c("GI50", "TGI","LCt", "AUCq"), 1))
      
      colheads <- colheads[types]
      if("T0" %in% types){ 
        cgroup = c("T0", paste("Model", model)) 
        n.cgroup =c(1, rep(length(types) - 1, 1))
      }else{
        cgroup = paste("Model", model)
        n.cgroup =rep(length(types) - 1, 1)
      }
      rownames(data) <- gsub("\\_", "\\\\textunderscore ", rownames(data))
      latex(data[, -c(1)], #booktabs = TRUE,, ncol(data)
            file = "",
            title = "Cell Line",
            colheads = colheads,
            cgroup = cgroup,
            n.cgroup = n.cgroup,
            rgroup = names(table(data[,splitvar])),
            n.rgroup = table(data[,splitvar]),,
            size = size,
            ...)
    }
    
    if(length(types) > 1 & is.null(splitvar) & length(drugs) == 1){
      
      colheads <-  
        c("Hours", rep(c("$GI_{50}$", "$TGI$" ,
                         paste('$LC_{',t,'}$',sep =''),
                         paste('$AUC_{',q,'}$',sep ='')), 1))
      
      names(colheads) <- c("T0", rep(c("GI50", "TGI","LCt", "AUCq"), 1))
      
      colheads <- colheads[types]
      if("T0" %in% types){ 
        cgroup = c("T0", paste("Model", model)) 
        n.cgroup =c(1, rep(length(types) - 1, 1))
      }else{
        cgroup = paste("Model", model)
        n.cgroup =rep(length(types) - 1, 1)
      }
      rownames(data) <- gsub("\\_", "\\\\textunderscore ", rownames(data))
      latex(data[, -c(1)], booktabs = TRUE,#, ncol(data)
            file = "",
            title = "Cell Line",
            colheads = colheads,
            cgroup = cgroup,
            n.cgroup = n.cgroup,
            size = size,
            ...)
    }
    
    if(length(types) > 1 & !is.null(splitvar) & length(drugs) > 1){
      
      colheads <-  
        rep(c("$T_0$", "$GI_{50}$", "$TGI$" ,
              paste('$LC_{',t,'}$',sep =''),
              paste('$AUC_{',q,'}$',sep ='')), length(drugs))
      
      names(colheads) <-  rep(c("T0", "GI50", "TGI","LCt", "AUCq"), length(drugs))
      
      colheads <- colheads[types]
      gsub("\\","\\textunderscore", data[, ])
      
      rownames(data) <- gsub("\\_", "\\\\textunderscore ", rownames(data))
      latex(data[, -c(1)], booktabs = TRUE, #, ncol(data)
            file = "",
            title = "Cell Line",
            colheads = colheads,
            cgroup = drugs,
            n.cgroup = rep(length(types), length(drugs)),
            rgroup = names(table(data[,splitvar])),
            n.rgroup = table(data[,splitvar]),
            size = size,
            ...)
    }
    if(file != ""){
      cat("\\end{document}\n")
      sink()
      if (compile.tex) {
        texi2dvi(file =  gsub(".pdf",".tex", tempfile),
                 pdf = TRUE, clean = clean)
        
        file.copy(basename(tempfile), file, overwrite = TRUE)
        file.remove(basename(tempfile))
        file.copy(gsub(".pdf",".tex",tempfile), gsub(".pdf", ".tex", file), 
                  overwrite = TRUE)
        if(png){
          #require(animation)
          file.png <- im.convert(file, output = gsub(".pdf", ".png", file),
                                 extra.opts="-density 300")
          file.png <- strsplit(file.png, paste(file, " ", sep = ""))[[1]][2]
          
          file.copy(eval(parse(text= file.png )) , gsub(".pdf", ".png", file), 
                    overwrite = TRUE)
        }           
      }
      
      
      
      if (!keep.tex) {
        file.remove(gsub(".pdf", ".tex", file))
      }
    }
  }


  tableCI <- 
  function(A.data, 
           model       = "G",
           dec         = 2,
           conf        = c(2.5, 97.5), 
           types       = c("T0", "GI50", "TGI","LCt", "AUCq"),
           type.order  = "GI50",
           dose.scale  = "mol/l", 
           dose.logfun = "log10",
           drugs        = 1,
           splitvar    = NULL
  ) {
    
    if(is.numeric(drugs))
      drugs <- names(A.data$data$iso.fits)[drugs]
    
    q <- A.data$call$doseResponseModel$AUC.q
    t <- A.data$call$doseResponseModel$t
    
    round2 <- function(x, dec = 2) {
      x.orig <- x
      x.vec <- vector()
      for(x in x.orig){
        x   <- round(x, dec)
        x2  <- strsplit(as.character(x), "\\.")
        nc  <- ifelse(is.na((x2[[1]][2])), 0, nchar(x2[[1]][2]))
        dif <- dec - nc
        n0  <-  paste(rep("0", dif), collapse  = "")
        if(is.na((x2[[1]][2]))){
          x.vec <- c(x.vec, paste(x2[[1]][1], ".", n0, collapse = "", sep = ""))
        }else{
          x.vec <- c(x.vec, paste(x2[[1]][1], ".",x2[[1]][2], n0, collapse = "", sep = ""))
        }
      }
      x.vec
    }
    
    sum.list <- CI(A.data, model= c("G"), conf = conf,
                   dose.scale = dose.scale, 
                   dose.logfun = dose.logfun)
    
    names <- vector()
    
    for(drug in drugs)
      names <- c(names, rownames(sum.list[[drug]]))
    names <- unique(names)
    
    data <- data.frame(cell.line = names)
    
    row.names(data) <- data$cell.line
    if(!is.null(splitvar)){
      data[, splitvar] = A.data$meta.list$additional[row.names(data), splitvar]      
    }
    for(drug in drugs){
      x1 <- sum.list[[drug]]
      
      xx <- data.frame(
        T0 = paste(round(x1$T0), " (", round(x1$T0.2.5), ";", 
                   round(x1$T0.97.5), ")", sep = ""),
        GI50  = paste(round2(x1$GI50, dec), " (", round2(x1$GI502.5, dec), ";", 
                      round2(x1$GI5097.5, dec), ")", sep = ""),
        TGI = paste(round2(x1$TGI, dec), " (", round2(x1$TGI2.5, dec), ";", 
                    round2(x1$TGI97.5, dec), ")", sep = ""),
        LCt = paste(round2(x1$LC50, dec), " (", round2(x1$LC502.5, dec), ";", 
                    round2(x1$LC5097.5, dec), ")", sep = ""),
        AUCq = paste(round(x1$AUC), " (", round(x1$AUC2.5), ";", 
                     round(x1$AUC97.5), ")", sep = "")
        
      )
      wh <- match(types, colnames(xx))
      
      colnames(xx) <- c("T0", "GI50", "TGI", 
                        paste("LC", t, sep = ""), 
                        paste("AUC", q, sep = ""))
      
      xx <- xx[, wh]
      row.names(xx) <- rownames(x1)
      
      
      xx <- xx[rownames(data), ]
      data <- cbind(data, xx)
      
    }
    data <- data[,-1]
    if(!is.null(splitvar)){     
      data <- data[order(data[,splitvar], data[,type.order], row.names(data)), ]
    }else{
      data <- data[order(data[,type.order], row.names(data)), ]
    }
    data
  }



DeltaPlot <- function(A.data, 
                      diff.var    = 1, 
                      drug        = 1, 
                      model       = "G", 
                      type        = "AUC", 
                      time        = "48",
                      names       = NULL,
                      dose.scale  = "mol/l",
                      dose.logfun = "log10",
                      col         = "#4F6E9F",
                      reverse     = FALSE, 
                      las         = 2, 
                      ylab        = NULL, 
                      main        = "Difference", 
                      ...){
  
  namevar <- A.data$auxiliary$passed.var$namevar
  correctionname <- A.data$auxiliary$passed.var$correctionname
  if(is.numeric(diff.var))
    diff.var <- A.data$auxiliary$passed.var$namecols[diff.var]
  
  if(is.numeric(drug))
    drug <- names(A.data$summary)[drug]
  
  q <- A.data$call$doseResponseModel$AUC.q
  t <- A.data$call$doseResponseModel$t
  
  if(model == "G"){
    if(is.null(ylab)){
      if(type == "GI50")
        ylab <- expression(paste(Delta, GI[50]))
      if(type == "AUC")
        ylab <- as.expression(substitute(paste(Delta, AUC[a]), list(a = q)))
      if(type == "TGI")
        ylab <- as.expression((paste(Delta, TGI)))
      if(type == "LC50")
        ylab <-  as.expression(substitute(paste(Delta, LC[a]), list(a = t)))
    }
  }
  if(model == "D"){
    if(is.null(ylab)){
      if(type == "GI50")
        ylab <- expression(paste(Delta, GI[50]))
      if(type == "AUC")
        ylab <- as.expression(substitute(paste(Delta, AUC[a]), list(a = q)))
      if(type == "TGI")
        ylab <- as.expression((paste(Delta, TGI)))
      if(type == "LC50")
        ylab <-  as.expression(substitute(paste(Delta, LC[a]), list(a = 50)))
    }
  }
  if(model == "R"){
    if(is.null(ylab)){
      if(type == "GI50")
        ylab <- expression(paste(Delta, GI[75]))
      if(type == "AUC")
        ylab <- as.expression(substitute(paste(Delta, AUC[a]), list(a = q)))
      if(type == "TGI")
        ylab <- expression(paste(Delta, GI[50]))
      if(type == "LC50")
        ylab <-  expression(paste(Delta, GI[25]))
    }
  }
  
  types <- unique(A.data$meta.list$metadata.long[, diff.var])
  
  
  if(is.na(types[1])){
    names.1 <- A.data$meta.list$metadata.long[is.na(A.data$meta.list$metadata.long[, diff.var]), namevar] %w/o% correctionname
  }else{
    names.1 <- A.data$meta.list$metadata.long[A.data$meta.list$metadata.long[, diff.var] == types[1] & 
                                                !is.na(A.data$meta.list$metadata.long[, diff.var]), namevar] %w/o% correctionname
  }
  
  if(is.na(types[2])){
    names.2  <- A.data$meta.list$metadata.long[is.na(A.data$meta.list$metadata.long[, diff.var]), namevar] %w/o% correctionname
  }else{
    names.2 <- A.data$meta.list$metadata.long[A.data$meta.list$metadata.long[, diff.var] == types[2] & 
                                                !is.na(A.data$meta.list$metadata.long[, diff.var]), namevar] %w/o% correctionname
  }
  
  names.1 <- unique(names.1)
  names.2 <- unique(names.2)
  if(!is.null(names)){
    
    grepl.v <- function(patterns, x){
      h1 <- vector()
      for(i in 1:length(x)){
        h2 <- vector()
        for(j in 1:length(patterns))
          h2[j] <- grepl(patterns[j], x[i])
        h1[i] <- any(h2)
      }
      return(h1)
    }
    
    names.1 <- names.1[grepl.v(names, names.1)]
    names.2 <- names.2[grepl.v(names, names.2)]
  }
  if(model == "G"){
    data <- A.data$summary[[drug]][[model]][[type]]
  }else{
    data <- A.data$summary[[drug]][[model]][[type]][[time]]
  }
  
  if(!grepl("AUC", type))       
    data <- DoseR:::concConvert(data, mol.mass = A.data$auxiliary$mol.data[drug, "mol.mass"], 
                                logfun.from=A.data$call$summary.DRdata$dose.logfun   ,         
                                from = A.data$call$summary.DRdata$dose.scale, 
                                to = dose.scale,
                                logfun.to = dose.logfun)
  
  
  data.1   <- data[, names.1]
  data.2   <- data[, names.2]
  
  data <- (data.1 - data.2 ) 
  if(reverse)
    data <- data * -1
  
  data <- data[, order(data[1,])]
  
  boxplot(data[-1, ,drop = FALSE], las = las, ylab = ylab, col =col,
          main = main, ...)
}


DeltaCI <- function(A.data, 
                    conf     = c(2.5, 97.5), 
                    diff.var = 1, 
                    drug     = 1, 
                    model    = "G", 
                    type     = "AUC", 
                    time     = "48",
                    dose.scale  = "mol/l",
                    dose.logfun = "log10",
                    reverse     = FALSE){
  
  if(max(conf) > 1)
    conf <- conf/100
  
  namevar <- A.data$auxiliary$passed.var$namevar
  correctionname <- A.data$auxiliary$passed.var$correctionname
  if(is.numeric(diff.var))
    diff.var <- A.data$auxiliary$passed.var$namecols[diff.var]
  
  if(is.numeric(drug))
    drug <- names(A.data$summary)[drug]
  
  
  types <- unique(A.data$meta.list$metadata.long[, diff.var])
  
  
  if(is.na(types[1])){
    names.1 <- A.data$meta.list$metadata.long[is.na(A.data$meta.list$metadata.long[, diff.var]), namevar] %w/o% correctionname
  }else{
    names.1 <- A.data$meta.list$metadata.long[A.data$meta.list$metadata.long[, diff.var] == types[1] & 
                                                !is.na(A.data$meta.list$metadata.long[, diff.var]), namevar] %w/o% correctionname
  }
  
  if(is.na(types[2])){
    names.2  <- A.data$meta.list$metadata.long[is.na(A.data$meta.list$metadata.long[, diff.var]), namevar] %w/o% correctionname
  }else{
    names.2 <- A.data$meta.list$metadata.long[A.data$meta.list$metadata.long[, diff.var] == types[2] & 
                                                !is.na(A.data$meta.list$metadata.long[, diff.var]), namevar] %w/o% correctionname
  }
  
  names.1 <- unique(names.1)
  names.2 <- unique(names.2)
  
  if(model == "G"){
    data <- A.data$summary[[drug]][[model]][[type]]
  }else{
    data <- A.data$summary[[drug]][[model]][[type]][[time]]
  }
  
  if(!grepl("AUC", type))       
    data <- DoseR:::concConvert(data, mol.mass = A.data$auxiliary$mol.data[drug, "mol.mass"], 
                                logfun.from=A.data$call$summary.DRdata$dose.logfun   ,         
                                from = A.data$call$summary.DRdata$dose.scale, 
                                to = dose.scale,
                                logfun.to = dose.logfun)
  
  
  
  data.1   <- data[, names.1]
  data.2   <- data[, names.2]
  
  data <- (data.1 - data.2 ) 
  if(reverse)
    data <- data * -1
  
  data <- data[, order(data[1,])]
  
  GI.ci <- t(apply(data[-1, ,drop = FALSE],
                   2, FUN = function(x)
                     quantile(x, c(conf[1], conf[2]), na.rm = TRUE ))) 
  
  return(cbind(data[1,], GI.ci))
}
