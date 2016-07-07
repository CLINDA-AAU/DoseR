
##
##
## Function for reading dbf data into R
##
##

readDBFDataShiny <- function(
  A.data,  # Absorbance, data.frame, list  
  update         = TRUE, 
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
  
  this.call$session <- NULL
  
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
  if(!is.null(shiny.input)){
    progressbar = "none"
  }
  withProgress(session, min = 0,
              max = nrow(data), {
                setProgress(message = 'Calculation in progress',
                            detail = 'This may take a while...')   
                
                 if(!progressbar == "none"){
                   on.exit(close(pb))
                   pb <- progressBar(title = "Reading dbf files", min = 0,
                                     max = nrow(data), width = 300, window = progressbar == "window" )
                 }
                 
                 #  pb <- DoseResponse:::progressBar(title = "Reading dbf files", min = 0,
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
                       protocols[[prot]][["setup"]] <- read.xls(filesetup, sheet = 1)
                       protocols[[prot]][["conc"]]  <- read.xls(filesetup, sheet = 2)
                     }else{
                       if(is64){
                         protocols[[prot]][["setup"]] <- 
                           gdata::read.xls(filesetup, sheet = 1,
                                           perl = "C:/perl/perl/bin/perl.exe")
                         protocols[[prot]][["conc"]]  <- 
                           gdata::read.xls(filesetup, sheet = 2,
                                           perl = "C:/perl/perl/bin/perl.exe")
                       }else{
                         protocols[[prot]][["setup"]] <- 
                           read.xls(filesetup, colNames = TRUE, sheet = 1)
                         protocols[[prot]][["conc"]]  <- 
                           read.xls(filesetup, colNames = TRUE, sheet = 2)
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
                   
                   setProgress(value = i, message="Reading the dBase files",
                               detail=paste( round(i/nrow(data)*100, 0),
                                           "% done"))
                   
                   
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
                     mol.data[drug.iter, ] <- 
                       c(drug.iter, info[[drug.iter]]$"Chemical data"["Molar mass", ])
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
                   row.names(new.meta) <- new.meta[,idvar]
                   
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
               })
              
}




##
##
##     Function for fitting the nonlinear function dA + b
##
##


bgModelShiny <- function(A.data, #..., 
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
drugColorCorrectionShiny <- 
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

bootstrapShiny <- function(A.data, update = TRUE, n.samples = 50, max.iter = 100,
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
    
    try(bs <- do.call("bgModel", as.list(call)))
    
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


doseResponseModelShiny <- function(
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
    models <- c(models, "R")
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

