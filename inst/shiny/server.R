
try(setwd(A.data))

"%w/o%" <- function(x, value) x[!x %in% value]



colscheme <- list()
colscheme$line <- list(col = c("#71965A", "#4F6E9F", "#9F9692", "#9D2441",
                               "#333333", "#662D91", "#71DEC0", "#F7931E"),
                       lty = 1, lwd = 1)

colscheme$point <- list(col = c("#71965A", "#4F6E9F", "#9F9692", "#9D2441",
                                "#333333", "#662D91", "#71DEC0", "#F7931E"),
                        pch = 1)

colscheme$grid <- list(col = "grey", lwd = par("lwd"), lty = "dotted")
colscheme$bar <-  list(col = "grey")
colscheme$bs <-   list(col = "grey", lty = 1, lwd = 0.5)

colscheme$boxplot <- list(col = NULL)


###########################################################
##  Shiny server
############################################################

A.data <<- list()
shinyServer(function(input, output, session) {
  
 # A.data <- list()
  
  current.protocol <- list()
  DBF.current.protocol <- list()
  ##########################################
  ## read an old A.data file
#   fileExt <- function(x) {
#     db <- grepl("\\.[^.]+\\.(gz|bz2|xz)$", x)
#     ans <- sub(".*\\.", "", x)
#     ans[db] <- sub(".*\\.([^.]+\\.)(gz|bz2|xz)$", "\\1\\2", 
#                    x[db])
#     ans
#   }
  
  observe({
    if(input$savefileButton == 0)
      return()
    data.file <- input$datafile
    if("createMetadata" %in% class(A.data) | "Absorbance" %in% class(A.data)){
      
      if(data.file != DoseR:::no.extension(A.data$auxiliary$passed.var$data.file)){
        
        xls.file <- paste(DoseR:::no.extension(A.data$auxiliary$passed.var$data.file), ".xls", sep = "")
        excel.2 <- paste(data.file, ".xls", sep = "")
        file.copy(from=xls.file, to=excel.2)
        A.data$auxiliary$passed.var$data.file <- paste(data.file, ".Rdata", sep = "")
      }
    
      A.data$auxiliary$shiny.input <- eval(input)
      old <- A.data
      save(old, file = paste(data.file, ".Rdata", sep = ""))
    }else{
      return()
    }
    
    
  })
  

  
  observe({
    if(input$changesavefileButton == 0)
      return()
    
    isolate({
      data.file <- input$datafile
      data.path <- 
        paste(strsplit(data.file, "/")[[1]][1:(length(strsplit(data.file, "/")[[1]])-1)], 
              collapse = "/")
      old.wd <- getwd()
      setwd(data.path)   
      file <- tclvalue(tkgetSaveFile(initialfile="Absorbance"))
      if(file != "")
        data.file <- file
      
      if(DoseR:::fileExt(data.file) %in% c("xls", "xlsx", "RData", "Rdata"))
        data.file <-  gsub(paste(".", DoseR:::fileExt(data.file), sep = ""), "", data.file)
      x <- input$controller
      
      updateTextInput(session, "datafile", value = paste(data.file, x, sep = ""))  
      setwd(old.wd)
    })
  })
  
  
  observe({
    if(input$oldfilebutton == 0){
      analysisfile <- "Choose file"
    }else{
      analysisfile <- tk_choose.files()#try(file.choose(), silent=TRUE)
        
      #if(is(analysisfile, 'try-error')){
      if(length(analysisfile) == 0 ){
        "choose file"
      }else{
        print(analysisfile)
        
        load(analysisfile)  
        
        x <- input$controller
        
        updateTextInput(session, "oldfile",  value = paste(analysisfile, x, sep = ""))
        
        A.data <<- old
        
        x <- input$controller
        y <- old$auxiliary$shiny.input
        
        updateTextInput(session, "datafile", value = paste(y$datafile, x, sep = ""))
        
        
        ############################################
        ##
        ##  update the input variables for createmetada function 
        ##
        ###########################################
        
        updateCheckboxInput(session, "CMDcreateMetaData",
                            value = y$CMDcreateMetaData)
        
        updateTextInput(session, "CMDdbfdir", 
                        value = paste(y$CMDdbfdir, x, sep = ""))
        
        updateTextInput(session, "CMDprotocoldir", 
                        value = paste(y$CMDprotocoldir, x, sep = ""))                    
        
        updateTextInput(session, "CMDcorrectionname", 
                        value = paste(y$CMDcorrectionname, x, sep = ""))
        
        
        l<-  list("Choose" = "non",
                  "Directly from file names" = "filenames", 
                  "From an existing Excel sheet" = "excel")
        
        selected <- names(l)[(l == y$CMDmetadatafromfile)]
        updateSelectInput(session,
                          "CMDmetadatafromfile", 
                          choices= l,                      
                          selected = selected)
        
        updateTextInput(session, "CMDmetadataexcel", 
                        value = paste(y$CMDmetadataexcel, x, sep = ""))
        
        
        updateTextInput(session, "CMDfileColnames",  
                        value = paste(y$CMDfileColnames, x, sep = ""))
        
        l <-   c("Not separated",
                 "Comma (,)"=',',
                 "Semicolon (;)"=';',
                 "Underscore (_)"='_')
        
        selected <- names(l)[(l == y$CMDsep)]
        
        updateRadioButtons(session, 'CMDsep', 
                           choices= l,                      
                           selected = selected)
        #         
        #         
        l <-  list("Long" = "long",
                   "Wide" = "wide")
        
        selected <- names(l)[(l == y$CMDmetadataformat)]
        updateRadioButtons(session, "CMDmetadataformat", 
                           choices= l,                      
                           selected = selected)
        
        updateCheckboxInput(session, "CMDeditcolnames", 
                            value = y$CMDeditcolnames)
        
        updateTextInput(session, "CMDnamevar",
                        value = paste(y$CMDnamevar, x, sep = ""))
        #       
        updateTextInput(session, "CMDdrugvar", 
                        value = paste(y$CMDdrugvar, x, sep = ""))
        
        updateTextInput(session, "CMDprotocolvar",
                        value = paste(y$CMDprotocolvar, x, sep = ""))
        #       
        updateTextInput(session,"CMDtimevar",
                        value = paste(y$CMDtimevar, x, sep = ""))
        #       
        updateTextInput(session, "CMDidentifier",
                        value = paste(y$CMDidentifier, x, sep = ""))
        
        updateCheckboxInput(session, "CMDspecialidentifier", 
                            value = y$CMDspecialidentifier)
        
        l <- c("Comma (,)"=',',
               "Semicolon (;)"=';',
               "Underscore (_)"='_')
        
        selected <- names(l)[(l == y$CMDidentifiersep)]
        updateRadioButtons(session, "CMDidentifiersep", 
                           choices= l,                      
                           selected = selected)
        
        l <-  list("Not a date"    = "empty",
                   "ddmmyy (150812)" = "%d%m%y",
                   "m/d/y (02/27/92)"= "%m/%d/%y",
                   "ddmmmyyyy (2jan1960)" = "%d%b%Y")
        
        selected <- names(l)[l == y$CMDdateformat]
        
        updateSelectInput(session, "CMDdateformat", 
                          choices= l,                      
                          selected = selected)
        
        
        updateTextInput(session, "CMDincubationvar", 
                        value = paste(y$CMDincubationvar, x, sep = ""))                    
        
        updateTextInput(session, "CMDdoublingvar", 
                        value = paste(y$CMDdoublingvar, x, sep = ""))
        
        l <- list("micro g/ml"    = "ug/ml",
                  "g/l"      = "g/l",
                  "mol/l"    = "mol/l")
        
        selected <- names(l)[(l == y$CMDcreateMetaDataunit)]
        
        updateSelectInput(session, "CMDcreateMetaDataunit", 
                          choices= l,                      
                          selected = selected)
        
        updateCheckboxInput(session, "CMDextrametadatafromfile",
                            value = y$CMDextrametadatafromfile)
        
        updateTextInput(session, "CMDadditionalmetadata",
                        value = paste(y$CMDadditionalmetadata, x, sep = ""))
        
        
        
        ############################################
        ##
        ##  Update the readDBFdata call
        ##
        ###########################################
        
        updateCheckboxInput(session, "DBFreadDBFData",
                            value = y$DBFreadDBFData)
        
        updateTextInput(session, "DBFdiscardlines",
                        value = paste(y$DBFdiscardlines, x, sep = ""))
        
        updateTextInput(session, "DBFabsorbanceid",
                        value = paste(y$DBFabsorbanceid, x, sep = ""))
        
        updateTextInput(session, "DBFwellid",
                        value = paste(y$DBFwellid, x, sep = ""))
        
        updateTextInput(session, "DBFremoverows",
                        value = paste(y$DBFremoverows, x, sep = ""))
        
        updateTextInput(session, "DBFremovecols",
                        value = paste(y$DBFremovecols, x, sep = ""))
        
        updateTextInput(session, "DBFdosevar",
                        value = paste(y$DBFdosevar, x, sep = ""))
        
        updateTextInput(session, "DBFadditivevar",
                        value = paste(y$DBFadditivevar, x, sep = ""))
        
        updateTextInput(session, "DBFcontrolval",
                        value = paste(y$DBFcontrolval, x, sep = ""))
        
        updateTextInput(session, "DBFbackgroundval",
                        value = paste(y$DBFbackgroundval, x, sep = ""))
        
        updateTextInput(session, "DBFmistakeval",
                        value = paste(y$DBFmistakeval, x, sep = ""))
        
        
        updateCheckboxInput(session, "DBFupdate",
                            value = y$DBFupdate)
        
        
        ############################################
        ##
        ##  Update the Background colour correction
        ##
        ###########################################
        
      
        updateNumericInput(session, "DCCoutlierTest", 
                           value = y$DCCoutlierTest)
        updateNumericInput(session, "DCCoutlierIter", 
                           value = y$DCCoutlierIter)
        updateNumericInput(session, "DCCvarpowerMin", 
                           value = y$DCCvarpowerMin)
        updateNumericInput(session, "DCCvarpowerIter", 
                           value = y$DCCvarpowerIter)
        
        
        l <-  list("Unrestricted" = "unrestricted",
                   "Restricted"   = "restricted")
        
        selected <- names(l)[l == y$DCCpararametrisation]
        
        updateSelectInput(session, "DCCpararametrisation", 
                          choices= l,                      
                          selected = selected)
        
        
        l <-  list("Sum"       = "sum",
                   "Helmert"   = "helmert", 
                   "Treatment" = "treatment")
        
        selected <- names(l)[l == y$DCCcontr]
        
        updateSelectInput(session, "DCCcontr", 
                          choices= l,                      
                          selected = selected)
        
        updateCheckboxInput(session, "DCCweights",
                            value = y$DCCweights)
        updateCheckboxInput(session, "DCCfittedA",
                            value = y$DCCfittedA)
        updateCheckboxInput(session, "DCCupdate",
                            value = y$DCCupdate)
        
        
        ############################################
        ##
        ##  Update the Normalisation
        ##
        ###########################################
        
        
        updateNumericInput(session, "BGNoutlierTest", 
                           value = y$BGNoutlierTest)
        updateNumericInput(session, "BGNoutlierIter", 
                           value = y$BGNoutlierIter)
        updateNumericInput(session, "BGNvarpowerMin", 
                           value = y$BGNvarpowerMin)
        updateNumericInput(session, "BGNvarpowerIter", 
                           value = y$BGNvarpowerIter)
        
        
        l <-  list("Unrestricted" = "unrestricted",
                   "Restricted"   = "restricted")
        
        selected <- names(l)[l == y$BGNpararametrisation]
        
        updateSelectInput(session, "BGNpararametrisation", 
                          choices= l,                      
                          selected = selected)
        
        
        l <-  list("Sum"       = "sum",
                   "Helmert"   = "helmert", 
                   "Treatment" = "treatment")
        
        selected <- names(l)[l == y$BGNcontr]
        
        updateSelectInput(session, "BGNcontr", 
                          choices= l,                      
                          selected = selected)
        
        updateCheckboxInput(session, "BGNweights",
                            value = y$BGNweights)
        updateCheckboxInput(session, "BGNfittedA",
                            value = y$BGNfittedA)
        updateCheckboxInput(session, "BGNupdate",
                            value = y$BGNupdate)
        
      
      
      
        ############################################
        ##
        ##  Update the Bootstrap
        ##
        ###########################################
        
        
        updateNumericInput(session, "Bootnsamples", 
                           value = y$Bootnsamples)
        updateNumericInput(session, "Bootnsamplesmax", 
                           value = y$Bootnsamplesmax)
        
        l <- list("Parametric"       = "parametric",
                  "Residual"   = "residual")
        
        selected <- names(l)[l == y$Boottype]
        
        updateSelectInput(session, "Boottype", 
                          choices= l,                      
                          selected = selected)  
                  
#                   
#                   
#         ############################################
#         ##
#         ##  Update the Bootstrap
#         ##
#         ###########################################
#         
#         
#     
#         l <- list(c("G", "R", "D", "DG", "RG"))
#         
#         selected <- names(l)[l == y$DRmodels]
#         
# #         updateCheckboxGroupInput(session, "DRmodels", 
# #                           choices= l,                      
# #                           selected = selected)  
# #         
#         
#         updateCheckboxGroupInput(session, "DRmodels", label = NULL,
#                                  choices = NULL, selected = NULL)
#         
        
        l <-  list("Restricted"       = "restricted",
                   "Non restricted"   = "nonrestricted")
        
        selected <- names(l)[l == y$DRparametrisation]
        
        updateSelectInput(session, "DRparametrisation", 
                          choices= l,                      
                          selected = selected)  
        
        updateNumericInput(session, "DRcut", 
                           value = y$DRcut)
        
        updateTextInput(session, "DRdoublingvar",
                        value = paste(y$DRdoublingvar, x, sep = ""))

        updateTextInput(session, "DRunit",
                        value = paste(y$DRunit, x, sep = ""))
      
        
        
        l <-  list("log base 10" = "log10",
                   "log base 2"  = "log2",
                   "log base e"  = "log",
                   "no log"      = "nolog")
        
        selected <- names(l)[l == y$DRlogfrom]
        
        updateSelectInput(session, "DRlogfrom", 
                          choices= l,                      
                          selected = selected) 
      
        
        l <-  list("log base 10" = "log10",
                   "log base 2"  = "log2",
                   "log base e"  = "log",
                   "no log"      = "nolog")
        
        selected <- names(l)[l == y$DRlogto]
        
        updateSelectInput(session, "DRlogto", 
                          choices= l,                      
                          selected = selected)  
        
        
        
        
        updateNumericInput(session, "DRt", 
                           value = y$DRt)
        
        updateNumericInput(session, "DRAUCq", 
                           value = y$DRAUCq)           
      }
      
      
      output$CMDprint <- renderPrint({
        input$oldfilebutton
        if("createMetaData" %in% class(A.data)){
           print.createMetaData(A.data)
        }else{
          print("Create the meta data")
        }
      })
        
      output$DBFreadDBFDataPrint <- renderPrint({
        input$oldfilebutton
        if("Absorbance" %in% class(A.data)){
          print.Absorbance(A.data)
        }else{
          print("Load the dBase files")
        }
      })
      
      
      output$DCCPrint <- renderPrint({
        input$oldfilebutton
        if("drugColorCorrection" %in% class(A.data)){
          print.drugColorCorrection(A.data)
        }else{
          print("Run the drug color correction program")
        }
      })
      
      output$BGNPrint <- renderPrint({
        input$oldfilebutton
        if("bgModel" %in% class(A.data)){
          print.bgModel(A.data)
        }else{
          print("Run the drug pre-processing program")
        }
      })
        
      
      output$BootPrint <- renderPrint({
        input$oldfilebutton
        if("bootstrap" %in% class(A.data)){
          print.bootstrap(A.data)
        }else{
          print("Run the drug bootstrap algorithm")
        }
      })
      
      output$DRPrint <- renderPrint({
        input$oldfilebutton
        if("doseResponseModel" %in% class(A.data)){
          print.doseResponseModel(A.data)
        }else{
          print("Run the drug bootstrap algorithm")
        }
      })

      
      #########################################
    }    
  })
  
  ############################################
  ##
  ##  Welcome screen
  ##
  ###########################################
  
  output$AnalysisFlow <- renderPrint({
    input$CMDCreateMetaDataButton
    input$oldfilebutton
    input$DBFreadDBFDataButton
    input$DCCanalysisButton
    input$BGNanalysisButton
    input$BootanalysisButton
    input$DRanalysisButton
    #print("amders")
    DoseR:::printAnalysisFlow(A.data)
    cat("A.data <-", A.data$call$bgModel$shiny.call)
  })
  
  output$Citation <- renderPrint({
    citation("DoseResponse")
  })
  

  
  ############################################
  ##
  ##  Create metadata
  ##
  ###########################################
  observe({
    input$CMDdbfdirbutton == 0
    isolate({
    dbfdir <- input$CMDdbfdir
    if(input$CMDdbfdirbutton == 0){
      dbfdir <- getwd()
    }else{
      dbfdir <-  try(DoseR:::choose.dir(old.dir = dbfdir))
    }
    x <- input$controller
    updateTextInput(session, "CMDdbfdir",  value = paste(dbfdir, x, sep = ""))
    })
  })
  
  observe({
    input$CMDprotocoldirbutton
    isolate({
    protocoldir <- input$CMDprotocoldir
    if(input$CMDprotocoldirbutton == 0){
      protocoldir <-  getwd()
                        
    }else{
      protocoldir <-  try(DoseR:::choose.dir(old.dir = protocoldir))
    }
    x <- input$controller
    updateTextInput(session, "CMDprotocoldir",  value = paste(protocoldir, x, sep = "")) 
    })
  })
  
  
  observe({
    if(input$CMDmetadataexcelbutton == 0){
      if(is.null(A.data$call$createMetaData$data)){
        metadataexcel <- "file"
      }else{
        if(class(A.data$call$createMetaData$data) == "character")
          metadataexcel <- A.data$call$createMetaData$data
      }
    }else{
      metadataexcel <- DoseR:::my.file.choose(old.file = metadataexcel)
    }
    x <- input$controller
    updateTextInput(session, "CMDmetadataexcel",  value = paste(metadataexcel, x, sep = "")) 
    
    if(file.exists(metadataexcel)){
    metadata <- 
      read.xls(metadataexcel,
               fileEncoding = "latin1",
               method = "tab",
               stringsAsFactors = FALSE) 
    
    output$CMDexample <- renderTable({
      metadata
    })}
  })
  
  
  
  observe({
    if(input$CMDadditionalmetadatabutton == 0){
      if(is.null(A.data$call$createMetaData$additional.metadata)){
        additional.metadata <- "file"
      }else{
        if(class(A.data$call$createMetaData$additional.metadata) == "character")
          additional.metadata <- A.data$call$createMetaData$additional.metadata
      }
    }else{
      additional.metadata <-  DoseR:::my.file.choose(old.file = input$CMDadditionalmetadata)
    }
    x <- input$controller
    updateTextInput(session, "CMDadditionalmetadata",  
                    value = paste(additional.metadata, x, sep = ""))  
  })
  
  
     
  
  observe({
    if(input$CMDCreateMetaDataButton == 0)
      return()
    
   # updateTabsetPanel(session, "panelCMDcreateMetaData", selected = input$Print)
    
    A.data <<- isolate({ 
      
      withProgress(session, min=1, max=15, {
        setProgress(message = 'Calculation in progress',
                    detail = 'This may take a while...')
    
        metadataexcel <- input$CMDmetadataexcel
        
    if(metadataexcel ==  "file"){
      data <- NULL
    }else{
      data <- metadataexcel
    }
        additional.metadata<-input$CMDadditionalmetadata
    if(additional.metadata ==  "file"){
      additional.metadata <- NULL
    }else{
      additional.metadata <- additional.metadata
    }
    
    if(input$CMDdateformat ==  "empty"){
      dateformat <- NULL
    }else{
      dateformat <- input$CMDdateformat
    }
    
    if(input$CMDdoublingvar == "not available"){
      doblingvar <- NULL
    }else{
      doblingvar <- input$CMDdoublingvar
    }
  
   
    do.call(createMetaData, list(
      data             = quote(data),
      data.file        = eval(input$datafile),
      namevar          = eval(input$CMDnamevar),
      drugvar          = eval(input$CMDdrugvar),
      protocolvar      = eval(input$CMDprotocolvar),
      identifier       = eval(input$CMDidentifier),
      timevar          = eval(input$CMDtimevar),
      correctionname   = eval(input$CMDcorrectionname),
      incubationvar    = eval(input$CMDincubationvar),
      doublingvar      = eval(doblingvar),
      format           = eval(input$CMDmetadataformat),
      dbf.path         = eval(input$CMDdbfdir),
      protocol.path    = eval(input$CMDprotocoldir),
      dbf.files        = NULL,
      file.extension = ".dbf",
      protocol.files   = NULL,
      colnames         = eval(strsplit(input$CMDfileColnames, ";")[[1]]),
      sep              = eval(input$CMDsep),
      namesep          = " ",
      identifiersep    = eval(input$CMDidentifiersep),
      date.format      = eval(dateformat),
      unit             = eval(input$CMDcreateMetaDataunit),
      additional.metadata = eval(additional.metadata),
      show.warnings    = TRUE,
      update           = eval(!input$CMDupdate),
      idcomb           = NULL,#c(namevar, drugvar, identifier) # vars used to generate id variable
      idvar            = "sampleid", 
      namecols         = NULL,
      dbf.file.var     = "dbf.file.name",
      shiny.input      = eval(input)))
    })
    })
   
  })
  
  output$CMDprint <- renderPrint({
    input$CMDCreateMetaDataButton
    A.data
  })
    
  output$CMDexample <- renderTable({
    
    additional.metadata<-input$CMDadditionalmetadata
    
    if(additional.metadata ==  "file"){
      additional.metadata <- NULL
    }else{
      additional.metadata <- additional.metadata
    }
    
    if(input$CMDdateformat ==  "empty"){
      dateformat <- NULL
    }else{
      dateformat <- input$CMDdateformat
    }
    
    if(input$CMDdoublingvar == "not available"){
      doblingvar <- NULL
    }else{
      doblingvar <- input$CMDdoublingvar
    }
    
    
    colnames <- eval(strsplit(input$CMDfileColnames, ";")[[1]])
    namevar  <- input$CMDnamevar
    drugvar <- input$CMDdrugvar
    identifier <- input$CMDidentifier
    timevar <- input$CMDtimevar
    protocolvar <- input$CMDprotocolvar
    sep <- input$CMDsep
    incubationvar    = eval(input$CMDincubationvar)
    dbf.file.var     = "dbf.file.name"
    unit             = input$CMDcreateMetaDataunit
    idvar            = "sampleid"
    format           = eval(input$CMDmetadataformat)
    namesep          = " "
    namevar2 <- paste(namevar, 2, sep = "")
    namecols         = NULL
    
    
    
    if(class(additional.metadata) == "character"){
      tr <- "error"
      if(file.exists(additional.metadata))  
        try(tr <- read.xls(additional.metadata ,
                           fileEncoding = "latin1",
                           method = "tab",
                           stringsAsFactors = FALSE))
      
      if(tr[1] != "error"){
        additional.metadata <- tr
      }else{
        additional.metadata = NULL
      }
    }
    
    dbf.files <- dir(input$CMDdbfdir, pattern = ".dbf", full.names=TRUE)
    ll <- min(20, length(dbf.files))
    if(length(dbf.files) > 0 & input$CMDmetadatafromfile == "filenames"){
      
      updateTabsetPanel(session, "panelCMDcreateMetaData", selected = "Example")
      
      
      dbf.files <-  dbf.files[1:ll]
      dbf.files4 <- basename(dbf.files)
      dbf.path2  <- dirname(dbf.files)
      
      dbf.files2 <- strsplit(dbf.files4,"\\.")
      
      ext2 <- ".dbf"
      for(i in ext2)
        dbf.files3 <- gsub(i, "", dbf.files4)
      
      spl <- strsplit(dbf.files3, sep)      
      
      
      met <- t(sapply(spl, t))
      
      if(ncol(met) == length(colnames))
      colnames(met) <- colnames
      if(nrow(met) == 1)
        met <- t(met)
      met <- as.data.frame(met)
      met
      colnamefit <- ncol(met) == length(colnames)
      if(colnamefit)
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
       
#        if(!is.null(date.format)){
#          if(!is.null(identifiersep)){
#            met[, "setupdate"] <- 
#              as.character(as.Date(unlist(lapply(strsplit(as.character(met[, identifier]), identifiersep), 
#                                    function(x) x[1])), format = date.format))
#            met[, "replicate"]<- 
#              unlist(lapply(strsplit(as.character(met[, identifier]), identifiersep), function(x) x[2]))
#          }else{
#            met[, "setupdate"] <- as.character(as.Date(as.character(met[, identifier]), format = date.format))
#          }
#        }
#       
#       
      #met[, "id"] <- paste(met[, namevar], met[, drugvar], met[, identifier], sep = " ")
      #if(identifier != "id")
      #  met <- met[, colnames(met) %w/o% identifier]
      
      ## the identifier is changed to id as this is used in following functions
      #identifier <- "id"
     # if(colnamefit){
      met[, dbf.file.var]  <- dbf.files3
      met[, "date"] <- as.character((file.info(dbf.files)$mtime))
      
      met[, "unit"] <- unit
      met[, "omit"] <- FALSE
      met[, "reason"] <- ""
      
      
      
      
      if(! incubationvar %in% colnames(met))
        met[, incubationvar] <- 2
      
      met.wide <- met
      
      if(colnamefit){
      if(format == "wide"){
        idcomb  <- c(namevar, drugvar, identifier)
      }else{
        idcomb  <- c(namevar, drugvar, identifier, timevar)
      }
      
      
      idcombwide  <- c(namevar, drugvar, identifier)
      
      idcomblong  <- c(namevar, drugvar, identifier, timevar)
      
      
      met[, idvar] <- apply(met[, idcomb], 1, function(x) paste(x, collapse = " "))
      
      met[, "idwide"] <- apply(met[, idcombwide], 1, function(x) paste(x, collapse = " "))
      met[, "idlong"] <- apply(met[, idcomblong], 1, function(x) paste(x, collapse = " "))
      
      met.wide
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
        rownames(additional.metadata) <- additional.metadata[, namevar]
        cols <-colnames(additional.metadata) %w/o% namevar
        met.wide[, cols] <-  additional.metadata[met.wide[, namevar2], cols]
      }
      if(format == "long"){
        met.wide[, timevar] <- as.numeric(met.wide[,timevar])
      }else{
        #timevars <- colnames(met.wide)[grepl(timevar, colnames(met.wide))]
        #met.wide[, timevars] <- apply(met.wide[,timevars ], 2, as.numeric)
      }
      }
      
      met.wide
      
    }else{
      return()
    }
  })
    

  
  
  output$CMDmetadataTable <- renderTable({
    input$CMDCreateMetaDataButton
    input$oldfilebutton
    
    
    data <-A.data$meta.list$metadata.new
    if(!is.null(data))
      for(i in 1:ncol(data)){
        if(class(data[,i])[1] == "POSIXct" | class(data[,i])[1] == "Date")
          data[,i] <- as.character(data[,i])
      }
    data
  })
  
  output$CMDmetaCorrection <- renderTable({
    input$CMDCreateMetaDataButton
    input$oldfilebutton
    data <-A.data$meta.list$metadata.correction
    if(!is.null(data))
      for(i in 1:ncol(data)){
        if(class(data[,i])[1] == "POSIXct" | class(data[,i])[1] == "Date")
          data[,i] <- as.character(data[,i])
      }
    data
  })
  
  output$CMDmetaadditional <- renderTable({
    input$CMDCreateMetaDataButton
    input$oldfilebutton
    A.data$meta.list$additional
  })
  
  output$CMDcreateMetaDataCall <- renderPrint({
    input$CMDCreateMetaDataButton
    input$oldfilebutton
    cat("A.data <-", A.data$auxiliary$shiny.calls$createMetaData)
  })
  

 
  ################################################
  ##
  ## Protocol Editor
  ##
  ################################################
  
  
  observe({  
    input$PRCusethecreatedbutton
    
    
    protocol.files <- dir(input$CMDprotocoldir, pattern = ".xls", full.names=TRUE)
    if(length(protocol.files)){
      l <- list()
      
      for(i in 1:length(protocol.files))
        l[[paste(basename(protocol.files[i]))]] <- protocol.files[i]
      
      updateSelectInput(session, "PRCchooseProtocol",
                        choices=l,
                        select=basename(protocol.files[1]))
      
      updateSelectInput(session, "PRCcreatechooseProtocol",
                        choices=l,
                        select=basename(protocol.files[1]))
      
      updateSelectInput(session, "DBFchooseprotocol",
                        choices=l,
                        select=basename(protocol.files[1]))           
    }
  })
  
  
  observe({
    if(input$PRCCreateeditProtocolButton == 0)
      return() 
  })
  

  observe({
    # if(class(A.data)[1] != "list")
    if(file.exists(input$PRCcreatechooseProtocol))
     # current.protocol <<- try(readProtocol(input$PRCcreatechooseProtocol),
      #                         silent = TRUE)
    
    current.protocol <<- readProtocol(input$PRCcreatechooseProtocol)
  })
  
  
  observe({ 
    current.protocol <<- 
      createProtocol(n.doses = input$PRCNumberofDoses, 
                     n.backgrounds = input$PRCNumberofBackgrounds,
                     n.controls = input$PRCNumberofControls,
                     fold = input$PRCdoseFold,
                     max = input$PRCMaxDose,
                     remove.edges = !input$PRCCreateRemoveEdges,
                     drug = input$PRCcreateProtocoldrug)
  })
  
  
  
   observe({
     if(input$PRCCreateeditProtocolButton == 0)
       return()
     
     protocolfile <<- editProtocol(current.protocol)
   })
  
  
  observe({
    if(input$PRCCreatereloadProtocolButton == 0)
      return()
    
    current.protocol <<- readProtocol(protocolfile)
    
  })
  

  output$PRCCurrentSetup <- renderTable({
    input$PRCCreatereloadProtocolButton
    input$CMDCreateMetaDataButton
    input$oldfilebutton
    input$oldfilebutton
    input$PRCcreatechooseProtocol
    n.doses = input$PRCNumberofDoses
    n.backgrounds = input$PRCNumberofBackgrounds
    n.controls = input$PRCNumberofControls
    fold = input$PRCdoseFold
    max = input$PRCMaxDose
    remove.edges = input$PRCCreateRemoveEdges
    
    current.protocol[["setup"]]
  })
  

  output$PRCCurrentConc <- renderTable({
    input$PRCCreatereloadProtocolButton
    input$CMDCreateMetaDataButton
    input$oldfilebutton
    input$PRCcreatechooseProtocol
    n.doses = input$PRCNumberofDoses
    n.backgrounds = input$PRCNumberofBackgrounds
    n.controls = input$PRCNumberofControls
    fold = input$PRCdoseFold
    max = input$PRCMaxDose
    remove.edges = input$PRCCreateRemoveEdges
    
    current.protocol[["conc"]]
  })
  
  
  
  observe({
    if(input$PRCusethecreatedbutton == 0)
      return() 
    A.data <<-
      isolate(
        useTempProtocol(A.data, protocol = current.protocol, 
                        input$PRCcreateProtocolName, 
                        DoseR:::no.extension.vec(basename(input$PRCCreatechoosedbf)))
      )
  })
  
  
  observe({  
    if(input$PRCusethecreatedbutton == 0)
      return()
    
    dbf.file <- input$PRCCreatechoosedbf
    
    if("createMetaData" %in% class(A.data)){
      file.info <- A.data$auxiliary$file.info
     # file.info$dbf
      
      protocolvar   <- A.data$call$createMetaData$protocolvar
      identifiervar <- A.data$call$createMetaData$identifier
      dbf.file.var  <-  A.data$call$createMetaData$dbf.file.var
      
      file.info <- A.data$auxiliary$file.info
      
      dbf.file.name <- rownames(file.info$dbf[
        file.info$dbf$file.full == dbf.file,,drop = FALSE])
      
      protocol <- A.data$meta.list$metadata.full[
        A.data$meta.list$metadata.full[, dbf.file.var] == dbf.file.name, protocolvar]
      
      ident <- A.data$meta.list$metadata.full[
        A.data$meta.list$metadata.full[, dbf.file.var] == dbf.file.name, identifiervar]
      
      
      protocol.files <- dir(input$CMDprotocoldir, pattern = ".xls", full.names=TRUE)
      if(length(protocol.files)){
        l <- list()
        
        for(i in 1:length(protocol.files))
          l[[paste(basename(protocol.files[i]))]] <- protocol.files[i]
      }
      
      updateSelectInput(session, "PRCcreatechooseProtocol",
                        choices = l,
                        select=paste(protocol,".xls", sep = ""))
      
      x <- input$controller
      new.name <- paste(protocol, ident, sep = "_")
      updateTextInput(session, "PRCcreateProtocolName",  value = paste(new.name, x, sep = ""))

      
    }  
  })
  
  ##########################################
## DBF file reader
##

  
  
  
  observe({
    if(input$DBFreadDBFDataButton == 0)
      return()
    
    
    A.data <<- isolate({ 
      
        do.call(readDBFDataShiny, 
                list(
                  A.data         = quote(A.data),  # Absorbance, data.frame, list  
                  update         = eval(!input$DBFupdate), 
                  discard.lines  = eval(as.numeric(strsplit(input$DBFdiscardlines, ";")[[1]])),   
                  well.id        = eval(input$DBFwellid), 
                  absorbance.id  = eval(input$DBFabsorbanceid),    
                  remove.rows    = eval(strsplit(input$DBFremoverows, ";")[[1]]),
                  remove.cols    = eval(as.numeric(strsplit(input$DBFremovecols,  ";")[[1]])), 
                  verbose        = FALSE,     
                  dosevar        = eval(input$DBFdosevar),
                  additivevar    = eval(input$DBFadditivevar),
                  controlval     = eval(input$DBFcontrolval),
                  backgroundval  = eval(input$DBFbackgroundval), 
                  mistakeval     = eval(input$DBFmistakeval),
                  progressbar    = "window",
                  shiny.input    = input,
                  session = session))
        
      })
    output$DBFreadDBFDataPrint <- renderPrint({
      A.data
      
    })
  })
  
  
  
  
  output$dataRaw <- renderDataTable({ 
    input$oldfilebutton
    input$DBFreadDBFDataButton
    
    if("Absorbance" %in% class(A.data)){
     return(A.data$data$raw.data)
    }else{
      return(NULL)
    }
  })
  
  observe({  
    input$PRCusethecreatedbutton
    dbf.files <- dir(input$CMDdbfdir, pattern = ".dbf", full.names=TRUE)
    if(length(dbf.files)){
      l <- list()
      for(i in 1:length(dbf.files))
        l[[paste(basename(dbf.files[i]))]] <- dbf.files[i]
      
      updateSelectInput(session, "DBFchoosedbf",
                        choices=l,
                        select=basename(dbf.files[1]))
      
      updateSelectInput(session, "PRCchoosedbf",
                        choices=l,
                        select=basename(dbf.files[1]))
      
      updateSelectInput(session, "PRCCreatechoosedbf",
                        choices=l,
                        select=basename(dbf.files[1]))
      
      updateSelectInput(session, "PRCEditchoosedbf",
                        choices=l,
                        select=basename(dbf.files[1]))
      
      
    }
  })
  
    
  output$PRCdbfeksample <- renderTable({
    input$CMDCreateMetaDataButton
    input$oldfilebutton
    file <- input$PRCEditchoosedbf
    if(class(file) == "character")
    if(file.exists(c(file)[1])){
      cell0 <- read.dbf(file, as.is=TRUE)
      rem <- as.numeric(eval(strsplit(input$DBFdiscardlines, ";")[[1]]))
      if(length(rem) > 0 )
        cell0 <- cell0[ -rem, ]
      cell0
      
    }
  })
  
  output$DBFreadDBFDataexsfilep <- renderPrint({
    input$CMDCreateMetaDataButton
    input$oldfilebutton
    file <- input$DBFchoosedbf
    file
  })
  
  
  output$DBFeksample <- renderTable({
    input$CMDCreateMetaDataButton
    input$oldfilebutton
    file <- input$DBFchoosedbf
    if(class(file) == "character")
    if(file.exists(file)){
      cell0 <- read.dbf(file, as.is=TRUE)
      rem <- as.numeric(eval(strsplit(input$DBFdiscardlines, ";")[[1]]))
      if(length(rem) > 0 )
        cell0 <- cell0[ -rem, ]
      cell0
      
    }
  })
  

  
  observe({
    # if(class(A.data)[1] != "list")
    if(file.exists(input$DBFchooseprotocol))
      DBF.current.protocol <<- try(readProtocol(input$DBFchooseprotocol),
                               silent = TRUE)
    
  })
  

  output$DBFprotocolSetup <- renderTable({
    input$CMDCreateMetaDataButton
    input$oldfilebutton
    input$DBFeditprotocolvar
    input$DBFchooseprotocol
    #try(readProtocol(input$DBFchooseprotocol),
    #    silent = TRUE)
    DBF.current.protocol[["setup"]]
  })
  
  output$DBFprotocolConc <- renderTable({
    input$CMDCreateMetaDataButton
    input$oldfilebutton
    input$DBFeditprotocolvar
    input$DBFchooseprotocol
    DBF.current.protocol[["conc"]]
  })
  
  
  
  output$DBFreadDBFDataCall <- renderPrint({
    input$CMDCreateMetaDataButton
    input$oldfilebutton
    input$DBFreadDBFDataButton
    #print("amders")
    cat("A.data <-", A.data$call$readDBFData$shiny.call)
  })
  
  
  observe({  
    input$oldfilebutton
    input$CMDCreateMetaDataButton
    
    if("Absorbance" %in% class(A.data)){
      drugvar <- A.data$auxiliary$passed.var$drugvar
      
      drugs <- sort(unique(A.data$data$raw.data[, drugvar]))
      l <- list()
      for(i in 1:length(drugs))
        l[[drugs[i]]] <- drugs[i]
      
      updateSelectInput(session, "Molareditdrug",
                        choices = l,
                        select  = drugs[1])
    }
    
  })
  
  
  output$molarmassTable <- renderTable({
    input$oldfilebutton
    input$DBFreadDBFDataButton
    
    A.data$auxiliary$mol.data
  })
  
  
  output$molarmassTable <- renderTable({
    if(input$Molareditdrugbutton == 0)
      return()
    if("Absorbance" %in% class(A.data)){
      A.data <<- isolate({ 
        
        do.call(changeMolarMass, 
                list(
                  A.data         = quote(A.data),  # Absorbance, data.frame, list  
                  drug           = input$Molareditdrug,
                  mol.mass       = input$Molareditmass,
                  shiny.input    = input,
                  session = session))
      })
      A.data$auxiliary$mol.data
    }
  })
  
  ##########################################
  ## Drug Color correction
  ##
#   output$DCCPrint <- renderPrint({
#     if(input$DCCanalysisButton == 0)
#       return() 
#     
#     withProgress(session, min = 0,
#                  max = 15, {
#                    setProgress(message = 'Calculation in progress',
#                                detail = 'This may take a while...')   
#     A.data <<- isolate({ 
#       if(input$DCCweights){
#         weights <- "fitted"
#       }else{
#         weights <- NULL
#       }
#         
#       do.call(drugColorCorrection, 
#               list(
#                 A.data          = quote(A.data),  # Absorbance, data.frame, list  
#                 update          = eval(!input$DCCupdate), 
#                 weights         = weights, 
#                 fitted.a        = input$DCCfittedA, # used to be FALSE
#                 contr           = input$DCCcontr,
#                 outlier.test    = input$DCCoutlierTest,
#                 outlier.iter    = input$DCCoutlierIter,
#                 varpower.min    = input$DCCvarpowerMin,
#                 varpower.iter   = input$DCCvarpowerIter,
#                 parametrisation = input$DCCpararametrisation               
#                 )
#     )
#     })
#                  })
#     A.data
#   })
  
  
  observe({
    if(input$DCCanalysisButton == 0)
      return() 
    
    A.data <<- isolate({ 
      if(input$DCCweights){
        weights <- "fitted"
      }else{
        weights <- NULL
      }
      
      withProgress(session, min=1, max=15, {
        setProgress(message = 'Calculation in progress',
                    detail = 'This may take a while...')
      
      do.call(drugColorCorrection, 
              list(
                A.data            = quote(A.data),  # Absorbance, data.frame, list  
                update          =  eval(!input$DCCupdate), 
                weights         = weights, 
                fitted.a        = input$DCCfittedA, # used to be FALSE
                contr           = input$DCCcontr,
                outlier.test    = input$DCCoutlierTest,
                outlier.iter    = input$DCCoutlierIter,
                varpower.min    = input$DCCvarpowerMin,
                varpower.iter   = input$DCCvarpowerIter,
                parametrisation = input$DCCpararametrisation,
                progressbar     = "window",
                shiny.input     = input,
                session = session             
              )
      ) 
      })
    })
    output$DCCPrint <- renderPrint({
    A.data
    })
  })
  
  output$DCCcall <- renderPrint({
    input$CMDCreateMetaDataButton
    input$oldfilebutton
    input$DBFreadDBFDataButton
    input$DCCanalysisButton
    
    #print("amders")
    cat("A.data <-", A.data$call$drugColorCorrection$shiny.call)
  })
  
  observe({  
    input$oldfilebutton
    input$DCCanalysisButton 
    
    if("drugColorCorrection" %in% class(A.data)){

      drugs <- names(A.data$drug.color.correct$background.list)
      l <- list()
      for(i in 1:length(drugs))
        l[[drugs[i]]] <- drugs[i]

      updateSelectInput(session, "DCCdrugselector1",
                        choices = l,
                        select  = drugs[1])
      updateSelectInput(session, "DCCdrugselector2",
                        choices = l,
                        select  = drugs[1])
    }
    
  })

  
  
  output$DCCcolorCorrection <- renderPlot({
    input$DCCdrugselector1
    input$oldfilebutton
    input$DCCanalysisButton
    #
    
    if("drugColorCorrection" %in% class(A.data)){
#       plot(density(rnorm(100)))
      plotdrugColorCorrection(
        set.par=TRUE,
        A.data = A.data, 
        drugs  = input$DCCdrugselector1,
        #file   = file.path(figure.output, "drugColorCorrection.pdf"),
        type.plot = c("output"), 
        pdfit = FALSE)
    }else{
      return()
    }
    
  })
  
  output$DCCcolorCorrectionRaw <- renderPlot({
    input$DCCdrugselector2
    input$oldfilebutton
    input$DCCanalysisButton
    #
    
    if("drugColorCorrection" %in% class(A.data)){
      #       plot(density(rnorm(100)))
      plotdrugColorCorrection(
        A.data = A.data, 
        set.par=TRUE,
        drugs  = input$DCCdrugselector2,
        #file   = file.path(figure.output, "drugColorCorrection.pdf"),
        type.plot = c("raw"), 
        pdfit = FALSE)
    }else{
      return()
    }
    
  })
  
  ##########################################
  ## Normalisation
  ##
  observe({
    if(input$BGNanalysisButton == 0)
      return() 
    
    A.data <<- isolate({ 
      if(input$BGNweights){
        weights <- "fitted"
      }else{
        weights <- NULL
      }
      withProgress(session, min=1, max=15, {
        setProgress(message = 'Calculation in progress',
                    detail = 'This may take a while...')
      do.call(bgModel, 
              list(
                A.data            = quote(A.data),  # Absorbance, data.frame, list  
                update          =  eval(!input$BGNupdate), 
                weights         = "fitted", 
                fitted.a        = input$BGNfittedA, # used to be FALSE
                contr           = input$BGNcontr,
                outlier.test    = input$BGNoutlierTest,
                outlier.iter    = input$BGNoutlierIter,
                varpower.min    = input$BGNvarpowerMin,
                varpower.iter   = input$BGNvarpowerIter,
                parametrisation = input$BGNpararametrisation,
                shiny.input     = input,
                progressbar     = "window",
                session         = session
                
              )
      )
      })
    })
    output$BGNPrint <- renderPrint({
    A.data
    })
  })
  
  

  output$BGNcall <- renderPrint({
    input$CMDCreateMetaDataButton
    input$oldfilebutton
    input$DBFreadDBFDataButton
    input$DCCanalysisButton
    input$BGNanalysisButton
    #print("amders")
    cat("A.data <-", A.data$call$bgModel$shiny.call)
  })
  
  
  observe({  
    input$oldfilebutton
    input$BGNanalysisButton 
    
    if("bgModel" %in% class(A.data)){
      drugvar <- A.data$auxiliary$passed.var$drugvar
      namevar <- A.data$auxiliary$passed.var$namevar
      
      
      drugs <- sort(unique(A.data$data$bc.mean[, drugvar]))
      l <- list()
      for(i in 1:length(drugs))
        l[[drugs[i]]] <- drugs[i]
      
      updateSelectInput(session, "BGNdrugselector1",
                        choices = l,
                        select  = drugs[1])
      updateSelectInput(session, "BGNdrugselector2",
                        choices = l,
                        select  = drugs[1])
      
      data.name <- A.data$data$bc.mean[A.data$data$bc.mean[, drugvar] == drugs[1],] 
      names <- sort(unique(data.name[, namevar]))
      
      l <- list()
      for(i in 1:length(names))
        l[[names[i]]] <- names[i]
      
      updateSelectInput(session, "BGNnameselector1",
                        choices = l,
                        select  = names[1])
      updateSelectInput(session, "BGNnameselector2",
                        choices = l,
                        select  = names[1])
      
    }
    
  })
  
  BGNPointsize1 <- 8
  observe({ 
    input$BGNPointsize1  
    BGNPointsize1 <<- input$BGNPointsize1   
  })
  
  
  observe({  
    input$BGNdrugselector1
    input$BGNdrugselector2
    input$BGNPointsize1  
    
    if("bgModel" %in% class(A.data)){
      drugvar <- A.data$auxiliary$passed.var$drugvar
      namevar <- A.data$auxiliary$passed.var$namevar
      

      drug <- input$BGNdrugselector1
      
      data.name <- A.data$data$bc.mean[A.data$data$bc.mean[, drugvar] == drug,] 
      names <- sort(unique(data.name[, namevar])) 
      
      l <- list()
      for(i in 1:length(names))
        l[[names[i]]] <- names[i]
        
      updateSelectInput(session, "BGNnameselector1",
                        choices = l,
                        select  = names[1])
      
      
      
      drug <- input$BGNdrugselector2
      
      data.name <- A.data$data$bc.mean[A.data$data$bc.mean[, drugvar] == drug,] 
      names <- sort(unique(data.name[, namevar])) 
      
      l <- list()
      for(i in 1:length(names))
        l[[names[i]]] <- names[i]
      

      updateSelectInput(session, "BGNnameselector2",
                        choices = l,
                        select  = names[1])
      
    }
    
  })
  
  BGNmyres1 <- function(){
    input$BGNPPI1
  }
  
  BGNmypointsize1 <- function(){
    input$BGNPointsize1
  }
  
  BGNmyheight1 <- function(){
    input$BGNheight1*300/2.54
  }
  
  BGNmywidth1 <- function(){
    input$BGNwidth1*300/2.54
  }
  
  output$BGNnormalisedData <- renderPlot({
    input$BGNdrugselector1
    input$oldfilebutton
    input$BGNanalysisButton
    
    if("bgModel" %in% class(A.data)){

      if(input$BGNtitle1 == "Standard"){
        title1 = NULL
      }else{
        title1 = input$BGNtitle1
      }
      
      if(input$BGNmain1 == "Standard"){
        main1 = NULL
      }else{
        main1 = input$BGNmain1
      }
      
      if(input$BGNxlab1 == "Standard"){
        xlab = NULL
      }else{
        xlab = input$BGNxlab1
      }
      
      plotbgModel(A.data   = A.data,
                  names    = input$BGNnameselector1,
                  drugs    = input$BGNdrugselector1,
                  main     = main1,
                  title.1  = title1,
                  title.2  = eval(strsplit(input$BGNtitle2, ";")[[1]]),
                  times    = NULL,
                  ylab     = input$BGNylab1,
                  xlab     = xlab,
                  cex      = input$BGNPointsize1/8,
                  cex.axis = input$BGNPointsize1/8,
                  cex.lab  = input$BGNPointsize1/8,
                  cex.main = input$BGNPointsize1/8,
                  cex.sub  = input$BGNPointsize1/8,
                  col      = colscheme$line$col, 
                  unit     = input$BGNunit1,
                  logfun   = input$BGNlogselector1,
                  #pointsize = 8,
                  pdfit = FALSE)

          
    }else{
      return()
    }
    
  }, height = BGNmyheight1, width = BGNmywidth1, 
                                         pointsize = 8, 
                                         res = 300)
  
  
  
  observe({
    input$BGNpdfit1
    
    isolate({
    if("bgModel" %in% class(A.data)){
      
      if(input$BGNtitle1 == "Standard"){
        title1 = NULL
      }else{
        title1 = input$BGNtitle1
      }
      
      if(input$BGNmain1 == "Standard"){
        main1 = NULL
      }else{
        main1 = input$BGNmain1
      }
      
      if(input$BGNxlab1 == "Standard"){
        xlab = NULL
      }else{
        xlab = input$BGNxlab1
      }
      tempfile <- tempfile( fileext = ".pdf")
      pdf(tempfile, 
          height = input$BGNheight1/2.54, 
          width = input$BGNwidth1/2.54,
          pointsize = input$BGNPointsize1)
      plotbgModel(A.data   = A.data,
                   names    = input$BGNnameselector1,
                   drugs    = input$BGNdrugselector1,
                   main     = main1,
                   title.1  = title1,
                   title.2  = eval(strsplit(input$BGNtitle2, ";")[[1]]),
                   times    = NULL,
                   ylab     = input$BGNylab1,
                   xlab     = xlab,
                   col      = colscheme$line$col, 
                   unit     = input$BGNunit1,
                   logfun   = input$BGNlogselector1,
                   #pointsize = 8,
                   pdfit = FALSE)
      dev.off()
      system(paste("open", tempfile))
    }else{
      return()
    }
    
    })
  })
  
  
  observe({
    input$BGNpngit1
    
    isolate({
      if("bgModel" %in% class(A.data)){
        
        if(input$BGNtitle1 == "Standard"){
          title1 = NULL
        }else{
          title1 = input$BGNtitle1
        }
        
        if(input$BGNmain1 == "Standard"){
          main1 = NULL
        }else{
          main1 = input$BGNmain1
        }
        
        if(input$BGNxlab1 == "Standard"){
          xlab = NULL
        }else{
          xlab = input$BGNxlab1
        }
        tempfile <- tempfile( fileext = ".png")
        png(tempfile, 
            height = input$BGNheight1, 
            width = input$BGNwidth1, units = "cm",
            pointsize = input$BGNPointsize1,
            res = input$BGNPPI1)
       
        plotbgModel(A.data    = A.data,
                    names    = input$BGNnameselector1,
                    drugs    = input$BGNdrugselector1,
                    main     = main1,
                    title.1  = title1,
                    title.2  = eval(strsplit(input$BGNtitle2, ";")[[1]]),
                    times    = NULL,
                    ylab     = input$BGNylab1,
                    xlab     = xlab,
                    col      = colscheme$line$col, 
                    unit     = input$BGNunit1,
                    logfun   = input$BGNlogselector1,
                    #pointsize = 8,
                    pdfit = FALSE)
        dev.off()
        system(paste("open", tempfile))
      }else{
        return()
      }
      
    })
  })
  
  
  output$BGNModelcheck <- renderPlot({
    input$BGNdrugselector2
    input$oldfilebutton
    input$BGNanalysisButton
    #
    
    if("bgModel" %in% class(A.data)){
      #       plot(density(rnorm(100)))
      plotbgModelresid(
        A.data = A.data,
        names  = input$BGNnameselector2,
        drugs  = input$BGNdrugselector2,
        col      = colscheme$line$col,
        #file   = file.path(figure.output, "drugColorCorrection.pdf"), 
        pdfit = FALSE)
    }else{
      return()
    }
  })
  
  
  ##########################################
  ## Bootstrap
  ##
  observe({
    if(input$BootanalysisButton == 0)
      return() 
    
    A.data <<- isolate({   
      do.call(bootstrap, 
              list(
                A.data = quote(A.data),   
                update =  eval(!input$Bootupdate), 
                n.samples = input$Bootnsamples, 
                max.iter = max(input$Bootnsamples + 10, input$Bootnsamplesmax),
                type = input$Boottype,
                shiny.input     = input,
                progressbar     = "window",
                session         = session           
              )
      ) 
    })
    output$BootPrint <- renderPrint({
    A.data
    })
  })
  
  output$dataBoot <- renderDataTable({ 
    input$oldfilebutton
    input$BootanalysisButton
    
    if("bootstrap" %in% class(A.data)){
      return(A.data$data$bs.mean)
    }else{
      return(NULL)
    }
  })
 
  
 
  
  
  
  output$bootcall <- renderPrint({
    input$oldfilebutton
    input$BootanalysisButton
    #print("amders")
    cat("A.data <-", A.data$call$bootstrap$shiny.call)
  })
  
  ##########################################
  ## Dose Response Models
  ##
  observe({
    if(input$DRanalysisButton == 0)
      return() 
    
    A.data <<- isolate({   
      do.call(doseResponseModel, 
              list(
                A.data           = quote(A.data),   
                update           = eval(!input$DRupdate), 
                models           = input$DRmodels,  
                dose.scale       = input$DRunit,
                dose.logfun.from = input$DRlogfrom,
                dose.logfun.to   = input$DRlogto,
                parametrisation  = input$DRparametrisation,  
                AUC.q            = input$DRAUCq,
                t                = input$DRt,
                doublingvar      = input$DRdoublingvar,
                cut              = input$DRcut, 
                shiny.input      = input,
                progressbar      = "window",
                session          = session,
                verbose = FALSE           
              )
      ) 
    })
    output$DRPrint <- renderPrint({
    A.data
    })
  })
  
  
  output$DRData <- renderDataTable({ 
    input$oldfilebutton
    input$DRanalysisButton
    
    if("doseResponseModel" %in% class(A.data)){
      return(A.data$data$DR.data)
    }else{
      return(NULL)
    }
  })
  
  output$DRCall <- renderPrint({
    input$oldfilebutton
    input$DRanalysisButton
    #print("amders")
    cat("A.data <-", A.data$call$doseResponseModel$shiny.call)
  })
  
  
  ##########################################
  ## Figures and Tables 
  ##

  

  observe({  
    input$oldfilebutton
    input$DRanalysisButton 
    
    if("doseResponseModel" %in% class(A.data)){
      drugvar <- A.data$auxiliary$passed.var$drugvar
      namevar <- A.data$auxiliary$passed.var$namevar
      
      
      drugs <- sort(unique(A.data$data$bc.mean[, drugvar]))
      l <- list()
      for(i in 1:length(drugs))
        l[[drugs[i]]] <- drugs[i]
      
      updateSelectInput(session, "DRdrugselector1",
                        choices = l,
                        select  = drugs[1])
      
      data.name <- A.data$data$bc.mean[A.data$data$bc.mean[, drugvar] == drugs[1],] 
      names <- sort(unique(data.name[, namevar]))
      
      l <- list()
      for(i in 1:length(names))
        l[[names[i]]] <- names[i]
      
      updateSelectInput(session, "DRnameselector1",
                        choices = l,
                        select  = names[1])
      
    }
    
  })
  
  DRPointsize1 <- 8
  observe({ 
    input$DRPointsize1  
    DRPointsize1 <<- input$DRPointsize1   
  })
  
  
  observe({  
    input$DRdrugselector1
    input$DRPointsize1  
    
    if("doseResponseModel" %in% class(A.data)){
      drugvar <- A.data$auxiliary$passed.var$drugvar
      namevar <- A.data$auxiliary$passed.var$namevar
      
      
      drug <- input$DRdrugselector1
      
      data.name <- A.data$data$bc.mean[A.data$data$bc.mean[, drugvar] == drug,] 
      names <- sort(unique(data.name[, namevar])) 
      
      l <- list()
      for(i in 1:length(names))
        l[[names[i]]] <- names[i]
      
      updateSelectInput(session, "DRnameselector1",
                        choices = l,
                        select  = names[1])

    }
    
  })
  
  
  
  output$DRconc1 <- renderUI({
    
    drugs <- input$DRdrugselector1
    names <- input$DRnameselector1
    
    call <- A.data$auxiliary$passed.var
    back <- call$backgroundval
    additivevar <- call$additivevar
    namevar  <- call$namevar
    idvar    <- call$idvar
    identifier <- call$identifier
    drugvar  <- call$drugvar
    dosevar <- call$dosevar
    timevar <- call$timevar
    
    if("bootstrap" %in% class(A.data)){
      data.drug <- A.data$data$bs.raw.data[
        A.data$data$bs.raw.data[,drugvar] == drugs[1], ]
    }else{
      data.drug <- A.data$data$bc.data[
        A.data$data$bc.data[,drugvar] == drugs[1], ]
    }
    
    data <- data.drug[data.drug[, namevar] == names[1] & 
                        data.drug[, additivevar] != back, ]
    
    conc.names2 <- unique(data$type)
    
    
   
    checkboxGroupInput("DRconcnames", "Choose doses", 
                       conc.names2,selected = conc.names2)
  })
  
  
  
  
  DRmyres1 <- function(){
    input$DRPPI1
  }
  
  DRmypointsize1 <- function(){
    input$DRPointsize1
  }
  
  DRmyheight1 <- function(){
    input$DRheight1*300/2.54
  }
  
  DRmywidth1 <- function(){
    input$DRwidth1*300/2.54
  }
  
  output$DRGrowthCurves <- renderPlot({
    input$DRdrugselector1
    input$oldfilebutton
    input$DRanalysisButton
    
    if("doseResponseModel" %in% class(A.data)){
      

      if(input$DRnrows == 0){
        DRnrows = NULL
      }else{
        DRnrows = input$DRnrows
      }
      if(input$DRncols == 0){
        DRncols = NULL
      }else{
        DRncols = input$DRncols
      }
    
      
      plot.growthModel(x = A.data,
                       time.points.used = c("all"),
                       drugs  = input$DRdrugselector1,
                       names  = input$DRnameselector1,
                       conc.names = input$DRconcnames,
                       ylim   = NULL,
                       xlim   = NULL,
                       #  nrows  = DRnrows,
                       # ncols  = DRncols,
                       # xlab   = input$DRxlab1,
                       # ylab   = input$DRylab1,
                       scale  = TRUE,
                       #  main   = input$DRmain1,
                       absorbance.CI       = input$DRabsorbanceCI,
                       absorbance.CI.col   = "#C6C6C5",
                       absorbance.CI.alpha = 80,
                       bootstrap.conf = input$DRbootstrapcurves,
                       barcol         = "grey",
                       bar.height     = 1.3,
                      plotgrid = input$DRplotgrid,
                       grid.col = "#C6C6C5",
                       grid.lty = 1,
                       grid.lwd = 1,
                       bs.col   = "#C6C6C5",
                       bs.lty   = 1,
                       bs.lwd   = 0.5,
                       bs.alpha = 80, 
                       line.col = c("#4F6E9F", "#71965A", 
                                    "#9F9692", "#9D2441", 
                                    "#333333", "#662D91", 
                                    "#71DEC0", "#F7931E")[6],
                       line.lty = 1,
                       line.lwd = 1,
                       line.alpha = "",
                       plot.data  = TRUE,                            
                       col.by.identifier = TRUE,
                       col.points = c("#71965A", "#4F6E9F", 
                                      "#9F9692", "#9D2441", 
                                      "#333333", "#662D91", 
                                      "#71DEC0", "#F7931E")[3],         
                       pch.points.outlier = 4,
                       pch.points = 1,
                       plot.all = TRUE,
                       plot.data.all = FALSE,
                       line.col.all = "#9F9692",
                       line.lty.all = 1,
                       line.lwd.all = 0.8,
                       line.col.C0   = "#71965A",
                       line.col.GI50 = "#4F6E9F",
                       line.col.TGI  = "#9D2441",
                       line.col.LC48 = "#333333",
                       log        = "",
                       pdfit      = FALSE,
                       pdf.width  = 6.6929, 
                       pdf.height = 6.6929,
                       pointsize = 8)
      
    }else{
      return()
    }
    
  }, height = DRmyheight1, width = DRmywidth1, 
                                         pointsize = 8, 
                                         res = 300)
  
  
  observe({
    input$DRpngit1
    
    isolate({
      if("doseResponseModel" %in% class(A.data)){
        
        
        if(input$DRnrows == 0){
          DRnrows = NULL
        }else{
          DRnrows = input$DRnrows
        }
        if(input$DRncols == 0){
          DRncols = NULL
        }else{
          DRncols = input$DRncols
        }
        
        tempfile <- tempfile(fileext = ".png")
       
        png(tempfile, 
            height = input$DRheight1, 
            width = input$DRwidth1, units = "cm",
            pointsize = input$DRPointsize1,
            res = input$DRPPI1)
        
        plot.growthModel(x = A.data,
                         time.points.used = c("all"),
                         drugs  = input$DRdrugselector1,
                         names  = input$DRnameselector1,
                         conc.names = input$DRconcnames,
                         ylim   = NULL,
                         xlim   = NULL,
                         # nrows  = DRnrows,
                         # ncols  = DRncols,
                         # xlab   = input$DRxlab1,
                         # ylab   = input$DRylab1,
                         scale  = TRUE,
                         #  main   = input$DRmain1,
                         absorbance.CI       = input$DRabsorbanceCI,
                         absorbance.CI.col   = "#C6C6C5",
                         absorbance.CI.alpha = 80,
                         bootstrap.conf = input$DRbootstrapcurves,
                         barcol         = "grey",
                         bar.height     = 1.3,
                         plotgrid = input$DRplotgrid,
                         grid.col = "#C6C6C5",
                         grid.lty = 1,
                         grid.lwd = 1,
                         bs.col   = "#C6C6C5",
                         bs.lty   = 1,
                         bs.lwd   = 0.5,
                         bs.alpha = 80, 
                         line.col = c("#4F6E9F", "#71965A", 
                                      "#9F9692", "#9D2441", 
                                      "#333333", "#662D91", 
                                      "#71DEC0", "#F7931E")[6],
                         line.lty = 1,
                         line.lwd = 1,
                         line.alpha = "",
                         plot.data  = TRUE,                            
                         col.by.identifier = TRUE,
                         col.points = c("#71965A", "#4F6E9F", 
                                        "#9F9692", "#9D2441", 
                                        "#333333", "#662D91", 
                                        "#71DEC0", "#F7931E")[3],         
                         pch.points.outlier = 4,
                         pch.points = 1,
                         plot.all = TRUE,
                         plot.data.all = FALSE,
                         line.col.all = "#9F9692",
                         line.lty.all = 1,
                         line.lwd.all = 0.8,
                         line.col.C0   = "#71965A",
                         line.col.GI50 = "#4F6E9F",
                         line.col.TGI  = "#9D2441",
                         line.col.LC48 = "#333333",
                         log        = "",
                         pdfit      = FALSE,
                         pdf.width  = 6.6929, 
                         pdf.height = 6.6929,
                         pointsize = 8)
        dev.off()
        system(paste("open", tempfile))
      }else{
        return()
      }
      
    })
  })
  
  
  
  
}) 


#old.locale <- Sys.getlocale()