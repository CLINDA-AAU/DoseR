library(shiny)

shinyUI(pageWithSidebar(
  
  headerPanel("Dose response analysis"),
  sidebarPanel(
    
    textInput("oldfile", "Load existing project", value = "file"),
    actionButton("oldfilebutton", "Choose"),
    
    textInput("datafile", "File for storing results.", 
              value = file.path(getwd(), paste("AbsorbanceProject", Sys.Date(), sep = "_"))),
    actionButton("savefileButton", "Save Current project"),
    actionButton("changesavefileButton", "Choose filename"),
    
    tags$hr(),
    
   
    
    #########################################################
    ##
    ##  Program for creating metadata all input start with CMD
    ##
    ##########################################################
    
    conditionalPanel(
      condition = "input.conditionedPanels == '1. Meta data'",
      
    h2("Create the metadata"),
          
      textInput("CMDdbfdir", "Choose the directory storing the dBase files:", value = "directory"),
      actionButton("CMDdbfdirbutton", "Choose"),
      
      textInput("CMDprotocoldir", 
                "Choose the directory storing the protocol files:", value = "directory"),
      
      actionButton("CMDprotocoldirbutton", "Choose"),
      
      tags$hr(),
      
      textInput("CMDcorrectionname", "The name used for drug colour correction", value = "Control"),
      
      selectInput(
        "CMDmetadatafromfile", "Establish meta data",
        list("Choose" = "non",
             "Directly from file names" = "filenames", 
             "From an existing Excel sheet" = "excel" 
        )),
      
      conditionalPanel(
        condition = "input.CMDmetadatafromfile == 'excel'",  
        textInput("CMDmetadataexcel", "Choose the Excel file with the meta data:", value = "file"),
        actionButton("CMDmetadataexcelbutton", "Choose")
      ),
      
      
      conditionalPanel(
        condition = "input.CMDmetadatafromfile == 'filenames'",  
        textInput("CMDfileColnames", "List the identifiers in the filenames separated by a semicolon", 
                  value = "namevar;drugvar;protocolvar;identifier;timevar"),   
        
        radioButtons('CMDsep', 'Separator in the file names',
                     c("Comma (,)"=',',
                       "Semicolon (;)"=';',
                       "Underscore (_)"='_'),
                     'Semicolon (;)')
      ),
      
      conditionalPanel(
        condition = "input.CMDmetadatafromfile != 'non'", 
        radioButtons("CMDmetadataformat", "Is the metadata in wide or long format:",
                     list("Long" = "long",
                          "Wide" = "wide"))
      ),
      
      conditionalPanel(
        condition = "input.CMDmetadatafromfile == 'filenames'",
        
        
        checkboxInput("CMDeditcolnames", "Use custom names for variables", value=FALSE)        
      ),
      
      conditionalPanel(
        condition = "input.CMDmetadatafromfile == 'excel' | input.CMDeditcolnames == true",  
        textInput("CMDnamevar", "Choose the column name for identifying the cell line:", 
                  value = "name"),
        
        textInput("CMDdrugvar", "Choose the column name for identifying the drug:", 
                  value = "drug"),
        
        textInput("CMDprotocolvar", "Choose the column name for identifying the used protocol:", 
                  value = "R.protocol"),
        
        textInput("CMDtimevar", "Choose the column name for identifying time incubating with drug:", 
                  value = "time"),
        
        textInput("CMDidentifier", "Choose the name for the column identifying each run:", 
                  value = "identifier"),
        
        checkboxInput("CMDspecialidentifier", 
                      "Settings regarding identifier variable")  ,      
        
        conditionalPanel(
          condition = "input.CMDspecialidentifier == true", 
          radioButtons("CMDidentifiersep", "Separator to split the identifier into a replica variable",
                       c("Not separated",
                         "Comma (,)"=',',
                         "Semicolon (;)"=';',
                         "Underscore (_)"='_'),
                       'Underscore (_)'),
          
          selectInput("CMDdateformat", "Are the identifier variable in date format",
                      list("Not a date"    = "empty",
                           "ddmmyy (150812)" = "%d%m%y",
                           "m/d/y (02/27/92)"= "%m/%d/%y",
                           "ddmmmyyyy (2jan1960)" = "%d%b%Y"),
                      selected="Not a date")),
        
        
        
        textInput("CMDincubationvar", "Choose the name for the column containg the incubation time",  value = "incubation"),
        
        
        textInput("CMDdoublingvar", "Choose the column name for doubling time:", 
                  value = "not available")
      ),
      
      
      tags$hr(),
      
      selectInput("CMDcreateMetaDataunit", "What is the unit for the concentrations supplied in the protocols",
                  list("micro g/ml"    = "ug/ml",
                       "g/l"      = "g/l",
                       "mol/l"    = "mol/l"),
                  select="micro g/ml"),
      tags$hr(),
      checkboxInput("CMDextrametadatafromfile", 
                    "Add an extra metadata data sheet with cell line information"),
      conditionalPanel(
        condition = "input.CMDextrametadatafromfile == true",  
        textInput("CMDadditionalmetadata", "Choose the Excel file with the additional meta data:", value = "file"),
        actionButton("CMDadditionalmetadatabutton", "Choose")
      ),
      
      tags$hr(),
      checkboxInput("CMDupdate", 
                    "Do you want to create the metadata from scratch"),  
      actionButton("CMDCreateMetaDataButton", "Create metadata")
    
    
    ),
    
    #########################################################
    ##
    ##  Program for Editng the protocols PRC
    ##
    ##########################################################
    
    conditionalPanel(
      condition = "input.conditionedPanels == '2. Edit protocols'",
      
      h2("Edit the protocols"),
      
      conditionalPanel(
        condition = "input.panelEditProtocols == 'Create new protocol'",
        
        helpText("You can create a new protocol",
                 "associated with one or more dbf files"),
        
        textInput("PRCcreateProtocolName", "Choose a name for the new protocol:", 
                  value = "Name"),
       
        helpText("Note: Choose the dbf files",
                 "that should be associated with the new protocol"),
        
        selectInput("PRCCreatechoosedbf", "Choose the dbf file",
                    list(non = " "),
                    select="non",
                    multiple=FALSE),
          helpText("By pushing this button the selected dbf file will use",
                   "the protocol currently displayed to the right."),
        actionButton("PRCusethecreatedbutton", "Use the current protocol"),
        tags$hr(),
        helpText("You can create a protocol based on the criteria shown below"),
        textInput("PRCcreateProtocoldrug", "Name of the drug", 
                  value = "drug"),
        numericInput("PRCNumberofDoses", "Choose the number of doses:",    value = 16),
        numericInput("PRCNumberofControls", "Choose the number of Controls:", value =  6),
        numericInput("PRCNumberofBackgrounds", "Choose the number of Controls:", value =  6),
        numericInput("PRCMaxDose" , "Choose maximum dose:",  value = 60),
        numericInput("PRCdoseFold", "Choose dilution fold:", value =  2),
        checkboxInput("PRCCreateRemoveEdges", "Should the edges be used")    
      )
      
      
    
#       conditionalPanel(
#         condition = "input.panelEditProtocols == 'Edit a protocol'",
#         helpText("If one of the dbf files is associated with an error",
#                  "or change in the setup you",
#                  "can create a new protocol for that specific dbf file"),
#         
#           helpText("By pushing this button an Excel spreadsheet containg ",
#                    "the current protocol will be opened. The protocol have been renamed",
#                    "and the dbf file renamed accordingly"),
#           actionButton("PRCusetheEditdbutton", "Edit the protocol")  
#       )
     ),
    
    
    #########################################################
    ##
    ##  Program for creating metadata all input start with DBF
    ##
    ##########################################################
    
    conditionalPanel(
      condition = "input.conditionedPanels == '3. Read data'",
      
      h2("Load the dbf files"),
      
      conditionalPanel(
        condition = "input.panelReadDbffiles == 'Print' | input.panelReadDbffiles == 'Call to R'| input.panelReadDbffiles == 'dbf viewer'| input.panelReadDbffiles == 'Protocol viewer'",
        
        
        h3("dbf file information"),
        textInput("DBFdiscardlines", "Should any lines of the dbf files be discarded:", 
                  value = "1;2"),
        
        textInput("DBFabsorbanceid", "Name of the column containing the absorbance values", value = "M1"),
        textInput("DBFwellid", "Name of the column containing well id's", value = "WELLNUM"),
        tags$hr(),
        h3("Protocol information"),
        textInput("DBFremoverows", "Should any rows be discarded", value = "A;H"),
        textInput("DBFremovecols", "Should any columns be discarded", value = "1;12"),
        
        textInput("DBFdosevar", "Name of the column containg dosage information", value = "Concentration"),
        textInput("DBFadditivevar", "Name of the column containing information about containment", value = "Additive"),
        
        textInput("DBFcontrolval", "Value of control variables", value = "Control"),
        textInput("DBFbackgroundval", "Value of background variables", value = "Background"),
        textInput("DBFmistakeval", "Value used to indicate mistakes", value = "X"),
        
        tags$hr(),
        checkboxInput("DBFupdate", 
                      "Do you want to load the data from scratch.", value=FALSE), 
        helpText("Note: This will force each element of the analysis" ,
                 "to run from scratch the entire analysis to run from scratch"),
        actionButton("DBFreadDBFDataButton", "Load the dbf files") 
        
      ),
      
      conditionalPanel(
        condition = "input.panelReadDbffiles == 'Change molar masses'",
        
        h3("Change Molar masses"),
        helpText("Note: Choose one of the available dbf files"),
        selectInput("Molareditdrug", "Select a drug",
                    list(non = " "),
                    select="non"),
        
        numericInput("Molareditmass", "Specify the molar mass:", value =  1000),
        actionButton("Molareditdrugbutton", "Edit the molar mass") 
      ) 
    ),
    
    #########################################################
    ##
    ##  Drug colour correction program
    ##
    ##########################################################
    
    conditionalPanel(
      condition = "input.conditionedPanels == '4. Drug color correction'",

    h2("Drug color correction"),
  
      helpText("How many standard deviations does a point have to be from the",
               "mean in order to be classified as an outlier"),
      numericInput("DCCoutlierTest", "number of standard deviations", value =  3),
      helpText("How many times should the outlier detection routine be run"),
      numericInput("DCCoutlierIter",
                   "number of outlier detection iterations", 
                   value =  2),
      
      helpText("Should the optimisation routine be restricted to",
               "only fit positive absorbance values",
               "It is recommended to use the unrestricted fit",
               "and afterwords adjust all absorbance measures below a cut point",
               "to that cut point"),
      selectInput(
        "DCCpararametrisation", "The parametrisation to be used",
        list("Unrestricted"       = "unrestricted",
             "Restricted"   = "restricted"
        )),
      tags$hr(),
      checkboxInput("DCCtech", "Technical adjustments"),
      conditionalPanel(
        condition = "input.DCCtech == true", 
        helpText("What contrast should for the estimation of ",
                 "the absorbance values. The sum contrast is highly recommended",
                 "since this corresponds to an average."),
        selectInput(
          "DCCcontr", "The contrast to be used",
          list("Sum"       = "sum",
               "Helmert"   = "helmert", 
               "Treatment" = "treatment" 
          )),
        helpText("Fit the weights of the heteroscedastic variance"),
        checkboxInput("DCCweights", "Fit weights", value=TRUE),
        helpText("Use only the fitted absorbance vales for estimating",
                 "the weights of the heteroscedastic variance.",
                 "The background is thus not used when the weights are fitted."),
        checkboxInput("DCCfittedA", "Discard background", value=FALSE),
        helpText("What conversion criterin should be used"),
        numericInput("DCCvarpowerMin",
                     "Conversion criterion", 
                     value =  10e-4),
        helpText("Maximum number of iterations aloud"),
        numericInput("DCCvarpowerIter",
                     "Conversion criterion", 
                     value =  50)
      ),
      tags$hr(),
      checkboxInput("DCCupdate", 
                    "Do you want to load the data from scratch.", value=FALSE), 
      helpText("Note: This will force each element of the analysis" ,
               "to run from scratch"),
      actionButton("DCCanalysisButton", "Perform the drug colour correction")  
    ), 
    
    
    #########################################################
    ##
    ##  Normalization program
    ##
    ##########################################################
    conditionalPanel(
      condition = "input.conditionedPanels == '5. Pre-processing'",
    
    h2("Normalisation"),

      conditionalPanel(
        condition = "input.BGNconditionedPanels == 'Print' | input.BGNconditionedPanels == 'Call to R'",
        
      helpText("How many standard deviations does a point have to be from the",
               "mean in order to be classified as an outlier"),
      numericInput("BGNoutlierTest", "number of standard deviations", value =  3),
      helpText("How many times should the outlier detection routine be run"),
      numericInput("BGNoutlierIter",
                   "number of outlier detection iterations", 
                   value =  2),
      
      
      
      helpText("Should the optimisation routine be restricted to",
               "only fit positive absorbance values",
               "It is recommended to use the unrestricted fit",
               "and afterwords adjust all absorbance measures below a cut point",
               "to that cut point"),
      selectInput(
        "BGNpararametrisation", "The parametrisation to be used",
        list("Unrestricted"       = "unrestricted",
             "Restricted"   = "restricted"
        )),
      tags$hr(),
      checkboxInput("BGNtech", "Technical adjustments"),
      conditionalPanel(
        condition = "input.BGNtech == true", 
        helpText("What contrast should for the estimation of ",
                 "the absorbance values. The sum contrast is highly recommended",
                 "since this corresponds to an average."),
        selectInput(
          "BGNcontr", "The contrast to be used",
          list("Sum"       = "sum",
               "Helmert"   = "helmert", 
               "Treatment" = "treatment" 
          )),
        helpText("Fit the weights of the heteroscedastic variance"),
        checkboxInput("BGNweights", "Fit weights", value=TRUE),
        helpText("Use only the fitted absorbance vales for estimating",
                 "the weights of the heteroscedastic variance.",
                 "The background is thus not used when the weights are fitted."),
        checkboxInput("BGNfittedA", "Discard background", value=FALSE),
        helpText("What conversion criterin should be used"),
        numericInput("BGNvarpowerMin",
                     "Conversion criterion", 
                     value =  10e-4),
        helpText("Maximum number of iterations aloud"),
        numericInput("BGNvarpowerIter",
                     "Conversion criterion", 
                     value =  50)
      ),
      tags$hr(),
      checkboxInput("BGNupdate", 
                    "Do you want to load the data from scratch.", value=FALSE), 
      helpText("Note: This will force each element of the analysis" ,
               "to run from scratch"),
      actionButton("BGNanalysisButton", "Perform the normalisation") 
      ),
      
      conditionalPanel(
        condition = "input.BGNconditionedPanels == 'Normalised data'",
        
      
      selectInput("BGNdrugselector1", "Choose a drug",
                  list(non = " "),
                  select="non"),
      selectInput("BGNnameselector1", "Choose a cell line",
                  list(non = " "),
                  select="non"),
      
      numericInput("BGNwidth1", "Choose the width (cm)",    
                   value = 17),
      numericInput("BGNheight1", "Choose the height (cm)",    
                   value = 21),
      numericInput("BGNPPI1", "Set the PPI",    
                   value = 300),
      numericInput("BGNPointsize1", "Set the point size",    
                   value = 8),
      textInput("BGNunit1", "The unit used for the x-axis",    
                value = "mol/l"),
      
      selectInput("BGNlogselector1", "log transform the x-axis",
                  list("log base 10" = "log10",
                       "log base 2"  = "log2",
                       "log base e"  = "log",
                       "no log"      = ""),
                  select="log base 10"),
      
      
      helpText("You can edit the settings of the plot"),
      checkboxInput("BGNeditplot", "Edit the plot"),
      conditionalPanel(
        condition = "input.BGNeditplot == true",
        
        helpText("Edit titles"),
        textInput("BGNtitle1", "Title of the figure", 
                  value = "Standarad"),
        textInput("BGNtitle2", "Title of the rows",  
                  value = "Raw;Color Corrected;Simple Correction;Model Correction"),
        textInput("BGNmain1", "Title above each plot", 
                  value = "Standard"),
        
        
        helpText("Edit labels"),
        textInput("BGNylab1", "The name of the y-axis", 
                  value = "Absorbance"),
        textInput("BGNxlab1", "The name of the x-axis", 
                  value = "Standard")
      )
    ),
      
      conditionalPanel(
        condition = "input.BGNconditionedPanels == 'Model check'",
      
        selectInput("BGNdrugselector2", "Choose a drug",
                    list(non = " "),
                    select="non"),
        selectInput("BGNnameselector2", "Choose a cell line",
                    list(non = " "),
                    select="non")
      )
    ),
    
    #########################################################
    ##
    ##  Bootstrap Routine
    ##
    ##########################################################
    conditionalPanel(
      condition = "input.conditionedPanels == '6. Bootstrap'",
      
    h2("Bootstrap"),
   
      helpText("How many bootstrap samples should be generated"),
      numericInput("Bootnsamples", "Bootstrap samples", value =  10),
      
      helpText("How many iterations is the maximum allowed"),
      numericInput("Bootnsamplesmax", "Bootstrap samples", value =  15),
      
      selectInput(
        "Boottype", "The bootstrap procedure to be used",
        list("Parametric"       = "parametric",
             "Residual"   = "residual"
        )),
    
      tags$hr(),
      checkboxInput("Bootupdate", 
                    "Do you want to create bootstrap samples from scratch.", 
                    value=FALSE), 
      actionButton("BootanalysisButton", "Generate bootstrap samples")  
    ),
    
    
    #########################################################
    ##
    ##  Dose response models 
    ##
    ##########################################################
    
    conditionalPanel(
      condition = "input.conditionedPanels == '7. Dose Response'",
      
      h2("Dose response models"),
      
      conditionalPanel(
        condition = "input.DRconditionedPanels == 'Print' | input.DRconditionedPanels == 'Call to R'",
        
        
        helpText("Which models should be estimated"),
        checkboxGroupInput("DRmodels", "Selected models", 
                           list("G", "R", "D", "DG", "RG"),
                           selected = c("G", "R", "D")),
        
        
        helpText("Should the estimation of the G model be restrictd.",
                 "This is used when the doubling time for the drug",
                 "treated cell is not allowed to grow faster than",
                 "the same cell line untreated"),
        
        selectInput(
          "DRparametrisation", "The parametrisation to be used",
          list("Restricted"       = "restricted",
               "Non restricted"   = "nonrestricted"
          )),
        
        helpText("Set the threshold.",
                 "All absorbance measures below this value is set to this value."),
        numericInput("DRcut", "Threshold", value =  0.025),
        
        helpText("Do the supplied data include an estimate of the doubling time.",
                 "For dose response experiments using only one plate.",
                 "The G model can be calculated through the R using a prespecified",
                 "doubling time."),
        textInput("DRdoublingvar", "Doubling variable name", value = "No"),
        
        
        helpText("What unit shall the summary statistics be calculated"),
        textInput("DRunit", "The unit",value = "mol/l"),
        
        helpText("The log transformation of the concentrations given in the protocols"),
        selectInput("DRlogfrom", "prior log transform",
                    list("log base 10" = "log10",
                         "log base 2"  = "log2",
                         "log base e"  = "log",
                         "no log"      = "nolog"),
                    select="no log"),
        helpText("The log transformation used for calculation of summary statistics"),
        
        selectInput("DRlogto", "Post log transform",
                    list("log base 10" = "log10",
                         "log base 2"  = "log2",
                         "log base e"  = "log",
                         "no log"      = "nolog"),
                    select="log base 10"),
        
        
        helpText("The LC50 value for the G model is based on a number of hours."),
        numericInput("DRt", "LC50 halving time", value =  48),
        
        helpText("The AUCq shoul be calculated with the following value for q"),
        numericInput("DRAUCq", "AUC q", value =  0),
        
        tags$hr(),
        checkboxInput("DRupdate", 
                      "Do you want to create estimate the models from scratch.", 
                      value=FALSE), 
        actionButton("DRanalysisButton", "Estimate the models")           
        
      ),
      
      
      conditionalPanel(
        condition = "input.DRconditionedPanels == 'Growth curves'",
        
        
        selectInput("DRdrugselector1", "Choose a drug",
                    list(non = " "),
                    select="non"),
        selectInput("DRnameselector1", "Choose a cell line",
                    list(non = " "),
                    select="non"),
        
        helpText("You can deselect some of the doses seen in the plot"),
        checkboxInput("DReditconc", "Edit doses"),
        conditionalPanel(
          condition = "input.DReditconc == true",
          uiOutput("DRconc1")
        ),
        
        numericInput("DRwidth1", "Choose the width (cm)",    
                     value = 17),
        numericInput("DRheight1", "Choose the height (cm)",    
                     value = 21),
        numericInput("DRPPI1", "Set the PPI",    
                     value = 300),
        numericInput("DRPointsize1", "Set the point size",    
                     value = 8),
        
        
        # col = 1:8,
        
        
        helpText("You can edit the settings of the plot"),
        checkboxInput("DReditplot", "Edit the plot"),
        conditionalPanel(
          condition = "input.DReditplot == true",
          
          helpText("Edit titles"),
          textInput("DRmain1", "Title of the figure", 
                    value = "Dose Response Growth Curves"),
          
          
          
          helpText("Edit labels"),
          textInput("DRylab1", "The name of the y-axis", 
                    value = "Absorbance"),
          textInput("DRxlab1", "The name of the x-axis", 
                    value = "Time (Hours)"),
          
          checkboxInput("DRbootstrapcurves", "Show bootstrap curves"),
          checkboxInput("DRabsorbanceCI", "Show CI of absorbance"),
          checkboxInput("DRplotgrid", "Show grid"),
          
          numericInput("DRnrows", "Set the number of rows",    
                       value = 0),
          numericInput("DRncols", "Set the number of cols",    
                       value = 0)
        )
        
      )
    ),
    
    progressInit()
  ),
  
  
  mainPanel(    
    
#     conditionalPanel(
#       condition = "input.CMDcreateMetaData == false",
#       tabsetPanel(
#         tabPanel("Current inputs",  verbatimTextOutput("current.inputs"))
#       )
#     ),
    
    h2("Analysis Flow"),
    tabsetPanel(
      tabPanel("Info"),
      tabPanel("1. Meta data"),
      tabPanel("2. Edit protocols"),
      tabPanel("3. Read data"),
      tabPanel("4. Drug color correction"),
      tabPanel("5. Pre-processing"),
      tabPanel("6. Bootstrap"),
      tabPanel("7. Dose Response"),
      id = "conditionedPanels"
    ),

    
    conditionalPanel(
      condition = "input.conditionedPanels == 'Info'",
      tabsetPanel(
        tabPanel("Welcome"),
        tabPanel("Analysis flow", verbatimTextOutput("AnalysisFlow")),
        tabPanel("Citation", verbatimTextOutput("Citation")),
        id="panelinfo")
    ),
   
    
    conditionalPanel(
      condition = "input.conditionedPanels == '1. Meta data'",
      tabsetPanel(
        tabPanel("Print", verbatimTextOutput("CMDprint")),
        tabPanel("Call to R",  verbatimTextOutput("CMDcreateMetaDataCall")),
        tabPanel("Example", tableOutput("CMDexample")),
        tabPanel("Metadata", tableOutput("CMDmetadataTable")),
        tabPanel("Metadata Correction", tableOutput("CMDmetaCorrection")),
        tabPanel("Metadata Additional", tableOutput("CMDmetaadditional")),     
        id="panelCMDcreateMetaData",
        selected = "Example")
    ),
    
    
    conditionalPanel(
      condition = "input.conditionedPanels == '2. Edit protocols'",
      tabsetPanel(
        tabPanel("Create new protocol",   
                 helpText("Note: You can choose one the existing protocols to use as base",
                          "that should be associated with the new protocol"),
                 
                 selectInput("PRCcreatechooseProtocol", "Choose a protocol",
                             list(non = " "),
                             select="non"),
                 
                 actionButton("PRCCreateeditProtocolButton",  "Edit the protocol") ,
                 actionButton("PRCCreatereloadProtocolButton", "Reload the protocol") ,
                 
                 tableOutput("PRCCurrentSetup"),
                 tableOutput("PRCCurrentConc"),
                 
                 tableOutput("PRCdbfeksample")),
        tabPanel("Track changes"),
        id="panelEditProtocols"
      )
    ),
    
    
    
    conditionalPanel(
      condition = "input.conditionedPanels == '3. Read data'",
      tabsetPanel(
        tabPanel("Print", verbatimTextOutput("DBFreadDBFDataPrint")),
        tabPanel("Print2", verbatimTextOutput("DBFreadDBFDataexsfilep")),
        tabPanel("Call to R",  verbatimTextOutput("DBFreadDBFDataCall")),
        tabPanel("dbf viewer",
                   helpText("Note: You can choose one of the available dbf files",
                            "and thereby observe which column contains the absorbance measures",
                            "and which column contains the well id."),
                   helpText("You can also see if the files contains columns that should be discarded.",
                            "If empty columns are not removed the program will crash." ),
                   selectInput("DBFchoosedbf", "Choose a dbf file",
                               list(non = " "),
                               select="non"),
                   # p(textOutput("DBFprotocolselected"))
                   tableOutput("DBFeksample")
        ),
        tabPanel("Protocol viewer",
                   helpText("You can choose one of the available protocols",
                            "and thereby observe which column contains the ..."),
                   selectInput("DBFchooseprotocol", "Choose a protocol",
                               list(non = " "),
                               select="non"),
                   tableOutput("DBFprotocolSetup"),
                   tableOutput("DBFprotocolConc")
        ),
        tabPanel("Raw data",  dataTableOutput("dataRaw")),
        tabPanel("Change molar masses", tableOutput("molarmassTable")),
        id="panelReadDbffiles"
      )
    ),
    
    conditionalPanel(
      condition = "input.conditionedPanels == '4. Drug color correction'",
      tabsetPanel(
        tabPanel("Print", verbatimTextOutput("DCCPrint")),
        tabPanel("Call to R",  verbatimTextOutput("DCCcall")),
        tabPanel("Drug Colour", 
                 selectInput("DCCdrugselector1", "Choose a drug",
                             list(non = " "),
                             select="non"),
                 plotOutput("DCCcolorCorrection")),
        tabPanel("Model check", 
                 selectInput("DCCdrugselector2", "Choose a drug",
                             list(non = " "),
                             select="non"),
                 plotOutput("DCCcolorCorrectionRaw")),
        id = "DCCconditionedPanels"
      )
    ),
    
    
    conditionalPanel(
      condition = "input.conditionedPanels == '5. Pre-processing'",
      tabsetPanel(
        tabPanel("Print", verbatimTextOutput("BGNPrint")),
        tabPanel("Call to R",  verbatimTextOutput("BGNcall")),
        tabPanel("Normalised data",
                 actionButton("BGNpdfit1", "Open as pdf"),
                 actionButton("BGNpngit1", "Open as png"),
                 plotOutput("BGNnormalisedData", height="auto")),
        tabPanel("Model check", 
                 plotOutput("BGNModelcheck")),
        id = "BGNconditionedPanels"
      )
    ),
    
    conditionalPanel(
      condition = "input.conditionedPanels == '6. Bootstrap'",
      tabsetPanel(
        tabPanel("Print", verbatimTextOutput("BootPrint")),
        tabPanel("Call to R",  verbatimTextOutput("bootcall")),
        tabPanel("Bootstrap data",  dataTableOutput("dataBoot")),
       id = "BootstrappedconditionedPanels"
      )
    ),
    
    
    conditionalPanel(
      condition = "input.conditionedPanels == '7. Dose Response'",
      tabsetPanel(
        tabPanel("Print", verbatimTextOutput("DRPrint")),
        tabPanel("Call to R",  verbatimTextOutput("DRCall")),
        tabPanel("Dose response data",  dataTableOutput("DRData")),
        tabPanel("Growth curves",
                 actionButton("DRpdfit1", "Open as pdf"),
                 actionButton("DRpngit1", "Open as png"),                 
                 plotOutput("DRGrowthCurves", height="auto")),
        id = "DRconditionedPanels"
      )
    )
  )
)
)