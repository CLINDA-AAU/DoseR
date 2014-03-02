
plotbgModelresid <- function(A.data,
                              names = NULL,
                              times = NULL,
                              drugs = NULL,
                              pdfit = FALSE,
                              file= file.path(getwd(), "resid.plots.pdf"),
                              pointsize = 12,
                              xlab = "Fitted values",
                              ylab = "Residuals",
                              main = NULL,
                              col = 1:8, 
                              width = 3.34, height =2.23097112){
  
  old.par <- par(no.readonly = TRUE)  
  on.exit(par(old.par))
  
  main.orig <- main
  
  namevar <- eval(A.data$call$readDBFData$namevar)
  drugvar <- eval(A.data$call$readDBFData$drugvar)
  timevar <- eval(A.data$call$readDBFData$timevar)
  
  if(is.null(names))
    names <- unique(A.data$data$bc.mean[, namevar])
  if(is.null(drugs))
    drugs <- unique(A.data$data$bc.mean[, drugvar])
  if(is.null(times))
    times <- unique(A.data$data$bc.mean[, timevar])
  
  if(pdfit)
    pdf(file, pointsize = pointsize, width = width, height = height)
  
  
  wh <-  A.data$data$bc.mean[, namevar] %in% names &
    A.data$data$bc.mean[, drugvar] %in% drugs &
    A.data$data$bc.mean[, timevar] %in% times
  
  levels <- unique(A.data$data$bc.mean[wh, "level"])
  
  levels <- levels[
    levels %in% names(A.data$fits$bgModel[ 
      rep(1:4, length(A.data$fits$bgModel)/4) == 1])]
  
  for(level.iter in levels){
    power <- A.data$fits$bgModel[[level.iter]]$modelStruct$varStruct[1]
   
    par(mfrow = c(1,2))
    #main <- gsub(c("Doxorubicin"), c("Adriamycin"), level.iter)
    if(is.null(main.orig)){
      main <- level.iter
#       main <- gsub(c(48), c(49), main)
#       main <- gsub(c(36), c(37), main)
#       main <- gsub(c(24), c(25), main)
#       main <- gsub(c(12), c(13), main)
#       main <- gsub(c( 0), c( 1), main)
    }
    
    plot(resid(A.data$fits$bgModel[[level.iter]])#,
         #center = FALSE)   
         ~
           fitted(A.data$fits$bgModel[[level.iter]]),
         main = main,
         col = col[A.data$fits$bgModel[[paste(level.iter, "col", sep = ":")]]],
         pch = A.data$fits$bgModel[[paste(level.iter, "pch", sep = ":")]],
         xlab = xlab,
         ylab = ylab)
    abline(0, 0)
    
    
    
    plot(resid(A.data$fits$bgModel[[level.iter]]) / fitted(A.data$fits$bgModel[[level.iter]])^power#,
         #center = FALSE)   
         ~
           fitted(A.data$fits$bgModel[[level.iter]]),
         main = main,
         col = col[A.data$fits$bgModel[[paste(level.iter, "col", sep = ":")]]],
         pch = A.data$fits$bgModel[[paste(level.iter, "pch", sep = ":")]],
         xlab = xlab,
         ylab = ylab)
    abline(0, 0)
    #grid()
    
  }
  if(pdfit)
    dev.off()
  
  invisible()
}


plotbgModel <- function(A.data, names = NULL, drugs = NULL, times = NULL,
                         pdfit = FALSE,
                         unit = "mol/l",
                         logfun = "log10",
                         plots = c("Raw", "Color Corrected", 
                                   "Simple Correction", "Model Correction"),
                         ylab = "Absorbance",
                         xlab = NULL,
                         main = NULL,
                         title.1 = NULL,
                         title.2 = c("Raw", "Color Corrected", 
                                     "Simple Correction", "Model Correction"),
                         figure.output = file.path(getwd(), "Normalisation"), 
                         col = c("#4F6E9F", "#71965A", 
                                 "#9F9692", "#9D2441", 
                                 "#333333", "#662D91", 
                                 "#71DEC0", "#F7931E"),
                         width = 6.69,
                         height = NULL,
                         pointsize = 5,
                         cex = 1,
                        set.par = TRUE,
                         ...) {
  
  if(set.par){
    old.par <- par(no.readonly = TRUE)
    on.exit(par(old.par))
  }
  plots.orig <- plots
  
  if(is.null(xlab) & logfun == "log10")
    xlab <- as.expression(substitute(
      paste("Concentration", " ",
            log[10],"(",m,")", sep = ""), list(m = unit)))
  
  if(is.null(xlab) & logfun == "log2")
    xlab <- as.expression(substitute(
      paste("Concentration", " ",
            log[2],"(",m,")", sep = ""), list(m = unit)))
  
  if(is.null(xlab) & logfun == "log")
    xlab <- as.expression(substitute(
      paste("Concentration", " ",
            log,"(",m,")", sep = ""), list(m = unit)))
  
  if(is.null(xlab) & logfun == "")
    xlab <- as.expression(substitute(
      paste("Concentration", " ",
            "(", m,")", sep = ""), list(m = unit)))
  
  plots2 <- c("absorbance.nc", "absorbance", "NB", "BC2")
  
  names(plots2) <- c("Raw", "Color Corrected", 
                     "Simple Correction", "Model Correction")
  
  plots <- plots2[plots]
  
  
  dir.create(figure.output, showWarnings=FALSE, 
             recursive=TRUE, mode="0777")
  
  call <- A.data$auxiliary$passed.var
  namevar     <- call$namevar
  idvar       <- call$idvar
  drugvar     <- call$drugvar
  dosevar     <- call$dosevar
  timevar     <- call$timevar
  backgroundval <- call$backgroundval
  controlval    <- call$controlval
  additivevar <- call$additivevar
  incubationsvar <- call$incubationvar
  identifier <- call$identifier
  
  if(is.null(drugs))
    drugs <- unique(A.data$data$bc.data[, drugvar]) 
  
  if(is.null(times))
    times <- unique(A.data$data$bc.data[, timevar]) 
  
  if(is.null(names))
    names <- unique(A.data$data$bc.data[, namevar]) 
  
  for(drug.iter in drugs) {
    
    temp.drug <- A.data$data$bc.data[A.data$data$bc.data[, drugvar] == drug.iter &
                                       A.data$data$bc.data[, timevar] %in% times, ]
    
    
    
    
    for(name.iter in intersect(unique(temp.drug[, namevar]), names)){
      
      
      temp.name <- temp.drug[temp.drug[, namevar] == name.iter, ]
      
      plots.n <- vector()
      for(i in 1:length(plots)){
        plots.n[i] <- (plots[i] %in% colnames(temp.name))
        if(plots.n[i])
          plots.n[i] <- !all(is.na(temp.name[, plots[i]]))
      }
      plots.new <- plots[plots.n]
      if(!("absorbance.nc" %in% plots.new) & "absorbance.nc" %in% plots2){
        names(plots.new) 
      }
      
      levels <- unique(temp.name[, "level"])
      
      if(pdfit) {
        if(is.null(height)){
          height2 = (width) / length(plots.new) * length(levels)
        }else{
          height2 <- height
        }
        
        pdf(file.path(figure.output, 
                      paste("background_correction", 
                            drug.iter, "_", name.iter, ".pdf",  sep = "")),
            width = width, height = height2,
            pointsize = pointsize)
      }
     
      bar.height = 10
      
      if(set.par){
        par(mfrow = c(length(levels), length(plots.new)))
        
        par(oma = c(1, 2, 4.2, 1))
        par(mar = c(4, 4, 4, 2))
      }
      ylim <- range(temp.name[, plots.new], na.rm = TRUE)
      
      for(level.iter in levels){
        
        n.data <- temp.name[temp.name[,"level"] == level.iter  ,]
        time.iter <- unique(n.data[, timevar])
        
        
        if(length(title.2) == 1){
          title.2.2 <- rep(title.2, length(plots))
        }else{
          title.2.2 <- title.2
        }
        
        title.2.new <- title.2.2[plots.n]
        
        if(!("absorbance.nc" %in% plots.new) & "absorbance.nc" %in% plots2){
          CC  <- names(plots.new)  == "Color Corrected"
          Raw <- names(plots) == "Raw"
          
          title.2.new[CC] <- title.2.2[Raw]
        }
        
        
        n.data[,dosevar][n.data[,dosevar] == 0] <-
          min(n.data[,dosevar][n.data[,dosevar]>0])/2
        
        
        #xlab <- "Concentration"
        
        pch <- rep(1, nrow(n.data))
        
        pch[n.data[, additivevar] %in% c(backgroundval, controlval)] <-
          ifelse(n.data[, additivevar][n.data[, additivevar] %in% c(backgroundval, controlval)] == controlval,
                 2, 3)
        
        pch[n.data$outlier == 1] <- 4 
        
        n.data[, dosevar] <- 
          concConvert(n.data[, dosevar], 
                      mol.mass = A.data$auxiliary$mol.data[drug.iter, "mol.mass"],       
                      from = A.data$auxiliary$passed.var$unit, 
                      to        = unit,
                      logfun.to = logfun)
        
        for(j in 1:length(plots.new)){
          
          if(is.null(main)){
            t <-  time.iter + n.data[1, incubationsvar] / 2
            t1 <- ifelse(t == 1, "hour", "hours") 
            mains <- paste("Time", t, t1)
          }
          
          log <- ifelse(grepl("log", logfun), "", "x")
          
          plot(n.data[, plots.new[j]] ~ n.data[, dosevar], log = log,
               col  = col[as.numeric(as.factor(n.data$plate))], las = 1,
               pch  = pch,
               ...,
               main = ifelse(is.null(main), mains, main), 
               xlab = xlab,
               ylab = ylab,
               ylim = ylim,
               cex = cex)
          
          if(level.iter == levels[1])
            mtext(title.2.new[j], line= 3, cex = 1.2*cex)
          if(j == 1){
            legend("bottomleft", pch = 1, col = col[1:length(n.data$plate)],
                   legend = levels(as.factor(n.data[, identifier])), 
                   bty = "n", cex = cex)
            legend("bottomright", pch = 1:3, 
                   legend = c("Drug", controlval, backgroundval), 
                   bty = "n", cex = cex)
          }
        } 
      }
      is.null(title.1)
      title.1 = paste(name.iter, " (", drug.iter,")", sep = "")
      mtext(title.1, outer = TRUE, line = 2.2, cex = 1.5*cex)
      if(pdfit)
        dev.off()
    }    
  }
  invisible()
}


plotdrugColorCorrection <- 
  function(A.data, file = file.path(getwd(), "drugColorCorrection.pdf") ,
           drugs = NULL, pdfit = FALSE,
           width = 7, height = 3.5, 
           type.plot = c("output", "raw"),
           xlab = NULL, ylab = NULL,
           main = NULL, 
           xlim = NULL,
           ylim = NULL,
           col = c("#71965A", "#4F6E9F", "#9F9692", "#9D2441",
                   "#333333", "#662D91", "#71DEC0", "#F7931E"),
           col.plate = TRUE,
           set.par = FALSE,
           ...){
    
    if(set.par){
      old.par <- par(no.readonly = TRUE)
      on.exit(par(old.par))
    }
    
    ylim.orig <- ylim
    xlim.orig <- xlim
    
    
    call        <- A.data$auxiliary$passed.var
    A           <- "absorbance"
    B           <- eval(call$backgroundval)
    idvar       <- eval(call$idvar)
    drugvar     <- eval(call$drugvar)
    namevar     <- eval(call$namevar)
    timevar     <- eval(call$timevar)
    dosevar     <- eval(call$dosevar)
    additivevar <- eval(call$additivevar)
    backgroundval <- eval(call$backgroundval)
    controlval  <- eval(call$controlval)
    plateset    <- "plateset"
    plate       <- "plate"
    type        <- "new.type"
    record      <- A.data$auxiliary$record
    data.file   <- eval(call$data.file)
    
    correction.data <- A.data$drug.color.correct$correction.data
    
    xlab.orig <- xlab
    ylab.orig <- ylab
    main.orig <- main
    
    if("output" %in% type.plot) {
      data <- correction.data$data$bc.data
      # drug <- "drug"
      
      if(is.null(drugs)){
        drugs <- unique(data[,drugvar])
      }else{
        drugs <- intersect(unique(data[,drugvar]), drugs)
      }
      if(pdfit)
        pdf(file, width = width, height = height)
      for(drugColor in drugs){
        
        data.color <- data[data[,timevar] == 0 & data[, backgroundval] != 1 &
                             data[, drugvar] == drugColor, ]
        data.color$up <- data.color$BC + 1.96 * data.color$Std.Error
        data.color$lo <- data.color$BC - 1.96 * data.color$Std.Error
        
        formula <- as.formula(paste("cbind(BC, lo, up ) ~",  dosevar))
        col.C <-
          aggregate(formula, FUN = mean,
                    data = data.color)
        
        col.C$BC2 <- col.C$BC -
          mean(data.color[data.color[, additivevar] == controlval , "BC"], na.rm = TRUE)
        
        #plot(col.C$Concentration, col.C$BC2, log= "x", main = drugColor,
        #     xlab = "Concentration",
        #     ylab = "Relative to Control")
        if(is.null(xlab.orig))
          xlab <- expression(paste("Concentration", " ",
                                   log[10],"(", m,"mol/ml)", sep = ""))
        
        if(is.null(main.orig))
          main <- paste(drugColor, "Background Colour")
        
        if(is.null(ylab.orig))
          ylab <- "Measured absorbance"
        
        col.C$Conc <- (col.C[, dosevar])
        col.C$Conc[col.C$Conc == 0] <-
          min(col.C$Conc[col.C$Conc != 0])/2
        
        col.C$Conc <- DoseR:::concConvert(col.C$Conc, 
                                  mol.mass = A.data$auxiliary$mol.data[drugColor, 2], 
                                  logfun.to = "log10") 
        
        if(is.null(xlim.orig))
          xlim  <- range(col.C$Conc)
        
        if(is.null(ylim.orig))
          ylim  <- range(c(col.C$lo, col.C$up))
  
        plot(range(col.C$Conc), range(c(col.C$lo, col.C$up)),
             ylim = ylim,
             xlim = xlim,
             main = main,
             xlab = xlab, type = "n",
             col = col[1],
             ylab = ylab, ...)
        
        for(c in col.C$Conc)
          segments(c, col.C$lo[col.C$Conc == c], c, 
                   col.C$up[col.C$Conc == c], col =col[1])
        #ylim = c(0,0.08))
        lines(col.C$Conc, col.C$BC, col = col[1])
        points(col.C$Conc, col.C$BC, pch = 
                 c(15, rep(19, length(col.C$Conc))), col = col[1])
        
      } 
      if(pdfit)
        dev.off()
    }
    
    
    
    
    if("raw" %in% type.plot) {
      correction.data <- A.data$drug.color.correct$correction.data
      if(is.null(drugs))
        drugs <- unique(correction.data$data$bc.data[, drugvar]) 
      
      if(pdfit)
        pdf(file, width = 14, height = 7)
      
      for(drug.iter in drugs) {
        
        for(i in seq(1,
                     length(unique(correction.data$data$bc.data[
                       correction.data$data$bc.data[,drugvar] == drug.iter,    
                       "level"])) , 2 )) {
          
          k1 <- unique(correction.data$data$bc.data[
            correction.data$data$bc.data[,drugvar] == drug.iter, "level"])[i]
          
          n.data <- correction.data$data$bc.data[
            correction.data$data$bc.data$level == k1,]
          
          
          ylim <- range(n.data[, c("absorbance", "NB", "BC2")], na.rm = TRUE)
          n.data <- correction.data$data$bc.data[
            correction.data$data$bc.data$level == k1  ,]
          n.data[, dosevar][n.data[, dosevar] == 0] <-
            min(n.data[, dosevar][n.data[, dosevar]>0])/2
          
          
         
          
          
          if(is.null(xlab.orig))
            xlab <- "Concentration"
          
          if(is.null(main.orig))
            main = paste(n.data[,namevar][1], " (", n.data[,drugvar][1],")", sep = "")
          
          if(is.null(ylab.orig))
            ylab <- c("Raw", "Simple Correction", "Model Correction")
          
          
          if(set.par)
          par(mfrow = c(1,3))
          
          pch <- rep(1, nrow(n.data))
          pch[n.data[, additivevar] %in% c(backgroundval, controlval)] <-
            ifelse(n.data[, additivevar][
              n.data[, additivevar] %in% 
                c(backgroundval, controlval)] == controlval, 2, 3)
          pch[n.data$outlier == 1] <- 4 
          
          if(col.plate)
            col <- col[as.numeric(as.factor(n.data$plate))]
          
          plot(n.data$absorbance ~ n.data[, dosevar], log = "x",
               col = col,
               pch = pch,
               main = main, xlab = xlab,
               ylab = ylab[1],
               ylim = ylim)
          
          legend("bottomleft", pch = 1, col = 1:length(n.data$plate),
                 legend = levels(as.factor(n.data$plate)), bty = "n")
          legend("bottomright", pch = 1:3, legend = c("Drug", controlval, backgroundval), bty = "n")
          
          
          plot(n.data$NB ~ n.data[, dosevar], log = "x",
               col = col,
               pch = pch,
               main = main, xlab = xlab,
               ylab = ylab[2],
               ylim= ylim )
          
          plot(n.data$BC2 ~ n.data[, dosevar], log = "x",
               col = col,
               pch = pch,
               main = main, xlab = xlab,
               ylab = ylab[3],
               ylim= ylim )
          
        }     
      }
      if(pdfit)
        dev.off()
    }
    invisible()
  }

## Function for plotting the result of the fitted growth model


plot.growthModel <- function(x,
                             ...,
                             time.points.used = c("all", "48"),
                             drugs  = NULL,
                             names  = NULL,
                             conc.names = NULL,
                             ylim   = NULL,
                             xlim   = NULL,
                             nrows  = NULL,
                             ncols  = NULL,
                             xlab   = "Time (Hours)",
                             ylab   = "Absorbance",
                             #scale  = TRUE,
                             main   = "Dose Response Growth Curves",
                             absorbance.CI       = FALSE,
                             absorbance.CI.col   = "#C6C6C5",
                             absorbance.CI.alpha = 80,
                             bootstrap.conf = TRUE,
                             barcol         = "grey",
                             bar.height     = 1.3,
                             plotgrid = FALSE,
                             grid.col = "#C6C6C5",
                             grid.lty = 1,
                             grid.lwd = 1,
                             bs.col   = "#C6C6C5",
                             bs.lty   = 1,
                             bs.lwd   = 0.5,
                             bs.alpha = 80, 
                             line.col =  "#662D91",
                             line.lty = 1,
                             line.lwd = 1,
                             line.alpha = "",
                             plot.data  = TRUE,                            
                             col.by.identifier = TRUE,
                             col.points = "#9F9692",         
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
                            figure.output = getwd(),
                             pdf.width  = 6.6929, 
                             pdf.height = 6.6929,
                             pointsize = 8){
  
  old.par <- par(no.readonly = TRUE)
  on.exit(par(old.par))
  
  A.data <- x
  call <- A.data$auxiliary$passed.var
  back <- call$backgroundval 
  additivevar <- call$additivevar
  namevar  <- call$namevar
  idvar    <- call$idvar
  identifier <- call$identifier
  drugvar  <- call$drugvar
  names.orig <- names
  if(is.null(drugs))
    drugs <- unique(A.data$data$bc.data[, drugvar])
  
  dosevar <- call$dosevar
  timevar <- call$timevar
  #time.points <- A.data$GM.time.points
  ylim.orig <- ylim
  xlim.orig <- xlim
  
  is.odd <- function(x) x %% 2 != 0
  is.even <- function(x) x %% 2 == 0
  in.interval <- function(x, interval){ 
    stopifnot(length(interval) == 2L) 
    interval[1] <= x & x <= interval[2] 
  } 
  
  for(drug.iter in drugs){
    #(drug.iter <- drugs[2])
    if("bootstrap" %in% class(A.data)){
      data.drug <- A.data$data$bs.raw.data[
        A.data$data$bs.raw.data[,drugvar] == drug.iter, ]
    }else{
      data.drug <- A.data$data$bc.data[
        A.data$data$bc.data[,drugvar] == drug.iter, ]
    }
    #  c(namevar, drugvar, idvar,"setupdate","outlier", "type", 
    #    "level","plate","Additive","Concentration", 
    #    "time", "NB", "BC", "BC2", "Std.Error")]
    data.drug$up <-  data.drug$BC + 1.96 * data.drug$Std.Error
    data.drug$lo <-  data.drug$BC - 1.96 * data.drug$Std.Error
    
    if(is.null(names.orig)){
      names <- unique(data.drug[, namevar])
    }else{
      names <- names.orig[names.orig %in% data.drug[, namevar]]
    }
    (name.iter <- names[2])
    
    for(name.iter in names){
      if(pdfit){
        pdf(file = file.path(figure.output, 
                             paste(drug.iter, "_",name.iter, 
                                   ".pdf", sep = "" )),useDingbats = FALSE,
            width = pdf.width, height = pdf.height, pointsize = pointsize)
      }else{
        if(!(drug.iter == drugs[1] & name.iter == names[1]))
        plot.new()
      }
      if("all" %in% time.points.used){
        sum <- A.data$fits$growthModel[[drug.iter]][[
          name.iter]][["BC"]][["all"]]$summary
      }else{
        sum <- A.data$fits$growthModel[[drug.iter]][[
          name.iter]][["BC"]][[paste(max(time.points.used))]]$summary
      }
      
      
      if(!is.null(sum)){
        
        data <- data.drug[data.drug[, namevar] == name.iter & 
                            data.drug[, additivevar] != back, ]
        
        if(is.null(conc.names)){
          conc.names2 <- unique(data$type)
        }else{
          conc.names2 <-conc.names
        }
        
        new.types <- data[match(conc.names2, data$type), "new.type"]
        
        t <- seq(min(data[,timevar], na.rm = TRUE), 
                 max(data[,timevar]), length.out = 200)
        sum <- sum[rownames(sum) %in% new.types, ]
        concs <- sum[, dosevar]
        
        data[, timevar] <- data[, timevar] + 1 
        
        if(is.null(ylim.orig)){   
          ylim = c(range(c(data$BC2, data$NB)))
          ylim  = c(range(c(data$BC2)))
        }
        
        if(is.null(xlim.orig))
          xlim <- range(data[,timevar]) + c(-2, 3)
        if(grepl("y", log)){
          ylim[ylim <= 0] <- 0.025 
          #          ylim <- log(ylim)
        }
        #         
        #         if(grepl("x", log))
        #           xlim <- log(xlim)
        
        add <- 0
        if(plot.all)
          add <- 1
        concs <- sum[, dosevar]
        if(is.null(nrows) & is.null(ncols)){
          nrows <- 1
          ncols <- ceiling((length(concs) + add)/nrows)
          while(ncols > nrows){
            nrows <- nrows + 1
            ncols <- ceiling((length(concs) + add)/nrows)
          }
          
          if(ncols*nrows != (length(concs) + add) & 
               (ncols+1)*(nrows -1)==(length(concs) + add)){
            ncols <- ncols + 1
            nrows <- nrows - 1
          }
        }
        
        n.plots <- (length(concs) + add)
        
        if(is.null(ncols))
          ncols <- ceiling((length(concs) + add)/nrows)
        if(is.null(nrows))
          nrows <- ceiling((length(concs) + add)/ncols)
        
        layout <- data.frame(col = rep(1:ncols, nrows),
                             row = rep(1:nrows, each = ncols))
        
        
        if(ncols*nrows < n.plots)
          stop("The plot does not include enough panels")
        par(mfrow = c(nrows, ncols))
        
        par(mar = c(0, 0, bar.height, 0))
        par(oma = c(4, 4, 4.2, 3))
        
        layout(matrix(1:(nrows * ncols), nrow = nrows, 
                      ncol = ncols, byrow=TRUE))
        
        which.packet <- 0
        
        
        if(FALSE) {  
          BS.iters <- c("BC", paste("BS:", 1:50, sep = "")  )
          MSE.list <- list()
          
          for(time.iter in time.points.used){
            mat <- matrix(0, nrow=length(concs), ncol = length(BS.iters))
            colnames(mat) <- BS.iters
            
            for(i in 1:length(concs)){
              for(BS.iter in BS.iters){
                c1 <- data[data[, dosevar]==sum[, dosevar][i] &
                             data$type != "B"  ,]
                
                c2 <- c1[!duplicated(c1[, timevar]), ]
                
                fit <- A.data$fits$growthModel[[drug.iter]][[
                  name.iter]][[BS.iter]][[time.iter]]$fit
                sum.res <- A.data$fits$growthModel[[drug.iter]][[
                  name.iter]][[BS.iter]][[time.iter]]$summary
                sum.res <- sum.res[rownames(sum), ]
                if(i == 1){
                  pred.res <- sum.res$N0[1]*2^(c2[, timevar]/sum.res$T0[i])
                }else{
                  pred.res <- sum.res$N0[1]*2^(c2[, timevar]/sum.res$TC[i])
                }
                
                mat[i, BS.iter]       <-  mean((c2[, BS.iters[1]] - pred.res)^2, na.rm = TRUE)
                
                
              }  
            }
            MSE.list[[time.iter]] <- mat 
          }
          
          plot( ((MSE.list[["48"]] - MSE.list[["all"]] ) /  MSE.list[["48"]])[16,]  ) 
          
          MSE.list[["48"]] / MSE.list[["all"]]  
        }
        for(i in 1:length(concs)){
          which.packet <- which.packet + 1
          
          
          
          
          matplot(x = xlim, y= ylim,
                  type = "n",
                  axes = FALSE,
                  log = log,
                  cex.axis = 0.7)#,
          #...)
          
          if(plotgrid)
            grid(col =  grid.col, lty = grid.lty, lwd = grid.lwd )
          box()
          par(xpd=NA)
          
          rect(grconvertX(0, from='nfc'), grconvertY(1, from='npc'),
               grconvertX(1, from='nfc'), grconvertY(1, from='nfc'),
               col = barcol)
          
          par(xpd = FALSE)
          c <- i-1
          if(i == 1){
            
            title <- as.expression(
              substitute(paste(C[c], " = ", a, ", ", T[0], " = ", b, sep = ""), 
                         list(a = signif(concs[i], 2),
                              b = round(sum$T0[i]),
                              C = substr(conc.names2[i], 1,1),
                              c = substr(conc.names2[i], 
                                         2,nchar(conc.names2[i])))))
            
            title(title, line = 0.6,font = 4)
          }else{
            title <- 
              as.expression(substitute(
                paste(C[c], " = ", a, ", ", T[c], " = ", b, sep = ""), 
                list(a = signif(concs[i], 2),
                     b = round(sum$TC[i]),
                     C = substr(conc.names2[i], 1,1),
                     c = substr(conc.names2[i], 
                                2,nchar(conc.names2[i])))))
            
            title(title, line = 0.6,
                  font = 4)
          }
          # rect(1,2,30,40, col = 1)
          
          
          (spare <- nrow(layout) -  n.plots -1)
          spare <- c(ncols - spare, ncols)
          if(layout[which.packet, 2] == nrows & 
               is.odd(layout[which.packet, 1])) # bottom
            axis(1)
          if(layout[which.packet, 2] == (nrows -1) & 
               is.odd(layout[which.packet, 1]) & 
               in.interval(layout[which.packet, 1], spare) ) # bottom
            axis(1)
          if(layout[which.packet, 1] == 1 & 
               is.odd(layout[which.packet, 2])) # left side
            axis(2, las = 2)
          if(layout[which.packet, 2] == 1 & 
               is.even(layout[which.packet, 1])) # top
            axis(3, outer = TRUE)
          if(layout[which.packet, 1] == ncols & 
               is.even(layout[which.packet, 2])) #right side
            axis(4, las = 2)
          if(layout[which.packet, 2] == nrows & 
               is.even(layout[which.packet, 2]) 
             & which.packet == n.plots) #right side
            axis(4, las = 2)
          
          # axis(4, labels = FALSE)
          
          #####################
          ## The data is made
          #####################
          
          c1 <- data[data[, dosevar]==sum[, dosevar][i] &
                       data$type != "B"  ,]
          
          c2 <- c1[c1$id == c1$id[1]  ,]
          c2 <- c2[!duplicated(c2[, timevar]), ]
          
                 
          
          #####################
          ## Confodence interval on the points
          #####################
          if(absorbance.CI)
            polygon(x = c((c2[, timevar]), rev(c2[, timevar])), 
                    y = c(c2$lo, rev(c2$up)), 
                    col= paste(absorbance.CI.col, 
                               absorbance.CI.alpha, sep = ""), border =0)
          
          
          if(bootstrap.conf & "bootstrap" %in% class(A.data)){
            if(length(bs.col) < length(time.points.used))
              bs.col <- rep(bs.col, length(time.points.used))
            if(length(bs.lty) < length(time.points.used))
              bs.lty <- rep(bs.lty, length(time.points.used))
            if(length(bs.lwd) < length(time.points.used))
              bs.lwd <- rep(bs.lwd, length(time.points.used))
            
            for(bs.iter in 1:A.data$call$bootstrap$n.samples)
              for(j in 1:length(time.points.used)){
                sum.res <- A.data$fits$growthModel[[drug.iter]][[
                  name.iter]][[paste("BS:", bs.iter, sep ="")]][[
                    paste(time.points.used[j])]]$summary
                sum.res <- sum.res[rownames(sum), ]
                if(i == 1){
                  pred.res <- sum.res$N0[1]*2^(t/sum.res$T0[i])
                }else{
                  pred.res <- sum.res$N0[1]*2^(t/sum.res$TC[i])
                }
                lines(t, pred.res, type = "l", 
                      lty = bs.lty[j], 
                      col = paste(bs.col[j], bs.alpha, sep = ""), 
                      lwd = bs.lwd[j])
              }
          }
          
          formula <- as.formula(paste("BC2 ~ ", timevar))
          if(plot.data){
            if(col.by.identifier){
              
              col <- col.points[
                as.numeric(as.factor(c1[, identifier]))]
              points(formula, data = c1, 
                     col = col, 
                     pch = ifelse(c1$outlier, 
                                  pch.points.outlier[1],
                                  pch.points[1]))
            }else{
              points(formula, data = c1, 
                     col = col.points[1], 
                     pch = ifelse(c1$outlier, 
                                  pch.points.outlier[1],
                                  pch.points[1]))
            }
          }
          
          
          if(length(line.col) < length(time.points.used))
            line.col <- rep(line.col, length(time.points.used))
          if(length(line.lty) < length(time.points.used))
            line.lty <- rep(line.lty, length(time.points.used))
          if(length(line.lwd) < length(time.points.used))
            line.lwd <- rep(line.lwd, length(time.points.used))
          
          
          for(j in 1:length(time.points.used)){            
            sum.res <- A.data$fits$growthModel[[drug.iter]][[
              name.iter]][["BC"]][[paste(time.points.used[j])]]$summary
            sum.res <- sum.res[rownames(sum), ]
            if(i == 1){
              pred.res <- sum.res$N0[1]*2^(t/sum.res$T0[i])
              col = paste(line.col.C0[j], line.alpha, sep = "")
            }else{
              pred.res <- sum.res$N0[1]*2^(t/sum.res$TC[i])
              col = paste(line.col[j], line.alpha, sep = "")
            }
            lines(t, pred.res, type = "l", 
                  col = col, 
                  lty = line.lty[j],
                  lwd = line.lwd)
          }     
        }
        
        if(plot.all){
          which.packet <- which.packet + 1
          matplot(x = xlim, y= ylim,
                  type = "n",
                  axes = FALSE,
                  cex.axis = 0.7)
          
          if( plotgrid)
            grid(col =  grid.col, lty = grid.lty, lwd = grid.lwd )
          box()
          par(xpd=NA)
          
          rect(grconvertX(0, from='nfc'), grconvertY(1, from='npc'),
               grconvertX(1, from='nfc'), grconvertY(1, from='nfc'),
               col = barcol)
          
          par(xpd = FALSE)
          title(expression(paste("All Concentrations")), line = 0.7, font = 3)
          # rect(1,2,30,40, col = 1)
          
          (spare <- nrow(layout) -  n.plots -1)
          spare <- c(ncols - spare, ncols)
          if(layout[which.packet, 2] == nrows & 
               is.odd(layout[which.packet, 1])) # bottom
            axis(1)
          if(layout[which.packet, 2] == (nrows -1) & 
               is.odd(layout[which.packet, 1]) & 
               in.interval(layout[which.packet, 1], spare) ) # bottom
            axis(1)
          if(layout[which.packet, 1] == 1 & 
               is.odd(layout[which.packet, 2])) # left side
            axis(2, las = 2)
          if(layout[which.packet, 2] == 1 & 
               is.even(layout[which.packet, 1])) # top
            axis(3, outer = TRUE)
          if(layout[which.packet, 1] == ncols & 
               is.even(layout[which.packet, 2])) #right side
            axis(4, las = 2)
          if(layout[which.packet, 2] == nrows & 
               is.even(layout[which.packet, 2]) 
             & which.packet == n.plots) #right side
            axis(4, las = 2)
          for(i in 1:length(concs)){
            c1 <- data[data$Concentration==sum$Concentration[i] &
                         data$type != "B"  ,]
            
            formula <- as.formula(paste("BC ~", timevar))
            
            if(plot.data.all){
              if(col.by.identifier &FALSE){
                col <- col.points[
                  as.numeric(as.factor(c1[, identifier]))]
                points(formula, data = c1, 
                       col = col, 
                       pch = ifelse(c1$outlier, 
                                    pch.points.outlier[1],
                                    pch.points[1]))
              }else{
                points(formula, data = c1, 
                       col = col.points[1], 
                       pch = ifelse(c1$outlier, 
                                    pch.points.outlier[1],
                                    pch.points[1]))
              }
            }
          }
          
          for(i in 1:length(concs)){
            c1 <- data[data[, dosevar]==sum[, dosevar][i] &
                         data$type != "B"  ,]
            
            formula <- as.formula(paste("BC ~", timevar))
            
            for(j in 1:length(time.points.used)){   
              sum.res <- A.data$fits$growthModel[[
                drug.iter]][[name.iter]][["BC"]][[
                  paste(time.points.used[j])]]$summary
              sum.res <- sum.res[rownames(sum), ]
              if(i == 1){
              }else{
                pred.res <- sum.res$N0[1]*2^(t/sum.res$TC[i])
                col <- paste(line.col.all[j], line.alpha, sep = "")
                lines(t, pred.res, type = "l", 
                      col = col, 
                      lty = line.lty[j],
                      lwd = line.lwd)
              }
              
              
            }
            
            pred.res <- sum.res$N0[1]*2^(t/sum.res$T0[i])
            col <- paste(line.col.C0[1], line.alpha, sep = "")
            
            lines(t, pred.res, type = "l", 
                  col = col, 
                  lty = line.lty[j],
                  lwd = line.lwd)
            
            pred.res <- sum.res$N0[1]*2^(t/(sum.res$T0[1]*2))
            lines(t, pred.res, type = "l", 
                  col = line.col.GI50, 
                  lty = line.lty[j],
                  lwd = line.lwd)
            
            pred.res <- rep(sum.res$N0[1], length(t))
            lines(t, pred.res, type = "l", 
                  col = line.col.TGI, 
                  lty = line.lty[j],
                  lwd = line.lwd)
            
            pred.res <- sum.res$N0[1]*2^(t/(-48))
            lines(t, pred.res, type = "l", 
                  col = line.col.LC48, 
                  lty = line.lty[j],
                  lwd = line.lwd)
            
          }
        }
        
        mtext(paste(main, "for", name.iter), outer = TRUE, line = 2.2, cex = 1.3)
        mtext(ylab, side = 2, outer = TRUE, line = 2.5)
        mtext(xlab, side = 1, outer = TRUE, line = 2.8)
        if(pdfit)
          dev.off()
      }
    }
  }
  invisible()
}



plot.growthModel.resid <- function(x,
                             ...,
                             time.points.used = c("all"),
                             drugs  = NULL,
                             names  = NULL,
                             xlab   = "Time (Hours)",
                             ylab   = "Absorbance",
                             main   = "Dose Response Growth Curves",
                             pdfit      = FALSE,
                             figure.output = getwd(),
                             pdf.width  = 6.6929, 
                             pdf.height = 6.6929,
                             pointsize  = 8){
  
  old.par <- par(no.readonly = TRUE)
  on.exit(par(old.par))
  
  A.data <- x
  call <- A.data$auxiliary$passed.var
  back <- call$backgroundval 
  additivevar <- call$additivevar
  namevar  <- call$namevar
  idvar    <- call$idvar
  identifier <- call$identifier
  drugvar  <- call$drugvar
  names.orig <- names
  if(is.null(drugs))
    drugs <- unique(A.data$data$bc.data[, drugvar])
  
  dosevar <- call$dosevar
  timevar <- call$timevar
  #time.points <- A.data$GM.time.points
  ylim.orig <- ylim
  xlim.orig <- xlim
  
  for(drug.iter in drugs){
    #(drug.iter <- drugs[2])
    if("bootstrap" %in% class(A.data)){
      data.drug <- A.data$data$bs.raw.data[
        A.data$data$bs.raw.data[,drugvar] == drug.iter, ]
    }else{
      data.drug <- A.data$data$bc.data[
        A.data$data$bc.data[,drugvar] == drug.iter, ]
    }
    #  c(namevar, drugvar, idvar,"setupdate","outlier", "type", 
    #    "level","plate","Additive","Concentration", 
    #    "time", "NB", "BC", "BC2", "Std.Error")]
    
    if(is.null(names.orig)){
      names <- unique(data.drug[, namevar])
    }else{
      names <- names.orig[names.orig %in% data.drug[, namevar]]
    }
    (name.iter <- names[2])
    
    for(name.iter in names){
      if(pdfit){
        pdf(file = file.path(figure.output, 
                             paste(drug.iter, "_",name.iter, 
                                   ".pdf", sep = "" )),useDingbats = FALSE,
            width = pdf.width, height = pdf.height, pointsize = pointsize)
      }else{
        if(!(drug.iter == drugs[1] & name.iter == names[1]))
          plot.new()
      }
      if("all" %in% time.points.used){
        sum <- A.data$fits$growthModel[[drug.iter]][[
          name.iter]][["BC"]][["all"]]$summary
      }else{
        sum <- A.data$fits$growthModel[[drug.iter]][[
          name.iter]][["BC"]][[paste(max(time.points.used))]]$summary
      }
      
      
      if(!is.null(sum)){
        
        data <- data.drug[data.drug[, namevar] == name.iter & 
                            data.drug[, additivevar] != back, ]
        
        
        conc.names2 <- unique(data$type)
        
        
        new.types <- data[match(conc.names2, data$type), "new.type"]
        
        t <- seq(min(data[,timevar], na.rm = TRUE), 
                 max(data[,timevar]), length.out = 200)
        sum <- sum[rownames(sum) %in% new.types, ]
        concs <- sum[, dosevar]
        
        # incubation
        data[, timevar] <- data[, timevar] + 1 
        
        
        
        #####################
        ## The data is made
        #####################
        
        #View(data[data$type != "B"  ,])
        
        resid <- vector()
        predicted <- vector()
        measured <- vector()
        
        c.number <- vector()
        
        for(i in 1:length(concs)){
          
          c1 <- data[data[, dosevar]==sum[, dosevar][i] &
                       data$type != "B"  ,]
          
          c2 <- c1[!duplicated(c1[, timevar]), ]
          
          fit <- A.data$fits$growthModel[[drug.iter]][[
            name.iter]][["BC"]][[paste(time.points.used)]]$fit
          sum.res <- A.data$fits$growthModel[[drug.iter]][[
            name.iter]][["BC"]][[paste(time.points.used)]]$summary
          sum.res <- sum.res[rownames(sum), ]
          if(i == 1){
            pred.res <- sum.res$N0[1]*2^(c2[, timevar]/sum.res$T0[i])
          }else{
            pred.res <- sum.res$N0[1]*2^(c2[, timevar]/sum.res$TC[i])
          }
          resid     <- c(resid, c2[,"BC"] - pred.res)
          predicted <- c(predicted, pred.res)
          measured  <- c(measured, c2[,"BC"])
          c.number  <- c(c.number, i)
          
        }
        resid.data <- data.frame(measured  = measured,
                                 predicted = predicted, 
                                 resid     = resid,
                                 std       = fit$sigma,
                                 c.number  = c.number)
        
        
      }
      
      
      
      plot(resid / resid.data$std ~ predicted,
           xlab = xlab, ylab = ylab)
      abline(h = 0)
      
    }
  }
}







plot.DRdata <- 
  function(x, ...,
           drug         = 1,
           ylim         = NULL, 
           xlim         = NULL, 
           main         = drug,
           plot.data    = FALSE,
           model        = "G", 
           type         = "AUC", 
           ylab         = NULL, 
           xlab         = NULL, 
           cex          = 1,
           names        = NULL, 
           times        = NULL, 
           col.scheme   = NULL,
           color.palette = "Dark2" ,
           n.colors     = NULL,
           reverse.col  = FALSE,
           col          = NULL,  
           lty          = NULL,
           pch          = NULL,
           lwd = 2,
           legend.place = "bottomleft",
           plot.order   = NULL,
           dose.scale   = "mol/l",
           dose.logfun  = "log10",
           legend       = TRUE, 
           n.columns     = 2, 
           legend.cex   = 1,
           use.col.order = FALSE,
           split        = "all"
           )
{
    
 #   old.par <- par(no.readonly = TRUE) 
 #   on.exit(par(old.par))
    
    A.data <- x
    ###############################
    ## Function starts
    ###############################
    
    if(is.numeric(drug))
      drug <- names(A.data$data$iso.fits)[drug]
    
    y.min <- NULL
    
    call <- A.data$auxiliary$passed.var
    
    namevar <- eval(call$namevar)
    timevar <- eval(call$timevar)
    drugvar <- eval(call$drugvar)
    dosevar <- eval(call$dosevar)
    
    drugvar <- eval(call$drugvar)
    
    ######################
    ## The isotonic regression results are 
    ## chosen for the specified model
    #############################
    
    iso.fit <- A.data$data$iso.fits[[drug]][[model]]
    
    ##############################
    ## The concentration is scaled 
    ## to fit the scale call
    ##############################
    
    mol <- A.data$auxiliary$mol.data[drug, "mol.mass"]
    
    iso.fit$conc.m <- 
      concConvert(iso.fit[, dosevar], mol.mass = mol, 
                  from = iso.fit$unit, to = dose.scale,
                  logfun.from = A.data$call$isoreg.DRdata$logfun.to,
                  logfun.to = dose.logfun)
    
    
    ##############################
    ## The normalised data are used to 
    ## fit the isotonic regession
    ##############################
    
    BC <- "BC"
    
    if(split == "all")
      iso.fit[, split] <- rep(split, nrow(iso.fit))
    
    split2 <- unique(iso.fit[, split])
    
    for(split.iter in split2){
      print(split.iter)
      
      data.drug <- iso.fit[iso.fit[, split] == split.iter, ]
      
      if(is.null(times))
        times <- sort(unique(data.drug[, timevar])) 
      
      cells <- unique(data.drug[, namevar])
      if(!is.null(names))
        cells <- cells[cells %in% names]
      
      if(model != "G"){
        data.drug <- data.drug[data.drug[, namevar] %in% cells & 
                                 data.drug[, timevar] %in% times, ]
      }else{
        data.drug <- data.drug[data.drug[, namevar] %in% cells, ]
      }
      
      ####################################################
      ##
      ## The yaxis scale is defined
      ##
      ####################################################
      wh.range <- ifelse(plot.data, model, paste(model, "iso", sep = "."))
      
      
      ## For model G
      if(grepl("G", model)){
        
        wh.x <- iso.fit[, "data"] == "BC" & 
          iso.fit[, namevar] %in% cells
        if(model != "G")        
          wh.x <- iso.fit[, "data"] == "BC" &
          iso.fit[,timevar] %in% times & 
          iso.fit[, namevar] %in% cells
        
        twoscale <- min(iso.fit[wh.x, wh.range]) < 0
        if(is.null(ylim)){
          maxGup <- max(iso.fit[wh.x, wh.range])
          if(twoscale){
            ylim <- c(-maxGup, maxGup)
          }else{
            ylim <- range(iso.fit[, wh.range])
          }    
        }else{
          y.min <- min(ylim)
          if(min(ylim) < 0){
            twoscale <- TRUE
            ylim <- c(-max(ylim), max(ylim))
          }
        } 
      }
      
      ## For other models
      if(!grepl("G", model)){
        wh.x <- iso.fit[, "data"] == "BC" &
          iso.fit[,timevar] %in% times & 
          iso.fit[, namevar] %in% cells
        
        twoscale <- FALSE
        if(is.null(ylim)) 
          ylim <- range(iso.fit[wh.x, wh.range])
      }
      
      if(is.null(xlim)) xlim <- range(data.drug[, "conc.m"])
      wh <-  paste(model, "iso", sep = ".")
      #########################################################
      ##
      ## The colours are defined
      ##
      #########################################################
      
      #########################
      ##  The order of the cell lines are determined
      ##
      
      col.sc <- iso.fit[wh.x, c(namevar, timevar)]
      col.sc <- col.sc[!duplicated(paste(col.sc[,timevar], 
                                         col.sc[,namevar])), , drop = FALSE]
      
      col.sc <- col.sc[order(col.sc[, timevar], col.sc[, namevar]), , 
                       drop = FALSE]
      
      if(is.null(plot.order)){
        data <- A.data$summary[[drug]][[model]][[type]]
        
        if(model != "G"){
          data.c <- cbind(t(data[[paste(times[1])]][
            1,, drop = FALSE]), time = times[1])
          
          if(length(times >1))
            for(time.iter in times[-1])
              data.c <- rbind(data.c, cbind(t(data[[paste(time.iter)]][
                1,, drop = FALSE]), time = time.iter))
          data.c <- as.data.frame(data.c)
          data.c$name <- row.names(data.c)
          row.names(data.c) <- 1:nrow(data.c)
          data.c <- data.c[data.c[, "time"] %in% times &
                             data.c[, namevar] %in% cells,]
          data.c <- data.c[order(data.c[, "time"], data.c[, namevar]),]
          
          col.sc[, type] <-  data.c[, "BC"]
          
          
          data.c <- reshape(data.c, idvar=namevar, direction="wide")
          row.names(data.c) <- data.c$name
          data.c <- data.c[,-1, drop =FALSE]
          data.c <- apply(data.c, 1, mean)
          data <- data.c[order(data.c)]
          
          
        }else{
          data  <- data[, order(data[1,,drop=FALSE]), drop=FALSE][1,]
          data <- data[names(data) %in% cells]
          col.sc[, type] <-  data[col.sc[,namevar]]
          
        }
        
        names2 <- names(data)
        
      }else{
        names.HS <- colnames(A.data$summary[[drug]][[model]][[type]])
        plot.order <- intersect(plot.order, names.HS)
        names2 <- plot.order
      }
      
      
      
      
      ############################
      ## The number of colors are determined
      ##
      if(!is.null(col)) {
        
        col.sc <- col.sc[order(col.sc[,timevar], col.sc[,type]),]
        
        if(length(col) != nrow(col.sc))
          col <- rep(col, each = ceiling(nrow(col.sc)/length(col) ) )
        if(is.null(lty))
          lty <- rep(1:ceiling(length(cells)/length(col)),
                     nrow(col.sc)/ceiling(length(cells)/
                                            length(col)))[1:nrow(col.sc)]
        
        if(length(lty) != nrow(col.sc))
          lty <- rep(lty, ceiling(nrow(col.sc)/length(lty) ) )[1:nrow(col.sc)]
        
        if(is.null(pch))
          pch <- lty       
        
        if(length(pch) != nrow(col.sc))
          pch <- rep(pch, ceiling(nrow(col.sc)/length(pch) ) )[1:nrow(col.sc)]
        
        #         col.sc$col <- col
        #         col.sc$lty <- lty
        #         col.sc$pch <- pch    
      }
      
      if(is.null(col)){
        
        if(length(cells) == 1 &  length(times) == 1){
          lty <- 1
          n.colors <- 1
        }
        
        if(length(cells) > 1 &  length(times) > 1){ 
          if(is.null(n.colors))
            n.colors <- length(cells)
          if(is.null(lty))
            lty <- 1:length(times)
        }
        
        if(length(cells) == 1  &  length(times) > 1){
          if(is.null(lty) & is.null(n.colors)){
            n.cells <- length(times)
            if(n.cells <= 4)
              lty <- 1
            if(n.cells > 4 & n.cells <= 8)
              lty <- 1:2
            if(n.cells > 8 & n.cells <= 12)
              lty <- 1:3
            if(n.cells > 12)
              lty <- 1:4
          }
          if(is.null(lty) & !is.null(n.colors))
            lty <- ceiling(n.colors/length(lty))   
          if(is.null(n.colors))
            n.colors <- ceiling(length(names2)/length(lty))          
        }
        
        if(length(cells) > 1  &  length(times) == 1){
          if(is.null(lty) & is.null(n.colors)){
            n.cells <- length(names2)
            if(n.cells <= 4)
              lty <- 1
            if(n.cells > 4 & n.cells <= 8)
              lty <- 1:2
            if(n.cells > 8 & n.cells <= 12)
              lty <- 1:3
            if(n.cells > 12)
              lty <- 1:4
          }
          if(is.null(lty) & !is.null(n.colors))
            lty <- ceiling(n.colors/length(lty))   
          if(is.null(n.colors))
            n.colors <- ceiling(length(names2)/length(lty))    
        }
        
        if(is.null(pch))
          pch <- lty     
        
        col <- brewer.pal(ifelse(n.colors < 3, 3, n.colors), 
                          color.palette)[1:n.colors]
        if(!is.null(col.scheme))
          col <- col.scheme[1:length(col)]
        if(reverse.col)
          col <- rev(col)
        
        
        if(length(cells) > 1 &  length(times) > 1){ 
          
          names(col) <- names2
          names(lty) <- times
          names(pch) <- times
          col.sc$col <- col[col.sc[, namevar]]
          col.sc$lty <- lty[paste(col.sc[, timevar])]
          col.sc$pch <- pch[paste(col.sc[, timevar])]
          
        }
        
        if(!(length(cells) > 1 & length(times) > 1)){
          col.s <- rep(col, each = ceiling(nrow(col.sc) / length(col) ))[1:nrow(col.sc)]
          names(col.s) <- names2
          lty.s <- rep(lty, ceiling(nrow(col.sc) / length(lty) ))[1:nrow(col.sc)]
          names(lty.s) <- names2
          pch.s <- rep(pch, ceiling(nrow(col.sc) / length(pch) ))[1:nrow(col.sc)]
          names(pch.s) <- names2
          col.sc$col <- col.s[col.sc[, namevar]]
          col.sc$lty <- lty.s[col.sc[, namevar]]
          col.sc$pch <- pch.s[col.sc[, namevar]]      
        }
        
        col.sc <- col.sc[order(col.sc[,timevar], col.sc[,type]),]
        col <- col.sc$col
        lty <- col.sc$lty
        pch <- col.sc$pch
      }
      
      
      #########################################################
      ##
      ## The labels are defined
      ##
      #########################################################
      
      
      ##################
      ## the ylab
      
      if(is.null(ylab)){
        if(grepl("G", model))
          ylab <- "Growth inhibition - G-model"
        if(model %in% c("D"))
          ylab <- "Growth inhibition - D-model"
        if(model %in% c("R"))
          ylab <- "Growth inhibition - R-model"
      }
      
      ####################
      ## The xlab
      if(is.null(xlab)){
        if(dose.logfun == "log10")
          xlab <- as.expression(substitute(
            paste(conc, " ", log[10], "(", a, ")", sep = ""), 
            list(a = dose.scale, conc = dosevar)))
        
        if(dose.logfun == "log2")
          xlab <- as.expression(substitute(
            paste(conc, " ",log[2], "(", a, ")", sep = ""), 
            list(a = dose.scale, conc = dosevar)))
        
        if(dose.logfun == "log")
          xlab <- as.expression(substitute(
            paste(conc, " ",log, "(", a, ")", sep = ""), 
            list(a = dose.scale, conc = dosevar)))
        
        if(!(dose.logfun %in% c("log10", "log2", "log")))
          xlab <- as.expression(substitute(
            paste(conc, " ", "(", a, ")", sep = ""), 
            list(a = dose.scale, conc = dosevar)))
      }
      
      #########################################################
      ##
      ## The plot is called
      ##
      #########################################################
      
      plot(0,0,xlim, ylim,col = 0,
              type = "n",
              ylab = ylab,
              xlab = xlab,
              ...,
              main = main,
              las = 1,
              cex.axis = 0.7,
              axes = FALSE)
      
      
      axis(1)
      graphics::box()   
      
      iter <- 0
      
      #######################################
      ##
      ## Calculating the new scale for the lower axis
      ##
      ########################################
      normalized   <- function(x, scale = 1) 
        ((x- min(x)) / (max(x)-min(x))) * scale
      normalizeToX <- function(x, y, scale = 1) 
        ((y- min(x)) / (max(x)-min(x)))*scale
      
      if(twoscale){
        y.norm <- c(0, max(ylim))
        #ylim <<- ylim
        a <- max(ylim)*c(0.25, 0.5, 0.75, 1)
        label.ovr <- paste(round(a))
        y.ovr <- (normalizeToX(y.norm,  a, scale = max(ylim)))#+max(ylim))
        
        
        y.norm <- data.drug[data.drug[, namevar] %in% cells, wh.range]
       
        if(is.null(y.min)) y.min <- min(y.norm)
        
        y.norm <- c(y.min, 0)
        
        a <- y.min*c(0.25, 0.5, 0.75, 1)
        dig <- 0# ifelse(y.min < 1/5, 2, 0)
        label <- paste("-1/", -1*round(1/a, dig), sep = "" )
        y.und <- -1*(normalizeToX(y.norm,  a, scale = -max(ylim))+max(ylim))
        
        #y.norm <<- data.drug[data.drug[, namevar] %in% cells, wh.range]
        
       
        if(max(ylim) < 120)  {
          axis(2, at = c(100, 75, 50, 25, 0, y.und),
               labels = c("100", "75", "50", "25", "0", label),
               las = 2)
        }else{
          axis(2, at = c(y.ovr, 0, y.und),
               labels = c(label.ovr, "0", label),
               las = 2)
        }
        
        
      }else{
        axis(2, las = 2)
      }
      
      ############################
      ## matrix for keeping order in legends
      legend.order <- data.frame(timevar = NA, namevar = NA, 
                                 col=NA, lty = NA, pch = NA)
      
      if(model == "G")
        times <- unique(data.drug[,timevar])[1]
      
      for(time.iter in times[times %in% unique(data.drug[,timevar])]) {
        if("G" != model){
          data.time <- data.drug[data.drug[, timevar] == time.iter,]
        }else{
          data.time <- data.drug
        }
        
        data <- A.data$summary[[drug]][[model]][[type]]
        
        if(model != "G")
          data    <- data[[paste(time.iter)]]
        
        #data    <- data[, order(data[1,])]
        data <- data[, order(data[1,,drop=FALSE]), drop=FALSE]
        
        if(is.null(plot.order)){
          names2 <- colnames(data)
        }else{
          names(plot.order) <- plot.order
          int <- intersect(cells, plot.order)
          plot.order <- plot.order[plot.order %in% int]  
          names2 <- plot.order
        }
        
        names(names2) <- names2
        int <- intersect(cells, names2)
        names2 <- names2[names2 %in% int]
        
        
        for(name.iter in names2[names2 %in% unique(data.time[,namevar])]) {
          
          iter <- iter + 1
          
          x <-  data.time[data.time[, namevar] == name.iter &
                            data.time[, "data"] == "BC" , "conc.m"]
          
          y <-  data.time[data.time[,namevar] == name.iter &
                            data.time[, "data"] == "BC",  wh]
          
          if(grepl("G", model))
            y[y<0] <- -1*(normalizeToX(y.norm, y[y < 0], 
                                       scale = -max(ylim))+max(ylim))
          
          
          lines(x, y, col = col[iter], lty=lty[iter])
          
          if(plot.data){    
            
            y <-  data.time[data.time[,namevar] == name.iter &
                              data.time[, "data"] == "BC"&
                              data.time[,"origin"] != "at",  wh.range]
            
            x <-  data.time[data.time[, namevar] == name.iter &
                              data.time[, "data"] == "BC"&
                              data.time[,"origin"] != "at" , "conc.m"]
            if(grepl("G", model))
              y[y<0] <- -1*(normalizeToX(y.norm, y[y < 0], 
                                         scale = -max(ylim))+max(ylim))
            points(x,  y, col = col[iter], pch=pch[iter])
          }
          
          legend.order <- 
            rbind( legend.order, c(time.iter, name.iter, 
                                   col[iter], lty=lty[iter], pch=lty[iter]))
        } # name.iter
      } # time.iter
    }
    
    ###################################################
    ##
    ## The legend is created
    ##
    ###################################################
    if(legend){
      legend.order <- legend.order[-1, ]
      if(length(times) > 1 & length(names2) == 1){
        if(plot.data)
          legend(legend.place, legend = times, ncol = n.columns,
                 #title.col = c("Sensitive", "Resistant"),
                 col = col, lty = lty, pch = lty, bty = "n", cex = legend.cex)
        
        if(!plot.data)
          legend(legend.place, legend = times, ncol = n.columns,
                 #title.col = c("Sensitive", "Resistant"),
                 col = col, lty = lty, bty = "n", cex = legend.cex)
      }
      
      if(length(times) == 1 & length(names2) > 1){
        if(plot.data)
          legend(legend.place, legend = names2, ncol = n.columns,
                 #title.col = c("Sensitive", "Resistant"),
                 col = col, lty = lty, pch = lty, bty = "n", cex = legend.cex)
        
        if(!plot.data)
          legend(legend.place, legend = names2, ncol = n.columns,
                 #title.col = c("Sensitive", "Resistant"),
                 col = col, lty = lty, bty = "n", cex = legend.cex)
      }  
      if(length(times) > 1 & length(names2) > 1){
        if(plot.data)
          legend(legend.place, legend = paste(legend.order[,2], legend.order[,1]), 
                 ncol = n.columns,
                 #title.col = c("Sensitive", "Resistant"),
                 col = legend.order[,"col"], lty = as.numeric(legend.order[,"lty"]), 
                 pch = as.numeric(legend.order[,"pch"]), bty = "n", cex = legend.cex)
        
        if(!plot.data)
          legend(legend.place, legend = paste(legend.order[,2], legend.order[,1]), 
                 ncol = n.columns,
                 #title.col = c("Sensitive", "Resistant"),
                 col = legend.order[,"col"], 
                 lty = as.numeric(legend.order[,"lty"]), bty = "n", cex = legend.cex)
      }  
    }
    invisible()
  }


DRdataBoxplot <- function(A.data, splitvar =NULL, #split = "all", 
                          #curve = FALSE, 
                          splitvar.col = NULL, col.all = NULL,
                          model = "G", 
                          type = "AUC",
                          #q.curve = NULL, 
                          names = NULL,
                          drug = 1,
                          time = NULL,
                          main = NULL,
                          ylab = NULL,
                          clean.names = NULL,
                          plot.order = NULL,
                          dose.scale = "mol/l",
                          dose.logfun = "log10",
                           ...){
  
  
 # old.par <- par(no.readonly = TRUE)  
 #  on.exit(par(old.par))
  
  if(is.numeric(drug))
    drug <- names(A.data$summary)[drug]
  
  ###################################
  ##
  ## Call is supplied
  ##
  ####################################
  call <- A.data$auxiliary$passed.var
  
  namevar <- eval(call$namevar)
  timevar <- eval(call$timevar)
  drugvar <- eval(call$drugvar)
  dosevar <- eval(call$dosevar)
  
  drugvar <- eval(call$drugvar)
  
  ####################################
  ## The dataset is taken out
  
  
  if(model != "G"){
    
    if(is.null(time))
      times <- names(A.data$summary[[drug]][[model]][[type]])
    data <- A.data$summary[[drug]][[model]][[type]][[paste(time)]]
    
    #      data.c <- cbind(t(data[[paste(times[1])]]), time = times[1])
    #           if(length(times >1))
    #             for(time.iter in times)
    #               data.c <- rbind(data.c, cbind(t(data[[paste(time.iter)]]), time = time.iter))
    #           data.c <- as.data.frame(data.c)
    #           data.c$name <- row.names(data.c)
    #           row.names(data.c) <- 1:nrow(data.c)
    #           data.c <- data.c[data.c[, "time"] %in% times,]
    #           data <- data.c[order(data.c[, "time"], data.c[, namevar]),]
    #     
    #    colnames(data)[colnames(data)=="time"] <- timevar
    
    
  }else{
    data <- A.data$summary[[drug]][[model]][[type]]
  }
  
  if(!is.null(names))
    data <- data[, colnames(data) %in% names]
  
  #data <- data[, order(data[1,])]
  data <- data[, order(data[1,,drop=FALSE]), drop=FALSE]
  if(!is.null(plot.order)){
    plot.order <- plot.order[plot.order %in% colnames(data)]
    data <- data[, plot.order]
  }
  
  if(!grepl("AUC", type))
    data <- concConvert(data, mol.mass = A.data$auxiliary$mol.data[drug, "mol.mass"], 
                        logfun.from=A.data$call$summary.DRdata$dose.logfun   ,         
                        from = A.data$call$summary.DRdata$dose.scale, 
                        to = dose.scale,
                        logfun.to = dose.logfun)
  
  if(is.null(ylab)){
    if(dose.logfun == "log10")
      ylab <- as.expression((substitute(
        paste("Concentration ", log[10], "(", a, ")", sep = ""), 
        list(a = dose.scale))))
    
    if(dose.logfun == "log2")
      ylab <- as.expression(substitute(
        paste("Concentration ", log[2], "(", a, ")", sep = ""), 
        list(a = dose.scale)))
    
    if(dose.logfun == "log")
      ylab <- as.expression(substitute(
        paste("Concentration ", log, "(", a, ")", sep = ""), 
        list(a = dose.scale)))
    
    if(!(dose.logfun %in% c("log10", "log2", "log")))
      ylab <- as.expression(substitute(
        paste("Concentration ", "(", a, ")", sep = ""), 
        list(a = dose.scale)))
    
    if(grepl("AUC", type))
      ylab = "Area under dose response curve"
  }
  
  if(is.null(main)){
    if(grepl("GI", type)){
      mp <- ifelse(grepl("G", model), "G", model)
      GIx <- strsplit(type, split="GI")[[1]][2]
      main <- as.expression(substitute(
        paste(a, " induced ", GI[b]^c, sep = ""), 
        list(a = drug, b = GIx, c = mp)))
    }
    if(grepl("TGI", type)){
      mp <- ifelse(grepl("G", model), "G", model)
      main <- as.expression(substitute(
        paste(a, " induced ", TGI^c, sep = ""), 
        list(a = drug, c = mp)))
    }   
    
    if(grepl("LC", type)){
      mp <- ifelse(grepl("G", model), "G", model)
      GIx <- strsplit(type, split="LC")[[1]][2]
      if(grepl("G", model)) GIx <- A.data$call$summary.DRdata$t
      main <- as.expression(substitute(
        paste(a, " induced ", LC[b]^c, sep = ""), 
        list(a = drug, b = GIx, c = mp)))
    }
    if(grepl("AUC", type)){
      mp <- ifelse(grepl("G", model), "G", model)
      q <- A.data$call$summary.DRdata$AUC.q
      main <- as.expression(substitute(
        paste(a, " induced ", AUC[b]^c, sep = ""), 
        list(a = drug, b = q, c = mp)))
    }
    
  }
  col <- NULL
  disease <- NULL
  if(!is.null(splitvar)){
    data.sp <- A.data$meta.list$metadata.new[, c(namevar, splitvar)]
    data.sp[!duplicated(paste(data.sp[,namevar], data.sp[,splitvar])),]
    disease <- data.sp[, splitvar]
    names(disease) <- data.sp[, namevar]
  }
  if(!is.null(disease)){
    data2 <- matrix(NaN, nrow = ncol(data), ncol = 2)
    data2 <- as.data.frame(data2)
    data2[, 1] <- data[1, ]
    data2[, 2] <- disease[colnames(data)]
    colnames(data2) <- c("summary", "disease")
    data.mean <- aggregate(summary ~ disease, FUN = mean, data = data2)
    rownames(data.mean) <- data.mean$disease
    data2$mean <- data.mean[data2$disease, "summary"]
    row.names(data2) <- colnames(data)
    
    data2 <- data2[order(data2[,3], data2[,2], data2[,1]),]
    data <- data[, rownames(data2)]
    col <- as.numeric((as.factor(data2$disease)))
    if(!is.null(splitvar.col))
      col = splitvar.col[data2$disease]
    names(col) <- NULL
  }
  if(!is.null(col.all))
    col = col.all
  
  p <- data[-1,, drop = FALSE]
  if(!is.null(clean.names))
    colnames(p) <- gsub(clean.names, "", colnames(p))
  #  if(model == "G" & type %in% c("AUC.iso", "AUC")) p <- p*100
  
  boxplot(p, las = 2, ylab = ylab, main = main, col = col, ... )
  invisible()
}

 

plotGrid <- function(A.data = A.data,
                     model  = "G",
                     type   = "AUC",
                     times  = NULL,
                     dose.scale   = "mol/l",
                     dose.logfun  = "log10",
                     drug   = 1,
                     name   = 2,
                     names  = NULL,
                     plot.order = NULL,
                     conc.names = NULL,
                     ylim   = NULL,
                     xlim   = NULL,
                     nrows  = NULL,
                     ncols  = NULL,
                     ylab   = NULL, 
                     xlab   = NULL, 
                     #scale  = TRUE,
                     main   = paste("Dose Response Curves for", drug),
                     absorbance.CI       = FALSE,
                     absorbance.CI.col   = "#C6C6C5",
                     absorbance.CI.alpha = 80,
                     bootstrap.conf = TRUE,
                     barcol         = "#71965A",
                     bar.height     = 1.3,
                     plotgrid = TRUE,
                     grid.col = "#C6C6C5",
                     grid.lty = 1,
                     grid.lwd = 1,
                     bs.col   = c("#333333", "#9D2441"),
                     bs.lty   = 1,
                     bs.lwd   = 0.5,
                     bs.alpha = 50, 
                     line.col = c("#333333", "#9D2441"),
                     line.lty = rep(1, 8),
                     line.lwd = rep(1,8),
                     line.alpha = "",                         
                     col.by.identifier = TRUE,
                     col.points = c("#333333", "#9D2441"), 
                     pch = 1,
                     plot.data = FALSE,
                     log        = "",
                     pdfit      = FALSE,
                     pdf.width  = 6.6929, 
                     pdf.height = 6.6929,
                     pointsize  = 8){
  
  old.par <- par(no.readonly = TRUE)
  on.exit(par(old.par))
  
  
  if(is.numeric(drug))
    drug <- names(A.data$data$iso.fits)[drug]
  
  y.min <- NULL
  
  call <- A.data$auxiliary$passed.var
  
  namevar <- eval(call$namevar)
  namevar.orig <- namevar
  if(name == 2 & paste(namevar, 2, sep = "") %in% 
       colnames(A.data$data$raw.data))
    namevar <- paste(namevar, 2, sep = "")
  timevar <- eval(call$timevar)
  drugvar <- eval(call$drugvar)
  dosevar <- eval(call$dosevar)
  
  drugvar <- eval(call$drugvar)
  
  ######################
  ## The isotonic regression results are 
  ## chosen for the specified model
  #############################
  
  iso.fit <- A.data$data$iso.fits[[drug]][[model]]
  
  ##############################
  ## The concentration is scaled 
  ## to fit the scale call
  ##############################
  
  mol <- A.data$auxiliary$mol.data[drug, "mol.mass"]
  
  iso.fit$conc.m <- 
    DoseR:::concConvert(iso.fit[, dosevar], mol.mass = mol, 
                from = iso.fit$unit, to = dose.scale,
                logfun.from = A.data$call$isoreg.DRdata$logfun.to,
                logfun.to = dose.logfun)
  
  ##############################
  
  ##############################
  ## The normalised data are used to 
  ## fit the isotonic regession
  ##############################
  
  BC <- "BC"
  
  
  data.drug <- iso.fit
  
  if(is.null(times))
    times <- sort(unique(data.drug[, timevar])) 
  
  cells <- unique(data.drug[, namevar])
  if(!is.null(names))
    cells <- cells[cells %in% names]
  
  if(model != "G"){
    data.drug <- data.drug[data.drug[, namevar] %in% cells & 
                             data.drug[, timevar] %in% times, ]
  }else{
    data.drug <- data.drug[data.drug[, namevar] %in% cells, ]
  }
  
  ####################################################
  ##
  ## The yaxis scale is defined
  ##
  ####################################################
  wh.range <- ifelse(plot.data, model, paste(model, "iso", sep = "."))
  
  ## For model G
  if(grepl("G", model)){
    
    wh.x <- iso.fit[, "data"] == "BC" & 
      iso.fit[, namevar] %in% cells
    if(model != "G")        
      wh.x <- iso.fit[, "data"] == "BC" &
      iso.fit[,timevar] %in% times & 
      iso.fit[, namevar] %in% cells
    
    twoscale <- min(iso.fit[wh.x, wh.range]) < 0
    if(is.null(ylim)){
      maxGup <- max(iso.fit[wh.x, wh.range])
      if(twoscale){
        ylim <- c(-maxGup, maxGup)
      }else{
        ylim <- range(iso.fit[, wh.range])
      }    
    }else{
      y.min <- min(ylim)
      if(min(ylim) < 0){
        twoscale <- TRUE
        ylim <- c(-max(ylim), max(ylim))
      }
    } 
  }
  
  ## For other models
  if(!grepl("G", model)){
    wh.x <- iso.fit[, "data"] == "BC" &
      iso.fit[,timevar] %in% times & 
      iso.fit[, namevar] %in% cells
    
    twoscale <- FALSE
    if(is.null(ylim)) 
      ylim <- range(iso.fit[wh.x, wh.range])
  }
  
  if(is.null(xlim)) xlim <- range(data.drug[, "conc.m"])
  wh <-  paste(model, "iso", sep = ".")
  
  #####################################
  
  
  ylim.orig <- ylim
  xlim.orig <- xlim
  
  
  #########################
  ##  The order of the cell lines are determined
  ##
  
  col.sc <- iso.fit[wh.x, c(namevar, timevar)]
  col.sc <- col.sc[!duplicated(paste(col.sc[,timevar], 
                                     col.sc[,namevar])), , drop = FALSE]
  
  col.sc <- col.sc[order(col.sc[, timevar], col.sc[, namevar]), , 
                   drop = FALSE]
  
  if(is.null(plot.order)){
    data <- A.data$summary[[drug]][[model]][[type]]
    
    if(model != "G"){
      data.c <- cbind(t(data[[paste(times[1])]][
        1,, drop = FALSE]), time = times[1])
      
      if(length(times >1))
        for(time.iter in times[-1])
          data.c <- rbind(data.c, cbind(t(data[[paste(time.iter)]][
            1,, drop = FALSE]), time = time.iter))
      data.c <- as.data.frame(data.c)
      data.c$name <- row.names(data.c)
      row.names(data.c) <- 1:nrow(data.c)
      data.c <- data.c[data.c[, "time"] %in% times &
                         data.c[, namevar] %in% cells,]
      data.c <- data.c[order(data.c[, "time"], data.c[, namevar]),]
      
      col.sc[, type] <-  data.c[, "BC"]
      
      
      data.c <- reshape(data.c, idvar=namevar, direction="wide")
      row.names(data.c) <- data.c$name
      data.c <- data.c[,-1, drop =FALSE]
      data.c <- apply(data.c, 1, mean)
      data <- data.c[order(data.c)]
      
      
    }else{
      #data   <- data[, order(data[1,])][1,]
      data <- data[, order(data[1,,drop=FALSE]), drop=FALSE][1,]
      data <- data[names(data) %in% cells]
      col.sc[, type] <-  data[col.sc[,namevar]]     
    }
    
    names2 <- names(data)
    
  }else{
    names.HS <- colnames(A.data$summary[[drug]][[model]][[type]])
    plot.order <- intersect(plot.order, names.HS)
    names2 <- plot.order
  }
  
  #######################
  
  is.odd <- function(x) x %% 2 != 0
  is.even <- function(x) x %% 2 == 0
  in.interval <- function(x, interval){ 
    stopifnot(length(interval) == 2L) 
    interval[1] <= x & x <= interval[2] 
  } 
  
  
  if(is.null(nrows) & is.null(ncols)){
    nrows <- 1
    ncols <- ceiling((length(names2))/nrows)
    while(ncols > nrows){
      nrows <- nrows + 1
      ncols <- ceiling((length(names2))/nrows)
    }
    
    if(ncols*nrows != (length(names2)) & 
         (ncols+1)*(nrows -1)==(length(names2))){
      ncols <- ncols + 1
      nrows <- nrows - 1
    }
  }
  
  n.plots <- (length(names2))
  
  if(is.null(ncols))
    ncols <- ceiling((length(names2))/nrows)
  if(is.null(nrows))
    nrows <- ceiling((length(names2))/ncols)
  
  layout <- data.frame(col = rep(1:ncols, nrows),
                       row = rep(1:nrows, each = ncols))
  
 
  par(mfrow = c(nrows, ncols))
  
  par(mar = c(0, 0, bar.height, 0))
  par(oma = c(5.2, 5, 4.2, 3.5))
  
  layout(matrix(1:(nrows * ncols), nrow = nrows, 
                ncol = ncols, byrow=TRUE))
  
  #########################################################
  ##
  ## The labels are defined
  ##
  #########################################################
  
  
  ##################
  ## the ylab
  
  if(is.null(ylab)){
    if(grepl("G", model))
      ylab <- "Growth inhibition - G-model"
    if(model %in% c("D"))
      ylab <- "Growth inhibition - D-model"
    if(model %in% c("R"))
      ylab <- "Growth inhibition - R-model"
  }
  
  ####################
  ## The xlab
  if(is.null(xlab)){
    if(dose.logfun == "log10")
      xlab <- as.expression(substitute(
        paste(conc, " ", log[10], "(", a, ")", sep = ""), 
        list(a = dose.scale, conc = dosevar)))
    
    if(dose.logfun == "log2")
      xlab <- as.expression(substitute(
        paste(conc, " ",log[2], "(", a, ")", sep = ""), 
        list(a = dose.scale, conc = dosevar)))
    
    if(dose.logfun == "log")
      xlab <- as.expression(substitute(
        paste(conc, " ",log, "(", a, ")", sep = ""), 
        list(a = dose.scale, conc = dosevar)))
    
    if(!(dose.logfun %in% c("log10", "log2", "log")))
      xlab <- as.expression(substitute(
        paste(conc, " ", "(", a, ")", sep = ""), 
        list(a = dose.scale, conc = dosevar)))
  }
  
  #############################################
  ##
  ## Scale Functions
  ##
  #############################################
  normalized   <- function(x, scale = 1) 
    ((x- min(x)) / (max(x)-min(x))) * scale
  normalizeToX <- function(x, y, scale = 1) 
    ((y- min(x)) / (max(x)-min(x)))*scale
  
  if(twoscale){
    y.norm <- c(0, max(ylim))
    #ylim <<- ylim
    a <- max(ylim)*c(0.25, 0.5, 0.75, 1)
    label.ovr <- paste(round(a))
    y.ovr <- (normalizeToX(y.norm,  a, scale = max(ylim)))#+max(ylim))
    
    y.norm <- data.drug[data.drug[, namevar] %in% cells, wh.range]

    if(is.null(y.min)) y.min <- min(y.norm)
    
    y.norm <- c( y.min, 0)
    
    a <- y.min*c(0.25, 0.5, 0.75, 1)
    dig <- 0# ifelse(y.min < 1/5, 2, 0)
    label <- paste("-1/", -1*round(1/a, dig), sep = "" )
    y.und <- -1*(normalizeToX(y.norm,  a, scale = -max(ylim))+max(ylim))
    
    if(max(ylim) < 120)  {
      at.axis = c(100, 75, 50, 25, 0, y.und)
      labels.axis = c("100", "75", "50", "25", "0", label)
    }else{
      at.axis = c(y.ovr, 0, y.und)
      labels.axis = c(label.ovr, "0", label)
    }
    
  }
  
  #############################################
  ##
  ## ordering of the data
  ##
  #############################################      
  
  data <- A.data$summary[[drug]][[model]][[type]]
  
  if(model != "G")
    data    <- data[[paste(time.iter)]]
  
  #data    <- data[, order(data[1,])]
  data <- data[, order(data[1,,drop=FALSE]), drop=FALSE]
  
  if(is.null(plot.order)){
    names2 <- colnames(data)
  }else{
    names(plot.order) <- plot.order
    int <- intersect(cells, plot.order)
    plot.order <- plot.order[plot.order %in% int]  
    names2 <- plot.order
  }
  
  names(names2) <- names2
  int <- intersect(cells, names2)
  names2 <- names2[names2 %in% int]
  
  which.packet <- 0
  
  
  
  for(name.iter in names2){
    iter <- 0
    which.packet <- which.packet + 1
    
    data.name <- data.drug[data.drug[, namevar] == name.iter,]
    
    matplot(x = xlim, y= ylim,
            type = "n",
            axes = FALSE,
            log = log,
            cex.axis = 0.7)#,
    #...)
    
    if(plotgrid)
      grid(col =  grid.col, lty = grid.lty, lwd = grid.lwd )
    graphics::box()
    par(xpd=NA)
    
    rect(grconvertX(0, from='nfc'), grconvertY(1, from='npc'),
         grconvertX(1, from='nfc'), grconvertY(1, from='nfc'),
         col = barcol)
    
    par(xpd = FALSE)
    
    
    title(name.iter, line = 0.4, font = 4)
    
    # rect(1,2,30,40, col = 1)
    
    
    (spare <- nrow(layout) -  n.plots -1)
    spare <- c(ncols - spare, ncols)
    if(layout[which.packet, 2] == nrows & 
         is.odd(layout[which.packet, 1])) # bottom
      axis(1)
    if(layout[which.packet, 2] == (nrows -1) & 
         is.odd(layout[which.packet, 1]) & 
         in.interval(layout[which.packet, 1], spare) ) # bottom
      axis(1)
    if(layout[which.packet, 1] == 1 & 
         is.odd(layout[which.packet, 2])) {# left side
      
      if(twoscale){
        at.axis
        axis(2, at = at.axis,
             labels = labels.axis,
             las = 2)
       # axis(2, at = c(100, 75, 50, 25, 0, y.und),
       #       labels = c("100", "75", "50", "25", "0", label),
       #       las = 2)
      }else{
        axis(2, las = 2)
      }  
    }
    if(layout[which.packet, 2] == 1 & 
         is.even(layout[which.packet, 1])) # top
      axis(3, outer = TRUE)
    if(layout[which.packet, 1] == ncols & 
         is.even(layout[which.packet, 2])){# left side
      
      if(twoscale){
        axis(4, at = at.axis,
             labels = labels.axis,
             las = 2)
       # axis(4, at = c(100, 75, 50, 25, 0, y.und),
       #       labels = c("100", "75", "50", "25", "0", label),
       #       las = 2)
      }else{
        axis(4, las = 2)
      }  
    }
    if(layout[which.packet, 2] == nrows & 
         is.even(layout[which.packet, 2]) 
       & which.packet == n.plots){# left side
      
      if(twoscale){
        axis(4, at = c(100, 75, 50, 25, 0, y.und),
             labels = c("100", "75", "50", "25", "0", label),
             las = 2)
      }else{
        axis(4, las = 2)
      }  
    }
    
    if(model == "G")
      times <- unique(data.name[, timevar])[1]
    
    
    if(bootstrap.conf & "bootstrap" %in% class(A.data)){
      
      for(time.iter in times[times %in% unique(data.name[, timevar])]) {
        if("G" != model){
          data.time <- data.name[data.name[, timevar] == time.iter,]
        }else{
          data.time <- data.name
        }
        
        names3 <- unique(data.time[, namevar.orig])
        
        if(length(bs.col) < length(times)*length(names3))
          bs.col <- rep(bs.col, length(times)*length(names3))
        if(length(bs.lty) < length(times)*length(names3))
          bs.lty <- rep(bs.lty, length(times)*length(names3))
        if(length(bs.lwd) < length(times)*length(names3))
          bs.lwd <- rep(bs.lwd, length(times)*length(names3))
        
        
        iter <- 0
        for(name.iter2 in names3){
          iter <- iter + 1
          for(bs.iter in 1:A.data$call$bootstrap$n.samples){
            bs.iter <- paste("BS:", bs.iter, sep ="")
            x <-  data.time[data.time[,namevar.orig] == name.iter2 &  
                              data.time[, "data"] == bs.iter, "conc.m"]
            
            y <-  data.time[data.time[,namevar.orig] == name.iter2 & 
                              data.time[, "data"] == bs.iter,  wh]
            
            if(grepl("G", model))
              y[y<0] <- -1*(normalizeToX(y.norm, y[y < 0], 
                                         scale = -max(ylim))+max(ylim))
            
            
            lines(x, y, type = "l", 
                  lty = bs.lty[iter], 
                  col = paste(bs.col[iter], bs.alpha, sep = ""), 
                  lwd = bs.lwd[iter])
            
          }
        }
      }
    }
    
    iter <- 0
    
    for(time.iter in times[times %in% unique(data.name[, timevar])]) {
      if("G" != model){
        data.time <- data.name[data.name[, timevar] == time.iter,]
      }else{
        data.time <- data.name
      }
      
      names3 <- unique(data.time[, namevar.orig])
      
      
      
      
      
      for(name.iter2 in names3){
        
        iter <- iter + 1
        
        x <-  data.time[data.time[,namevar.orig] == name.iter2 &  
                          data.time[, "data"] == "BC", "conc.m"]
        
        y <-  data.time[data.time[,namevar.orig] == name.iter2 & 
                          data.time[, "data"] == "BC",  wh]
        
        if(grepl("G", model))
          y[y<0] <- -1*(normalizeToX(y.norm, y[y < 0], 
                                     scale = -max(ylim))+max(ylim))
        
        
        lines(x, y, col = line.col[iter], lty=line.lty[iter])
        
        if(plot.data){    
          
          y <-  data.time[data.time[,namevar.orig] == name.iter2 &
                            data.time[, "data"] == "BC"&
                            data.time[,"origin"] != "at",  wh.range]
          
          x <-  data.time[data.time[, namevar.orig] == name.iter2 &
                            data.time[, "data"] == "BC"&
                            data.time[,"origin"] != "at" , "conc.m"]
          if(grepl("G", model))
            y[y<0] <- -1*(normalizeToX(y.norm, y[y < 0], 
                                       scale = -max(ylim))+max(ylim))
          points(x,  y, col = col.points[iter], pch=pch[iter])
        }
      }
    } # name.iter
  } #
  mtext(main, outer = TRUE, line = 2.2, cex = 1.3)
  mtext(ylab, side = 2, outer = TRUE, line = 3.1)
  mtext(xlab, side = 1, outer = TRUE, line = 2.8)
  invisible()
}
