shadedsigma <- function(clust_class, data){
  plotdata <- data[names(clust_class)[which(clust_class == 1)], ]
  comp <- which(stats::complete.cases(plotdata))
  meandata <- apply(plotdata[comp, ], 2, mean)
  sddata <- apply(plotdata[comp, ], 2, stats::sd)
  sdupdata <- meandata + sddata
  sddowndata <- meandata - sddata
  statdata <- data.frame(sddowndata, meandata, sdupdata)
  return(statdata)
}

startzero <- function(shadedata, withs = TRUE){
  datazero <- shadedata
  if (withs) {
    for (i in 1:3) {
      datazero[, i] <- datazero[, i] - rep(shadedata[1, 2], nrow(shadedata))
    }
  } else {
    datazero <- datazero - rep(shadedata[1], nrow(shadedata))
  }
  return(datazero)
}

plotshaded <- function(shadedata, color, x = c(1:nrow(shadedata))){
  scales::polygon(c(x, rev(x)), c(shadedata[, 1], rev(shadedata[, 2])),
                  col = alpha(color, 0.3), border = NA)
  scales::polygon(c(x, rev(x)), c(shadedata[, 2], rev(shadedata[, 3])),
                  col = alpha(color, 0.3), border = NA)
  for (i in 1:3) {
    graphics::lines(x, shadedata[, i], type = "l", col = color)
  }
}

mysubstr <- function(char, start, end, sub = "_"){
 subresult <- substr(char, start, end)
 resultst <- character()
  if (!is.na(start)) {
    if (start < 1) {
      addsub <- 1 - start
      start <- 1
      resultst <- paste(c(rep(sub, addsub), subresult), collapse = "")
    } else if (end > nchar(as.character(char))){
      addsub <- end - nchar(as.character(char))
      end <- nchar(as.character(char))
      resultst <- paste(c(subresult, rep(sub, addsub)), collapse = "")
    }else{
      resultst <- subresult
    }
  }
  return(resultst)
}

plotlines <- function(plotnames, data, clust.class = c(1:length(plotnames)), 
                      colors = rainbow(length(plotnames)), 
                      scaled = TRUE, label = "none"){
  if (scaled == TRUE) {
    yrange <- c(0, 1)
    ind <- grep("scaled", colnames(data))
    plotdata <- as.matrix(data[plotnames, ind])
    ylabel <- "scaled logfc"
  }else{
    ind <- grep("log2..Median", colnames(data))
    ind <- grep("mean", colnames(data))
    plotdata <- cbind(rep(0, length(plotnames)),
                      as.matrix(data[plotnames, ind]))
    colnames(plotdata)[1] <- "0_mean"
    yrange <- range(plotdata, na.rm = TRUE)
    ylabel <- "logfc"
  }
  xrange <- range(1, ncol(plotdata))
  comp <- which(stats::complete.cases(plotdata))
  graphics::plot(xrange, yrange, type = "n", xaxt = "n",
                 xlab = "Time", ylab = ylabel)
  for (i in 1:length(comp)) {
    graphics::lines(c(1:ncol(plotdata)), plotdata[comp[i], ], type = "l",
                    col = colors[clust.class[comp[i]]])
    if (label == "all") {
      mtext(side = 4, at = plotdata[comp[i], ncol(plotdata)],
            text = data[rownames(plotdata)[comp[i]], "Gene.Name"],
            col = colors[clust.class[comp[i]]], line = 0.3, las = 2)
    }
  }
  # axis(1,at=c(1,2,3,4,5,6,7,8,9),
  # labels=c('0','15s','30s','1min','2min','5min','10min','20min','60min'))
  axis(1, at = c(1:ncol(plotdata)),
       labels = sapply(strsplit(colnames(plotdata), split = "_"), "[[", 1))
} 
