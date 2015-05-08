#' Create cluster assignment vector from pvpick object.
#'
#' \code{pvclust.clust} returns a named vector with cluster assignments
#'
#' Function pvclust.clust takes an object created by pvpick from the pvclust 
#' package and transforms it into a named vector with cluster assignments.
#'
#' @param pvclust_pick Object created by pvclust::pvpick
#' @param data Data used for clustering
#' @return named vector with cluster assignments for all instances that were 
#'clustered
#'
#' @examples
#' data(mtcars)
#' data <- t(mtcars)
#' result <- pvclust(data)
#' pvclust_pick <- pvpick(result, alpha = 0.95)
#' clust_class <- pvclust.clust(pvclust_pick, data)
#'
pvclust.clust <- function(pvclust_pick, data){
  clust_assg <- numeric(length = ncol(data))
  names(clust_assg) <- colnames(data)
  for (i in 1:length(pvclust_pick$clusters)) {
    clust_names <- pvclust_pick$clusters[[i]]
    clust_assg[match(clust_names, names(clust_assg))] <- i
  }
  return(clust_assg)
}

#' Find the lowest pvalue in clustering.
#'
#' \code{most.stable} returns highest au that still has lowest pvalue
#'
#' Function most.stable takes an object created by pvclust
#'  and finds the highest au that has lowest pvalue
#'
#' @param result Object created by pvclust
#' @return numeric of the highest au
#'
#' @examples
#' data(mtcars)
#' data <- t(mtcars)
#' result <- pvclust(data)
#' mostab <- most.stable(result)
#'
most.stable <- function(result){
  mostab <- which.max(round(result$edges[1:(nrow(result$edges) - 1), 1], 2))
  mult_mostab <- which(round(result$edges[1:(nrow(result$edges) - 1), 1], 2)
                       >= round(result$edges[mostab, 1], 2))
  ch_mult <- min(result$edges[mult_mostab, 1])
  return(ch_mult)
}

#' Create cluster assignment vector from pvpick object.
#'
#' \code{pvclust.clust} returns a named vector with cluster assignments
#'
#' Function pvclust.clust takes an object created by pvpick from the pvclust 
#' package and transforms it into a named vector with cluster assignments.
#'
#' @param pvclust_pick Object created by pvclust::pvpick
#' @param data Data used for clustering
#' @return named vector with cluster assignments for all instances that were 
#'clustered
#'
#' @examples
#' data(mtcars)
#' data <- t(mtcars)
#' result <- pvclust(data)
#' pvclust_pick <- pvpick(result, alpha = 0.95)
#' clust_class <- pvclust.clust(pvclust_pick, data)
#'
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
