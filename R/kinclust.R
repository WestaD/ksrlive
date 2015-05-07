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

#' Find the lowest pvalue in clustering.
#'
#' \code{kinclust} returns highest au that still has lowest pvalue
#'
#' Function most.stable takes an object created by pvclust
#'  and finds the highest au that has lowest pvalue
#'
#' @param result Object created by pvclust
#' @return result 
#' @return most_stable
#' @return clustering
#' @return data
#'
#' @examples
#' data(mtcars)
#' data <- t(mtcars)
#' result <- pvclust(data)
#' mostab <- most.stable(result)
#'
kinclust <- function(data_list, cl=NULL, ...){
  # insert check for complete data
  # insert check for more than 2 substrates
  clust_class_all <- list()
  result_all <- list()
  mostab_all <- list()
  data <- list()
  #set_seed(seed)
  for (i in 1:length(data_list)) {
    print(i)
    data[[i]] <- data_list[[i]]
    data_t <- t(data[[i]])
    if (!is.null(cl)) {
      result_all[[i]] <- pvclust::parPvclust(cl = cl, data = data_t, ...)
    }else{
      result_all[[i]] <- pvclust::pvclust(data = data_t, ...)
    }
    mostab_all[[i]] <- ksrlive::most.stable(result_all[[i]])
    pvclust_pick <- pvclust::pvpick(result_all[[i]], alpha = mostab_all[[i]],
                               pv = "au", type = "geq", max.only = TRUE)
    clust_class_all[[i]] <- ksrlive::pvclust.clust(pvclust_pick, data_t[[i]])
  }
  names(clust_class_all) <- names(data_list)

  outlist <- list(result = result_all, most_stable = mostab_all,
                clustering = clust_class_all, data = data)
  return(outlist)
}
