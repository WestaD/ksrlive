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
#' library(pvclust)
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
#' library(pvclust)
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

#' site specific kinase substrate relationship dataset
#'
#' This dataset contains all site specific kinase relationships from 
#' PhosphoSitePlus, PhosphoElm, HPRD and PhosphoPoint.
#'
"phosphonetwork.df"

#' Time course data for phosphorylation sites
#'
#' This dataset contains time course data of phosphorylation sites after insulin
#' stimulation.
#'
#' @source Humphrey et al., Cell Metabolism, 2013
"data_kin"

#' List containing time course data for phosphorylation sites
#'
#' This dataset contains time course data of phosphorylation sites after insulin
#' stimulation.
#'
#' @source Humphrey et al., Cell Metabolism, 2013
"data_list"