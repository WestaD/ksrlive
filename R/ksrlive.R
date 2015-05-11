#' Identify site specific kinase substrate relationships using dynamic data.
#'
#' @description Using this package you can combine known site specific 
#' kinase substrate relationships with experimental data and determine active 
#' kinases and their substrates.
#' 
#' @author Westa Domanova
#' @docType package
#' @name ksrlive
NULL

#' Create a kinase substrate list from a data frame
#'
#' \code{KSR.list} returns list of kinase substrate relationships
#'
#' The function KSR.list creates a list of kinase substrate relationships from
#' a data frame and can combine kinase families into one list.
#'
#' @param db data frame of kinase substrate relationships with substrate 
#' identifier in the first column and kinase identifier in the second column.
#' Data base substrate/kinase identifiers should correspond to data 
#' substrate/kinase identifiers.
#' @param kinasefamilies list of kinase identifiers that have to be combined, 
#' one list per kinase family
#' @param exclusive logical, if TRUE only substrates exclusive to the kinase
#' will be included in the list 
#' @return list with kinase substrate relationships, with the kinase
#' identifiers as the list names
#'
#' @examples
#' data(phosphonetwork)
#' test_db <- do.call(paste, 
#' c(phosphonetwork.df[ ,c("SUB_ACC_ID.human", "MODSITE_SEQ.human")],
#' sep = "_"))
#' test_db <- data.frame(substrate = test_db, 
#' kinase = phosphonetwork.df[ ,"KIN_ACC_ID.human"], 
#' stringsAsFactors = FALSE)
#' kin_list <- KSR.list(test_db)
#' # using family lists
#' fam <- list(akt = c("P31749", "P31751"), 
#' rsk = c("Q15418", "P23443"))
#' kin.fam <- KSR.list(test_db, kinasefamilies = fam)
#' kin.fam.exc <- KSR.list(test_db, kinasefamilies = fam, exclusive = TRUE)

KSR.list <- function(db, kinasefamilies = NULL, exclusive = FALSE){
  temp <- split(db[,1], f = db[,2])
  # out<-temp
  out_cl <- lapply(temp, unique) ## delete duplicates
  out_cl <- lapply(out_cl, function(x){as.character(na.omit(x))}) ## delete NAs
  # delete empty lists
  full <- which(sapply(out_cl, function(x){length(x) > 0}))
  if (length(full) > 0) {
    out_cl <- out_cl[full]
    # combine kinasefamilies together
    if (is.null(kinasefamilies)) {
      out_fam <- out_cl
    }else{
      out_fam <- lapply(kinasefamilies, 
                        function(x){unique(unlist(out_cl[unlist(x)]))})
      names(out_fam) <- sapply(kinasefamilies, "[[", 1)
      
      out_cl <- out_cl[-which(names(out_cl) %in% unlist(kinasefamilies))]
      out_fam <- append(out_cl, out_fam)
    }
    if (!exclusive) {
      out_final <- out_fam
    }else{
      fam_df <- data.frame(sub = unlist(out_fam), 
                           kin = names(out_fam)[rep(seq_along(out_fam), 
                                                    lapply(out_fam, length))],
                           stringsAsFactors = FALSE)
      substr <- split(fam_df[ , 2], f = fam_df[ , 1])
      substr_cl <- lapply(substr, unique) ## delete duplicates
      ### find exclusive substrates
      sub_kinases <- sapply(substr_cl, length)
      ## only one kinase
      onekin <- which(sub_kinases == 1)
      # twokin<-which(sub.kinases==2)
      fam_df<-fam_df[fam_df[ , 1] %in% names(substr_cl)[onekin], ]
      temp <- split(fam_df[ , 1], f = fam_df[ , 2])
      out_final <- lapply(temp, unique) ## delete duplicates
    }
  }else{
    print("No lists available")
    out_final <- NULL
  }
  return(out_final)
}
#' Integrate exclusive and complete clustering results 
#'
#' \code{clust.expand} returns a list of clustering assignments
#'
#' The function clust.expand takes two objects created by kinclust, one 
#' clustered using exclusive substrates and the other one all substrates and 
#' expands the core found by clustering exclusive substrates using a pvalue
#' threshold. The core sites 
#' can be tested for differential regulation if a list of differentially 
#' regulated sites is included (recommended). 
#'
#' @param kin_clust list of kinase substrate relationships with only exclusive
#' substrates
#' @param kin_clust_all list of all available kinase substrate relationships
#' @param thre threshold of au value for clustering, for example au = 0.95 
#' corresponds to a pvalue of 0.05
#' @param index numeric vector of indices for which to do the integration
#' @param diff_list character vector of names of differentially regulated
#' substrates
#' @return expand_clust_list integrated clustering assignments
#' @return noclust numeric vector of indices where a cluster could not be 
#' determined
#'
#' @examples
#' data(phosphonetwork)
#' # create identifier for substrate
#' test_db <- do.call(paste, 
#' c(phosphonetwork.df[ ,c("SUB_ACC_ID.human", "MODSITE_SEQ.human")],
#' sep = "_"))
#' test_db <- data.frame(substrate = test_db, 
#' kinase = phosphonetwork.df[ ,"KIN_ACC_ID.human"], 
#' stringsAsFactors = FALSE)
#' 
#' data(datalist)
#' 
#' # create identifier in data
#' nam_map <- do.call(paste, 
#' c(data_kin[, c("Uniprot.human", "Motif.human")], sep = "_"))
#' ind_map <- match(test_db[ ,"substrate"], nam_map)
#' test_db <- data.frame(test_db, data_name = rownames(data_kin)[ind_map],
#' stringsAsFactors = FALSE)
#' nona <- which(complete.cases(test_db))
#' kin_data <- KSR.list(test_db[nona, c(3, 2)])
#' fam <- list(akt = c("P31749", "P31751"))
#' kin_data_fam <- KSR.list(test_db[nona, c(3, 2)], kinasefamilies = fam)
#' kin_data_fam_exc <- KSR.list(test_db[nona, c(3, 2)], kinasefamilies = fam,
#'  exclusive = TRUE)
#'  
#' # get data for exclusive substrates
#' scaled_ind <- grep("scaled", colnames(data_kin))
#' substrate_profiles <- lapply(kin_data_fam_exc, 
#'                              function(x){data_kin[unlist(x), scaled_ind]})
#' # can only cluster things with more than 2 profiles
#' havesub2 <- which(sapply(substrate_profiles, nrow) > 2)
#' 
#' kin_clust <- kinclust(data_list = substrate_profiles, 
#' nboot = 100, 
#' method.dist = "euclidean",
#' method.hclust = "average")
#' 
#' # get data for all substrates
#' scaled_ind <- grep("scaled", colnames(data_kin))
#' substrate_profiles <- lapply(kin_data_fam, 
#'                              function(x){data_kin[unlist(x), scaled_ind]})
#' havesub <- match(names(kin_data_fam_exc)[havesub2], names(kin_data_fam))
#' kin_clust_all <- kinclust(data_list = substrate_profiles[havesub], 
#' nboot = 100, 
#' method.dist = "euclidean",
#' method.hclust = "average")
#'
#' expand_all_list <- clust.expand(kin_clust, kin_clust_all)

clust.expand <- function(kin_clust, kin_clust_all, thre = 0.95,
                         index = 1:length(kin_clust[[1]]), diff_list = NULL){
  # reorder using names
  ord <- match(names(kin_clust[[3]]), names(kin_clust_all[[3]]))
  kin_clust_all <- lapply(c(1:4), function(x){kin_clust_all[[x]][ord]})
  clust_class <- kin_clust[[3]]
  data <- kin_clust_all[[4]]
  expand_clust_list <- list()
  noclust <- numeric()
  for (i in 1:length(clust_class)) {
    # find all possible cluster cores
    imp_clu <- names(table(clust_class[[i]])[-1])
    clust_mem <- lapply(imp_clu, function(x){
      names(clust_class[[i]])[which(clust_class[[i]] == as.numeric(x))]
    })
    # clusters with 0.95 cut off when clustering all
    kin_pvclust <- pvclust::pvpick(kin_clust_all[[1]][[i]], alpha = thre,
                                   pv = "au", type = "geq", max.only = T)
    kin_class <- ksrlive::pvclust.clust(kin_pvclust, t(data[[i]]))
    # if any cluster member is not in a cluster when clustering all
    # then remove it from the member list
    if (any(kin_class[unlist(clust_mem)] == 0)) {
      ind <- lapply(clust_mem, function(x){
        which(kin_class[x] == 0)
      })
      clust_mem <- lapply(c(1:length(ind)), function(x){
        if (length(ind[[x]]) != 0){
          clust_mem[[x]][-ind[[x]]]
        }else{
          clust_mem[[x]]
        }
      })
    }
    # test whether the core is differentially regulated
    if (!is.null(diff_list)) {
      if (any(!(unlist(clust_mem) %in% diff_list))) {
        ind <- lapply(clust_mem, function(x){
          which(!(unlist(x) %in% diff_list))
        })
        clust_mem <- lapply(c(1:length(ind)), function(x){
          if (length(ind[[x]]) != 0) {
            clust_mem[[x]][-ind[[x]]]
          }else{
            clust_mem[[x]]
          }
        })
      }
    }
    if (all(sapply(clust_mem, length) == 0)) {
      message(paste(i, " does not have a cluster", sep = ""))
      noclust <- append(noclust, i)
      expand_clust_list[[i]] <- NULL
      next
    }
    expand <- lapply(clust_mem, function(x){
      kin_class[unlist(x)]
    })
    expand_t <- lapply(expand, table)
    if (any(sapply(expand_t, length) == 1)) {
      ind <- which(sapply(expand_t, length) == 1)
      cluster <- as.numeric(names(unlist(expand_t[ind])))
      expand_clust <- kin_class
      expand_clust[! (expand_clust %in% cluster)] <- 0
      for (j in c(1:length(cluster))){
        expand_clust[expand_clust == cluster[j]] <- j
      }
    }else{
      warning("Core in different clusters")
      warning(i)
      cluster <- as.numeric(names(expand_t))
      expand_clust <- kin_class
      expand_clust[! (expand_clust %in% cluster)] <- 0
      expand_clust[! (expand_clust == 0)] <- 1
    }
    expand_clust_list[[i]] <- expand_clust
  }
  outlist <- list(expand_clust_list = expand_clust_list, noclust = noclust)
  return(outlist)
}

#' Do pvclust clustering for data in list 
#'
#' \code{kinclust} returns list of clustering results for a given kinase
#'substrate list
#'
#' This is the main function of the package. The function kinclust takes a list of kinase substrate relationships, 
#' runs clustering using pvclust on every list and returns a list with
#' the resulting pvclust objects and cluster assignments. 
#'
#' @param data_list list containing time course data for substrates of a kinase 
#' where rows correspond to substrates and a column to the observation at a 
#' time point
#' @param cl cluster object created by the parallel package
#' @param ... arguments used in pvclust
#' @return result list of pvclust objects
#' @return most_stable list of aus used for cluster determination
#' @return clustering list of named vectors with cluster assignment
#' @return data list of data frames used for clustering
#'
#' @examples
#' data(datalist)
#' kin_clust <- kinclust(data_list = data_list, 
#' nboot = 1000, 
#' method.dist = "euclidean",
#' method.hclust = "average")
#'
#'\dontrun{
#' # using a cluster
#' library(parallel)
#' cl <- makeCluster(detectCores())
#' kin_clust <- kinclust(data_list = data_list, 
#' cl = cl,
#' nboot = 1000, 
#' method.dist = "euclidean",
#' method.hclust = "average",
#' init.rand=TRUE, iseed=21)
#' stopCluster(cl)
#' }
kinclust <- function(data, kin_list, cl = NULL, ...){
  # insert check for complete data
  data_list <- lapply(kin_list, function(x){data[unlist(x), ]})
  # can only cluster things with more than 2 profiles
  havesub2 <- which(sapply(data_list, nrow) > 2)
  if(length(havesub2) > 0) {
    data_list <- data_list[havesub2]
    clust_class_all <- list()
    result_all <- list()
    mostab_all <- list()
    data_save <- list()
    #set_seed(seed)
    for (i in 1:length(data_list)) {
      print(i)
      data_save[[i]] <- data_list[[i]]
      data_t <- t(data_save[[i]])
      if (!is.null(cl)) {
        result_all[[i]] <- pvclust::parPvclust(cl = cl, data = data_t, 
                                               r = seq(0.5, 1.5, by = 0.25), 
                                               ...)
      }else{
        result_all[[i]] <- pvclust::pvclust(data = data_t, 
                                            r = seq(0.5, 1.5, by = 0.25),
                                            ...)
      }
      mostab_all[[i]] <- ksrlive::most.stable(result_all[[i]])
      pvclust_pick <- pvclust::pvpick(result_all[[i]], alpha = mostab_all[[i]],
                                      pv = "au", type = "geq", max.only = TRUE)
      clust_class_all[[i]] <- ksrlive::pvclust.clust(pvclust_pick, data_t)
    }
    names(clust_class_all) <- names(data_list)
  }else{
    print("Not enough substrates for clustering, please ensure that more than two 
          substrates are available.")
  }
  
  outlist <- list(result = result_all, most_stable = mostab_all,
                  clustering = clust_class_all, data = data_save)
  return(outlist)
}
