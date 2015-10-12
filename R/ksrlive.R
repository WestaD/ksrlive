#' Identify site specific kinase substrate relationships using dynamic data.
#'
#' @description Using this package you can combine known site specific 
#' kinase substrate relationships with dynamic experimental data and determine active 
#' kinases and their substrates.
#' 
#' @author Westa Domanova
#' @docType package
#' @name ksrlive
NULL

#' Create a kinase substrate relationship list from a data frame
#'
#' \code{KSR.list} returns a list of kinase substrate relationships
#'
#' The function KSR.list creates a list of kinase substrate relationships from
#' a data frame and can combine kinase families into one list. Substrates occuring 
#' in multiple lists can be deleted. 
#'
#' @param db data frame of kinase substrate relationships with substrate 
#' identifier in the first column and kinase identifier in the second column.
#' @param kinasefamilies list of kinase identifiers that have to be combined, 
#' one list per kinase family, list will be named after first family member
#' @param exclusive logical, if TRUE only substrates exclusive to the kinase
#' will be included in the list 
#' @return list with kinase substrate relationships, with the kinase
#' identifiers as the list names
#'
#' @examples
#' data(phosphonetwork)
#' data(datalist)
#' # create identifier for substrate in database
#' test_db <- do.call(paste, 
#'                    c(phosphonetwork.df[ ,c("SUB_ACC_ID.human", "MODSITE_SEQ.human")],
#'                      sep = "_"))
#' test_db <- data.frame(substrate = test_db, 
#'                       kinase = phosphonetwork.df[ ,"KIN_ACC_ID.human"], 
#'                       stringsAsFactors = FALSE)
#' 
#' # create identifier for substrate in data
#' nam_map <- do.call(paste, 
#'                    c(data_kin[, c("Uniprot.human", "Motif.human")], sep = "_"))
#' ind_map <- match(test_db[ ,"substrate"], nam_map)
#' test_db <- data.frame(test_db, data_name = rownames(data_kin)[ind_map],
#'                       stringsAsFactors = FALSE)
#'                       
# first column has to be substrate id, second kinase id
#' kin_data <- KSR.list(test_db[, c(3, 2)]) 
#' # Akt1 and Akt2 belong to the same kinase family, combine their substrates 
#' # into one list and name the list after the first family member
#' fam <- list(akt = c("P31749", "P31751"))
#' kin_data_fam <- KSR.list(test_db[, c(3, 2)], kinasefamilies = fam)
#' 
#' # only include phosphosites appearing once
#' kin_data_fam_exc <- KSR.list(test_db[, c(3, 2)], kinasefamilies = fam,
#'                              exclusive = TRUE)

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
#' threshold. The core sites can be tested for differential regulation if a list 
#' of differentially regulated sites is included (recommended). 
#'
#' @param kin_clust list of kinase substrate relationships with only exclusive
#' substrates
#' @param kin_clust_all list of all available kinase substrate relationships
#' @param thre threshold of au value for clustering, for example au = 0.95 
#' corresponds to a pvalue of 0.05
#' @param index numeric vector of indices to perform the integration for
#' @param diff_list character vector of names of differentially regulated
#' substrates
#' @return expand_clust_list integrated clustering assignments
#' @return noclust numeric vector of indices where a cluster could not be 
#' determined
#'
#' @examples
#' data(phosphonetwork)
#' data(datalist)
#' # create identifier for substrate in database
#' test_db <- do.call(paste, 
#'                    c(phosphonetwork.df[ ,c("SUB_ACC_ID.human", "MODSITE_SEQ.human")],
#'                      sep = "_"))
#' test_db <- data.frame(substrate = test_db, 
#'                       kinase = phosphonetwork.df[ ,"KIN_ACC_ID.human"], 
#'                       stringsAsFactors = FALSE)
#' 
#' # create identifier for substrate in data
#' nam_map <- do.call(paste, 
#'                    c(data_kin[, c("Uniprot.human", "Motif.human")], sep = "_"))
#' ind_map <- match(test_db[ ,"substrate"], nam_map)
#' test_db <- data.frame(test_db, data_name = rownames(data_kin)[ind_map],
#'                       stringsAsFactors = FALSE)
#'                       
#' # first column has to be substrate id, second kinase id
#' kin_data <- KSR.list(test_db[, c(3, 2)]) 
#' # Akt1 and Akt2 belong to the same kinase family, combine their substrates 
#' # into one list and name the list after the first family member
#' fam <- list(akt = c("P31749", "P31751"))
#' kin_data_fam <- KSR.list(test_db[, c(3, 2)], kinasefamilies = fam)
#' 
#' # only include phosphosites appearing once
#' kin_data_fam_exc <- KSR.list(test_db[, c(3, 2)], kinasefamilies = fam,
#'                              exclusive = TRUE)
#'                              
#' scaled_ind <- grep("scaled", colnames(data_kin))
#' 
#' # clustering using exclusive substrates
#' kin_clust <- kinclust(data = data_kin[ , scaled_ind],
#'                       kin_list = kin_data_fam_exc,
#'                       nboot = 100,
#'                       method.dist = "euclidean",
#'                       method.hclust = "average")
#' 
#' # clustering using all substrates
#' kin_clust_all <- kinclust(data = data_kin[ , scaled_ind],
#'                           kin_list = kin_data_fam,
#'                           nboot = 100,
#'                           method.dist = "euclidean",
#'                           method.hclust = "average")
#' 
#' expand_all_list <- clust.expand(kin_clust, kin_clust_all)
#' expand_all <- expand_all_list$expand_clust_list

clust.expand <- function(kin_clust, kin_clust_all, thre = 0.95, 
                         diff_list = NULL){
  # reorder using names
  ord <- match(names(kin_clust), names(kin_clust_all))
  kin_clust_all <- kin_clust_all[ord]
  clust_class <- lapply(kin_clust, "[[", 3)
  data <- lapply(kin_clust, "[[", 4)
  expand_clust_list <- list()
  noclust <- numeric()
  for (i in 1:length(clust_class)) {
    # find all possible cluster cores
    imp_clu <- names(table(clust_class[[i]])[-1])
    clust_mem <- lapply(imp_clu, function(x){
      names(clust_class[[i]])[which(clust_class[[i]] == as.numeric(x))]
    })
    # clusters with 0.95 cut off when clustering all
    kin_pvclust <- pvclust::pvpick(kin_clust_all[[i]][[1]], alpha = thre,
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

#' Create a kinase substrate relationship list from a data frame
#'
#' \code{KSR.list} returns a list of kinase substrate relationships
#'
#' The function KSR.list creates a list of kinase substrate relationships from
#' a data frame and can combine kinase families into one list. Substrates occuring 
#' in multiple lists can be deleted. 
#'
#' @param db data frame of kinase substrate relationships with substrate 
#' identifier in the first column and kinase identifier in the second column.
#' @param kinasefamilies list of kinase identifiers that have to be combined, 
#' one list per kinase family, list will be named after first family member
#' @param exclusive logical, if TRUE only substrates exclusive to the kinase
#' will be included in the list 
#' @return list with kinase substrate relationships, with the kinase
#' identifiers as the list names
#'
#' @examples
#' data(phosphonetwork)
#' data(datalist)
#' # create identifier for substrate in database
#' test_db <- do.call(paste, 
#'                    c(phosphonetwork.df[ ,c("SUB_ACC_ID.human", "MODSITE_SEQ.human")],
#'                      sep = "_"))
#' test_db <- data.frame(substrate = test_db, 
#'                       kinase = phosphonetwork.df[ ,"KIN_ACC_ID.human"], 
#'                       stringsAsFactors = FALSE)
#' 
#' # create identifier for substrate in data
#' nam_map <- do.call(paste, 
#'                    c(data_kin[, c("Uniprot.human", "Motif.human")], sep = "_"))
#' ind_map <- match(test_db[ ,"substrate"], nam_map)
#' test_db <- data.frame(test_db, data_name = rownames(data_kin)[ind_map],
#'                       stringsAsFactors = FALSE)
#'                       
# first column has to be substrate id, second kinase id
#' kin_data <- KSR.list(test_db[, c(3, 2)]) 
#' # Akt1 and Akt2 belong to the same kinase family, combine their substrates 
#' # into one list and name the list after the first family member
#' fam <- list(akt = c("P31749", "P31751"))
#' kin_data_fam <- KSR.list(test_db[, c(3, 2)], kinasefamilies = fam)
#' 
#' # only include phosphosites appearing once
#' kin_data_fam_exc <- KSR.list(test_db[, c(3, 2)], kinasefamilies = fam,
#'                              exclusive = TRUE)


### add random data
random.data <- function(data, back_data = NULL, n = 50, random.seed = NULL){
      if(n <= nrow(data)){
            n <- 0
            message("Data has more observations than given maximum.")
            return(NULL)
      }else{
            if(!is.null(back_data)){
                  indata <- match(rownames(data), rownames(back_data))
                  back_data <- back_data[-indata, ]
                  if(!is.null(random.seed)){
                        set.seed(random.seed)
                  }
                  random_data <- back_data[sample(rownames(back_data), n - nrow(data)), ]
                  rownames(random_data) <- paste(rownames(random_data), rep("_random", nrow(random_data)), sep = "")
            }else{
                  if(!is.null(random.seed)){
                        set.seed(random.seed)
                  }
                  random_data <- t(replicate(n-nrow(data), runif(ncol(data), min = min(data), max = max(data))))
                  rownames(random_data) <- paste(c(1:nrow(random_data)), rep("_random", nrow(random_data)), sep = "")
                  colnames(random_data) <- colnames(data)
            }
            return(random_data)
      }
}

#' Create a kinase substrate relationship list from a data frame
#'
#' \code{KSR.list} returns a list of kinase substrate relationships
#'
#' The function KSR.list creates a list of kinase substrate relationships from
#' a data frame and can combine kinase families into one list. Substrates occuring 
#' in multiple lists can be deleted. 
#'
#' @param db data frame of kinase substrate relationships with substrate 
#' identifier in the first column and kinase identifier in the second column.
#' @param kinasefamilies list of kinase identifiers that have to be combined, 
#' one list per kinase family, list will be named after first family member
#' @param exclusive logical, if TRUE only substrates exclusive to the kinase
#' will be included in the list 
#' @return list with kinase substrate relationships, with the kinase
#' identifiers as the list names
#'
#' @examples
#' data(phosphonetwork)
#' data(datalist)
#' # create identifier for substrate in database
#' test_db <- do.call(paste, 
#'                    c(phosphonetwork.df[ ,c("SUB_ACC_ID.human", "MODSITE_SEQ.human")],
#'                      sep = "_"))
#' test_db <- data.frame(substrate = test_db, 
#'                       kinase = phosphonetwork.df[ ,"KIN_ACC_ID.human"], 
#'                       stringsAsFactors = FALSE)
#' 
#' # create identifier for substrate in data
#' nam_map <- do.call(paste, 
#'                    c(data_kin[, c("Uniprot.human", "Motif.human")], sep = "_"))
#' ind_map <- match(test_db[ ,"substrate"], nam_map)
#' test_db <- data.frame(test_db, data_name = rownames(data_kin)[ind_map],
#'                       stringsAsFactors = FALSE)
#'                       
# first column has to be substrate id, second kinase id
#' kin_data <- KSR.list(test_db[, c(3, 2)]) 
#' # Akt1 and Akt2 belong to the same kinase family, combine their substrates 
#' # into one list and name the list after the first family member
#' fam <- list(akt = c("P31749", "P31751"))
#' kin_data_fam <- KSR.list(test_db[, c(3, 2)], kinasefamilies = fam)
#' 
#' # only include phosphosites appearing once
#' kin_data_fam_exc <- KSR.list(test_db[, c(3, 2)], kinasefamilies = fam,
#'                              exclusive = TRUE)


clustering <- function(tightclust, data){
      clustnum <- unique(tightclust$cluster)
      clustnum <- clustnum[-which(clustnum == -1)]
      clust <- lapply(c(1:length(clustnum)), function(x){which(tightclust$cluster == x)})
      
      ### delete all clusters that are only random data
      cluster <- lapply(clust, function(x){
            as.character(na.omit(rownames(data)[x]))
      })
      big <- which(sapply(cluster, length) > 1)
      ### create cluster assignment vector
      clust_assg <- numeric(length = nrow(data))
      names(clust_assg) <- rownames(data)
      if(length(big) > 0){
            for(i in 1:length(big)){
                  clust_assg[cluster[[big[i]]]] <- i
            }
      }else{
            message("Does not have a cluster.")
      }
      return(clust_assg)
}

#' Create a kinase substrate relationship list from a data frame
#'
#' \code{KSR.list} returns a list of kinase substrate relationships
#'
#' The function KSR.list creates a list of kinase substrate relationships from
#' a data frame and can combine kinase families into one list. Substrates occuring 
#' in multiple lists can be deleted. 
#'
#' @param db data frame of kinase substrate relationships with substrate 
#' identifier in the first column and kinase identifier in the second column.
#' @param kinasefamilies list of kinase identifiers that have to be combined, 
#' one list per kinase family, list will be named after first family member
#' @param exclusive logical, if TRUE only substrates exclusive to the kinase
#' will be included in the list 
#' @return list with kinase substrate relationships, with the kinase
#' identifiers as the list names
#'
#' @examples
#' data(phosphonetwork)
#' data(datalist)
#' # create identifier for substrate in database
#' test_db <- do.call(paste, 
#'                    c(phosphonetwork.df[ ,c("SUB_ACC_ID.human", "MODSITE_SEQ.human")],
#'                      sep = "_"))
#' test_db <- data.frame(substrate = test_db, 
#'                       kinase = phosphonetwork.df[ ,"KIN_ACC_ID.human"], 
#'                       stringsAsFactors = FALSE)
#' 
#' # create identifier for substrate in data
#' nam_map <- do.call(paste, 
#'                    c(data_kin[, c("Uniprot.human", "Motif.human")], sep = "_"))
#' ind_map <- match(test_db[ ,"substrate"], nam_map)
#' test_db <- data.frame(test_db, data_name = rownames(data_kin)[ind_map],
#'                       stringsAsFactors = FALSE)
#'                       
# first column has to be substrate id, second kinase id
#' kin_data <- KSR.list(test_db[, c(3, 2)]) 
#' # Akt1 and Akt2 belong to the same kinase family, combine their substrates 
#' # into one list and name the list after the first family member
#' fam <- list(akt = c("P31749", "P31751"))
#' kin_data_fam <- KSR.list(test_db[, c(3, 2)], kinasefamilies = fam)
#' 
#' # only include phosphosites appearing once
#' kin_data_fam_exc <- KSR.list(test_db[, c(3, 2)], kinasefamilies = fam,
#'                              exclusive = TRUE)

clust.expand <- function(clust, clust_all, diff_list = NULL){
      # find all possible cluster cores
      imp_clu <- names(table(clust))
      if(any(imp_clu == 0)){
            imp_clu <- imp_clu[-which(imp_clu == 0)]
      }
      clust_mem <- lapply(imp_clu, function(x){
            names(clust)[which(clust == as.numeric(x))]
      })
      # if any cluster member is not in a cluster when clustering all
      # then remove it from the member list
      if (any(clust_all[unlist(clust_mem)] == 0)) {
            ind <- lapply(clust_mem, function(x){
                  which(clust_all[x] == 0)
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
      # if no cluster found
      if (all(sapply(clust_mem, length) == 0)) {
            message("Does not have a cluster")
            expand_clust <- NA
      }else{
            expand <- lapply(clust_mem, function(x){
                  clust_all[unlist(x)]
            })
            expand_t <- lapply(expand, table)
            cluster <- lapply(expand_t, function(x){as.numeric(names(unlist(x)))})
            expand_clust <- clust_all
            expand_clust[! (expand_clust %in% unlist(cluster))] <- 0
            hascl <- which(sapply(cluster, length) > 0)
            clustnum_vec <- c(1:length(hascl))
            clust_ind <- list()
            for(i in 1:length(hascl)){
                  clust_ind[[i]] <- which(expand_clust %in% cluster[[hascl[i]]])
            }
            for(i in 1:length(clust_ind)){
                  expand_clust[clust_ind[[i]]] <- clustnum_vec[i]
            }
            
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
  
  # reformat into more intuitive object
  le <- length(outlist[[1]])
  outlist_ref <- lapply(c(1:le), function(x){lapply(outlist, "[[", x)})
  names(outlist_ref) <- names(kin_list)[havesub2]
  return(outlist_ref)
=======
      return(expand_clust)
}
