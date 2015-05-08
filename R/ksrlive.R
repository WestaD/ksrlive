#' Do pvclust clustering for a list 
#'
#' \code{kinclust} returns list of clustering results for a given kinase
#'substrate list
#'
#' This is the main function of the package. The function kinclust takes a list of kinase substrate relationships, 
#' runs clustering using pvclust on every list and returns a list with
#' the resulting pvclust objects and cluster assignments. 
#'
#' @param data_list list of kinase substrate relationships
#' @param cl cluster object created by the parallel package
#' @param ... arguments used in pvclust
#' @return result list of pvclust objects
#' @return most_stable list of aus used for cluster determination
#' @return clustering list of named vectors with cluster assignment
#' @return data list of data frames used for clustering
#'
#' @examples
#' jjj
#'
KSR.list<-function(db, kinasefamilies=NULL, exclusive=FALSE){
  temp<-split(db[,1], f=db[,2])
  # out<-temp
  out.cl<-lapply(temp,unique) ## delete duplicates
  # combine kinasefamilies together
  if(is.null(kinasefamilies)){
    out.fam<-out.cl
  }else{
    out.fam<-lapply(kinasefamilies, function(x){unique(unlist(out.cl[unlist(x)]))})
    names(out.fam)<-sapply(kinasefamilies, "[[", 1)
    
    out.cl<-out.cl[-which(names(out.cl) %in% unlist(kinasefamilies))]
    out.fam<-append(out.cl, out.fam)
  }
  if(!exclusive){
    out.final<-out.fam
  }else{
    fam.df<-data.frame(sub=unlist(out.fam), kin=names(out.fam)[rep(seq_along(out.fam), lapply(out.fam, length))],stringsAsFactors=FALSE)
    substr<-split(fam.df[,2], f=fam.df[,1])
    substr.cl<-lapply(substr,unique) ## delete duplicates
    ### find exclusive substrates
    sub.kinases<-sapply(substr.cl, length)
    ## only one kinase
    onekin<-which(sub.kinases==1)
    # twokin<-which(sub.kinases==2)
    fam.df<-fam.df[fam.df[,1] %in% names(substr.cl)[onekin],]
    temp<-split(fam.df[,1], f=fam.df[,2])
    out.final<-lapply(temp,unique) ## delete duplicates
  }
  return(out.final)
}
#' Do pvclust clustering for a list 
#'
#' \code{kinclust} returns list of clustering results for a given kinase
#'substrate list
#'
#' This is the main function of the package. The function kinclust takes a list of kinase substrate relationships, 
#' runs clustering using pvclust on every list and returns a list with
#' the resulting pvclust objects and cluster assignments. 
#'
#' @param data_list list of kinase substrate relationships
#' @param cl cluster object created by the parallel package
#' @param ... arguments used in pvclust
#' @return result list of pvclust objects
#' @return most_stable list of aus used for cluster determination
#' @return clustering list of named vectors with cluster assignment
#' @return data list of data frames used for clustering
#'
#' @examples
#' jjj
#'
clust.expand <- function(kin_clust, kin_clust_all, thre = 0.95,
                         index = 1:length(kin_clust[[1]]), diff_list = NULL){
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
    kin_class <- ksrlive::pvclust.clust(kin_pvclust, data[[i]])
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
  # delete all sites that are in multiple lists
  # df <- data.frame(sites=names(unlist(expand.clust.list)),
  # clust=unlist(expand.clust.list), index=rep(seq_along(expand.clust.list),
  # sapply(expand.clust.list,length)), stringsAsFactors=F)
  # convert list to data frame
  # dup <- which(duplicated(df[,1])) # find duplicates
  # dup.all <- unlist(lapply(dup, function(x){which(df[,1]==df[x,1])}))
  # if(length(dup.all)>0){ # save which ones are duplicated in a data frame
  #     df.dup <- df[dup.all,]
  # }else{
  #     df.dup <- NULL
  # }
  
  # df <- df[-dup.all,]
  # expand.clust.list <- split(df[,2], df[,3])
  # clust.list.names <- split(df[,1], df[,3])
  # expand.clust.list <- lapply(c(1:length(expand.clust.list)),
  # function(x){setNames(expand.clust.list[[x]], clust.list.names[[x]])})
  outlist <- list(expand_clust_list, noclust)
  return(outlist)
}

#' Do pvclust clustering for a list 
#'
#' \code{kinclust} returns list of clustering results for a given kinase
#'substrate list
#'
#' This is the main function of the package. The function kinclust takes a list of kinase substrate relationships, 
#' runs clustering using pvclust on every list and returns a list with
#' the resulting pvclust objects and cluster assignments. 
#'
#' @param data_list list of kinase substrate relationships
#' @param cl cluster object created by the parallel package
#' @param ... arguments used in pvclust
#' @return result list of pvclust objects
#' @return most_stable list of aus used for cluster determination
#' @return clustering list of named vectors with cluster assignment
#' @return data list of data frames used for clustering
#'
#' @examples
#' jjj
#'
kinclust <- function(data_list, cl = NULL, ...){
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
