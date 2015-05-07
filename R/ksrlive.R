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

