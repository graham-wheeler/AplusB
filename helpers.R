# Here is helpers.R

row <- function(...) {
  tags$div(class="row", ...)
}

col <- function(width, ...) {
  tags$div(class=paste0("span", width), ...)
}


threep3ABCDE<-function(truep, A, B, C, D, E, dose = NULL){
  if (!is.null(dose) & length(dose) != length(truep))
    stop("Length of 'dose' must be the same as the length of 'truep'.")
  path.mat<-prob<-ssize<-mtd<-dlt.no<-NULL
  exp<-0
  doses <- length(truep)
  mcohort <- 2 * doses
  mcplus1 <- mcohort + 1
  pmat <- as.data.frame(matrix(NA, nrow = 1, ncol = 3 * mcplus1))
  colnames(pmat) <- c("stop", paste(c("d", "tox", "ssize"),
                                    rep(1:mcohort, each = 3)), paste("d", mcplus1))
  pmat[1, 1:2] <- c(0,  1)
  pmat <- pmat[rep(seq_len(nrow(pmat)), rep(A+1, nrow(pmat))),
               ]
  pmat[, "tox 1"] <- 0:A
  pmat[, "ssize 1"] <- A
  pmat[pmat[, "tox 1"] <= D, "d 2"] <- 1
  pmat[pmat[, "tox 1"] < C, "d 2"] <- 2
  pmat[pmat[, "tox 1"] > D, "stop"] <- 1
  
  stopped.pmat <- pmat[pmat$stop == 1, ]
  if(dim(stopped.pmat)[1]>0){
    path.mat<-stopped.pmat # Edit from published threep3 - keeps record of all paths generated
    dose.mat <- stopped.pmat[, grep("d", names(stopped.pmat))]
    tox.mat <- stopped.pmat[, grep("tox", names(stopped.pmat))]
    prob <- apply(matrix(dbinom(as.matrix(tox.mat), A, truep[as.matrix(dose.mat)]),
                         nrow = nrow(dose.mat)), 1, prod, na.rm = T)
    ssize <- A * apply(!is.na(stopped.pmat[, grep("d", names(stopped.pmat))]),
                       1, sum)
    last.cohort <- apply(!is.na(stopped.pmat[, grep("d", names(stopped.pmat))]),
                         1, sum)
    last.drug.column <- paste("d", last.cohort)
    last.drug <- sapply(1:nrow(stopped.pmat), function(j) {
      stopped.pmat[j, last.drug.column[j]]
    })
    previous.drug.column <- paste("d", last.cohort - 1)
    previous.drug <- sapply(1:nrow(stopped.pmat), function(j) {
      ifelse(previous.drug.column[j] == "d 0", 0, stopped.pmat[j,
                                                               previous.drug.column[j]])
    })
    last.tox.column <- paste("tox", last.cohort)
    last.tox <- sapply(1:nrow(stopped.pmat), function(j) {
      stopped.pmat[j, last.tox.column[j]]
    })
    mtd <- rep(NA, nrow(stopped.pmat))
    mtd[last.tox <= D & previous.drug == last.drug] <- last.drug[last.tox <=
                                                                   D & previous.drug == last.drug] - 1
    mtd[last.tox <= D & previous.drug != last.drug] <- last.drug[last.tox <=
                                                                   D & previous.drug != last.drug]
    mtd[last.tox < C] <- last.drug[last.tox < C]
    mtd[last.tox > D] <- last.drug[last.tox > D] - 1
    exp <- sapply(1:doses, function(j) {
      sum(A * (stopped.pmat[, grep("d", names(stopped.pmat))] ==
                 j) * prob/ssize, na.rm = T)
    })
    dlt.no <- apply(stopped.pmat[, grep("tox", names(stopped.pmat))],
                    1, sum, na.rm = T)
  }
  for (i in 3:mcplus1) {
    cat(paste(round(100 * i/mcplus1), "% complete\n", sep = ""))
    dd <- as.character(paste("d", i))
    td <- as.character(paste("tox", i))
    sd <- as.character(paste("ssize", i))
    dc <- as.character(paste("d", i - 1))
    tc <- as.character(paste("tox", i - 1))
    sc <- as.character(paste("ssize", i - 1))
    db <- as.character(paste("d", i - 2))
    tb <- as.character(paste("tox", i - 2))
    sb <- as.character(paste("ssize", i - 2))
    pmat1 <- pmat[rep(which(pmat[, "stop"] == 0 & pmat[,tb]<=D & pmat[,tb]>=C & pmat[,dc]==pmat[,db]), each = B+1), ]
    #pmat2 <- pmat[rep(which(pmat[, "stop"] == 0 & pmat[,tb]<C), each = A+1), ]
    pmat3 <- pmat[rep(which(pmat[,"stop"] == 0 & pmat[,dc]!=pmat[,db]), each=A+1),]
    #browser()
    if(dim(pmat1)[1]>0){
      pmat1[, tc] <- 0:B
      pmat1[, sc] <- B
      pmat1[pmat1[,dc]==pmat1[,db] & pmat1[, tc] + pmat1[, tb] <= E & pmat1[, dc] + 1 <= doses, dd] <- pmat1[pmat1[,dc]==pmat1[,db] & pmat1[, tc] + pmat1[, tb] <= E & pmat1[, dc] + 1 <= doses, dc] + 1
      pmat1[pmat1[,dc]!=pmat1[,db] & pmat1[, tc]>=C & pmat1[, tc] <= D & pmat1[, dc] <= doses, dd] <- pmat1[pmat1[,dc]!=pmat1[,db] & pmat1[, tc]>=C & pmat1[, tc] <= D & pmat1[, dc] <= doses, dc]
    }else{pmat1<-NULL}
    #if(dim(pmat2)[1]>0){
    #pmat2[, tc] <- 0:A
    #pmat2[, sc] <- A
    ###pmat1[pmat1[, tc] == 1 & pmat1[, tb] == 0, dd] <- pmat1[pmat1[, tc] == 1 & pmat1[, tb] == 0, dc]
    #pmat2[pmat2[, tc] < C & pmat2[, dc] + 1 <= doses, dd] <- pmat2[pmat2[, tc] < C  & pmat2[, dc] + 1 <= doses, dc] + 1
    #pmat2[pmat2[, tc] <= D & pmat2[, tc] >= C, dd] <- pmat2[pmat2[, tc] <= D & pmat2[, tc] >= C, dc]
    #}else{pmat2<-NULL}
    if(dim(pmat3)[1]>0){
      pmat3[, tc] <- 0:A
      pmat3[, sc] <- A
      #pmat1[pmat1[, tc] == 1 & pmat1[, tb] == 0, dd] <- pmat1[pmat1[, tc] == 1 & pmat1[, tb] == 0, dc]
      pmat3[pmat3[, tc] < C & pmat3[, dc] + 1 <= doses, dd] <- pmat3[pmat3[, tc] < C  & pmat3[, dc] + 1 <= doses, dc] + 1
      pmat3[pmat3[, tc] <= D & pmat3[, tc] >= C, dd] <- pmat3[pmat3[, tc] <= D & pmat3[, tc] >= C, dc]
    }else{pmat3<-NULL}
    
    #pmat<-rbind(pmat1,pmat2,pmat3)
    pmat<-rbind(pmat1,pmat3)
    
    ###        pmat[pmat[, tc] == 1 & pmat[, tb] == 1 & pmat[, dc] - 1 >= 1, dd] <- pmat[pmat[, tc] == 1 & pmat[, tb] == 1, dc] - 1
    #        pmat[pmat[, tc] == 1 & pmat[, "desc"] == 0 & pmat[, tb] ==
    #            1 & pmat[, dc] - 1 >= 1, "desc"] <- 1
    ###        pmat[pmat[, tc] > 1 & pmat[, dc] - 1 >= 1, dd] <- pmat[pmat[, tc] > 1 & pmat[, dc] - 1 >= 1, dc] - 1
    #        pmat[pmat[, tc] > 1 & pmat[, dc] - 1 >= 1, "desc"] <- 1
    ###pmat<-rbind(pmat1,pmat2)
    # HERE NOW!
    excluding.dd <- names(pmat)[grepl("d ", names(pmat)) &
                                  names(pmat) != dd]
    #        cnt <- apply(pmat[!is.na(pmat[, dd]), dd] == pmat[!is.na(pmat[,
    #            dd]), excluding.dd], 1, sum, na.rm = T)
    cnt <- ifelse(sum(as.numeric(!is.na(pmat[, dd])))==0, 0, apply(pmat[!is.na(pmat[, dd]), dd] == pmat[!is.na(pmat[,
                                                                                                                    dd]), excluding.dd], 1, sum, na.rm = T)  )
    pmat[!is.na(pmat[, dd]), dd][cnt > 1] <- NA
    pmat[is.na(pmat[, dd]), "stop"] <- 1
    stopped.pmat <- pmat[pmat$stop == 1, ]
    if(dim(stopped.pmat)[1]>0){
      path.mat<-rbind(path.mat,stopped.pmat) # Edit from published threep3 - keeps record of all paths generated
      dose.mat <- stopped.pmat[, grep("d", names(stopped.pmat))[1:(i -
                                                                     1)]]
      tox.mat <- stopped.pmat[, grep("tox", names(stopped.pmat))[1:(i -
                                                                      1)]]
      ssize.mat <- stopped.pmat[, grep("ssize", names(stopped.pmat))[1:(i -
                                                                          1)]]
      prob.new <- apply(matrix(dbinom(as.matrix(tox.mat), as.matrix(ssize.mat),
                                      truep[as.matrix(dose.mat)]), nrow = nrow(dose.mat)),
                        1, prod, na.rm = T)
      prob <- c(prob, prob.new)
      
      ssize.new <- apply(ssize.mat,1,sum)#(rep(3 * (i - 1), nrow(stopped.pmat))
      ssize <- c(ssize, ssize.new)
      last.drug <- stopped.pmat[, dc]
      previous.drug <- stopped.pmat[, db]
      last.tox <- stopped.pmat[, tc]
      previous.tox <-stopped.pmat[,tb]
      last.ssize <- stopped.pmat[, sc]
      previous.ssize <-stopped.pmat[,sb]
      mtd.new <- rep(NA, nrow(stopped.pmat))
      #browser()
      #HERE!
      
      mtd.new[last.tox < C & last.drug!=previous.drug] <- last.drug[last.tox < C & last.drug!=previous.drug]
      mtd.new[last.tox + previous.tox <= E & last.drug==previous.drug] <- last.drug[last.tox + previous.tox <= E & last.drug==previous.drug]
      #mtd.new[last.tox == 1 & previous.drug == last.drug] <- last.drug[last.tox ==
      #1 & previous.drug == last.drug] - 1
      #mtd.new[last.tox == 1 & previous.drug != last.drug] <- last.drug[last.tox ==
      #1 & previous.drug != last.drug]
      mtd.new[last.tox > D & last.drug!=previous.drug] <- last.drug[last.tox > D & last.drug!=previous.drug] - 1
      mtd.new[last.tox + previous.tox > E & last.drug==previous.drug] <- last.drug[last.tox + previous.tox > E & last.drug==previous.drug] - 1
      
      mtd <- c(mtd, mtd.new)
      #exp <- exp + sapply(1:doses, function(j) {
      #sum(3 * (stopped.pmat[, grep("d", names(stopped.pmat))] == j) * prob.new/ssize.new, na.rm = T)
      #})
      exp <- exp + sapply(1:doses, function(j) { sum((stopped.pmat[, grep("ssize", names(stopped.pmat))])*(stopped.pmat[, grep("d", names(stopped.pmat))] == j) * prob.new/ssize.new, na.rm = T) })
      dlt.no <- c(dlt.no, apply(stopped.pmat[, grep("tox",
                                                    names(stopped.pmat))], 1, sum, na.rm = T))
    }
  }
  path.mat<-cbind(path.mat,prob,ssize,mtd) # Edit from published threep3 - keeps record of all paths generated
  obj <- list(prob = prob, ssize = ssize, mtd = mtd, exp = exp,
              dlt.no = dlt.no, truep = truep, dose = dose, path.mat = path.mat)
  class(obj) <- "threep3"
  return(obj)
}

threep3ABCDE.desc<-function (truep, A, B, C, D, E, dose = NULL)
{
  if (!is.null(dose) & length(dose) != length(truep))
    stop("Length of 'dose' must be the same as the length of 'truep'.")
  path.mat<-prob<-ssize<-mtd<-dlt.no<-NULL
  exp<-0
  doses <- length(truep)
  mcohort <- 2 * doses
  mcplus1 <- mcohort + 1
  pmat <- as.data.frame(matrix(NA, nrow = 1, ncol = 3 * mcplus1 +
                                 1))
  colnames(pmat) <- c("stop", "desc", paste(c("d", "tox", "ssize"),
                                            rep(1:mcohort, each = 3)), paste("d", mcplus1))
  pmat[1, 1:3] <- c(0, 0, 1)
  pmat <- pmat[rep(seq_len(nrow(pmat)), rep(A+1, nrow(pmat))),
               ]
  pmat[, "tox 1"] <- 0:A
  pmat[,"ssize 1"] <- A
  pmat[pmat[, "tox 1"] <= D, "d 2"] <- 1
  pmat[pmat[, "tox 1"] < C, "d 2"] <- 2
  pmat[pmat[, "tox 1"] > D, "stop"] <- 1
  pmat[pmat[, "tox 1"] > D, "desc"] <- 1
  stopped.pmat <- pmat[pmat$stop == 1, -2]
  if(dim(stopped.pmat)[1]>0){
    path.mat<-stopped.pmat # Edit from published threep3 - keeps record of all paths generated
    dose.mat <- stopped.pmat[, grep("d", names(stopped.pmat))]
    tox.mat <- stopped.pmat[, grep("tox", names(stopped.pmat))]
    prob <- apply(matrix(dbinom(as.matrix(tox.mat), A, truep[as.matrix(dose.mat)]),
                         nrow = nrow(dose.mat)), 1, prod, na.rm = T)
    ssize <- A * apply(!is.na(stopped.pmat[, grep("d", names(stopped.pmat))]),
                       1, sum)
    last.cohort <- apply(!is.na(stopped.pmat[, grep("d", names(stopped.pmat))]),
                         1, sum)
    last.drug.column <- paste("d", last.cohort)
    last.drug <- sapply(1:nrow(stopped.pmat), function(j) {
      stopped.pmat[j, last.drug.column[j]]
    })
    previous.drug.column <- paste("d", last.cohort - 1)
    previous.drug <- sapply(1:nrow(stopped.pmat), function(j) {
      ifelse(previous.drug.column[j] == "d 0", 0, stopped.pmat[j,
                                                               previous.drug.column[j]])
    })
    last.tox.column <- paste("tox", last.cohort)
    last.tox <- sapply(1:nrow(stopped.pmat), function(j) {
      stopped.pmat[j, last.tox.column[j]]
    })
    mtd <- rep(NA, nrow(stopped.pmat))
    mtd[last.tox <=D & previous.drug == last.drug] <- last.drug[last.tox <=
                                                                  D & previous.drug == last.drug] - 1
    mtd[last.tox <= D & previous.drug != last.drug] <- last.drug[last.tox <=
                                                                   D & previous.drug != last.drug]
    mtd[last.tox < C] <- last.drug[last.tox < C]
    mtd[last.tox > D] <- last.drug[last.tox > D] - 1
    exp <- sapply(1:doses, function(j) {
      sum(A * (stopped.pmat[, grep("d", names(stopped.pmat))] ==
                 j) * prob/ssize, na.rm = T)
    })
    dlt.no <- apply(stopped.pmat[, grep("tox", names(stopped.pmat))],
                    1, sum, na.rm = T)
  }
  for (i in 3:mcplus1) {
    cat(paste(round(100 * i/mcplus1), "% complete\n", sep = ""))
    dd <- as.character(paste("d", i))
    td <- as.character(paste("tox", i))
    sd <- as.character(paste("ssize", i))
    dc <- as.character(paste("d", i - 1))
    tc <- as.character(paste("tox", i - 1))
    sc <- as.character(paste("ssize", i - 1))
    db <- as.character(paste("d", i - 2))
    tb <- as.character(paste("tox", i - 2))
    sb <- as.character(paste("ssize", i - 2))
    pmat1 <- pmat[rep(which(pmat[, "stop"] == 0 & pmat[, tb] >= C & pmat[,tb]<=D & pmat[, "desc"] == 0 & pmat[,dc]==pmat[,db]), each = B+1), ]
    pmat2 <- pmat[rep(which(pmat[, "stop"] == 0 & pmat[, "desc"] == 0 & pmat[,dc]!=pmat[,db]), each = A+1),]
    if(dim(pmat1)[1]>0){
      pmat1[, tc] <- 0:B
      pmat1[, sc] <- B
      pmat1[pmat1[, tc] + pmat1[,tb] <=E & pmat1[,dc]+1 <=doses & pmat1[,"desc"]==0,dd]<-pmat1[pmat1[, tc] + pmat1[,tb] <=E & pmat1[,dc]+1 <=doses & pmat1[,"desc"]==0,dc]+1
      pmat1[pmat1[, tc] + pmat1[,tb] >E & pmat1[,dc] -1>=1 & pmat1[,"desc"]==0,dd]<-pmat1[pmat1[, tc] + pmat1[,tb] >E & pmat1[,dc] -1 >=1 & pmat1[,"desc"]==0,dc]-1
      pmat1[pmat1[, tc] + pmat1[,tb] >E & pmat1[,dc] -1>=1 & pmat1[,"desc"]==0,"desc"]<-1
    }else{pmat1<-NULL}
    if(dim(pmat2)[1]>0){
      pmat2[,tc] <- 0:A
      pmat2[, sc] <- A
      pmat2[pmat2[,tc] < C & pmat2[,dc]+1<=doses & pmat2[,"desc"]==0, dd]<-pmat2[pmat2[,tc] < C & pmat2[,dc]+1<=doses & pmat2[,"desc"]==0, dc]+1
      pmat2[pmat2[,tc] >= C & pmat2[,tc] <=D & pmat2[,dc]<=doses & pmat2[,"desc"]==0, dd]<-pmat2[pmat2[,tc] >= C & pmat2[,tc]<=D & pmat2[,dc]<=doses & pmat2[,"desc"]==0, dc]
      pmat2[pmat2[,tc] > D & pmat2[,dc]-1>=1 & pmat2[,"desc"]==0, dd]<-pmat2[pmat2[,tc] > D & pmat2[,dc]-1>=1 & pmat2[,"desc"]==0, dc]-1
      pmat2[pmat2[,tc] > D & pmat2[,dc]<=doses & pmat2[,"desc"]==0, "desc"]<-1
    }else{pmat2<-NULL}
    
    pmat3 <- pmat[rep(which(pmat[, "stop"] == 0 & pmat[, "desc"] == 1), each=B+1),]
    if(dim(pmat3)[1]>0){
      
      pmat3[,tc]<-0:B
      pmat3[,sc]<-B
      otherc.drug<- sapply(1:nrow(pmat3), function(j) which(pmat3[,grep("d ", names(pmat3))][j,1:which(colnames(pmat3[,grep("d ", names(pmat3))])==db)]==pmat3[j,dc]))
      for(j in 1:length(otherc.drug)){if(length(otherc.drug[[j]])==0){otherc.drug[[j]]<-NA}}
      otherc.drug<-unlist(otherc.drug)
      for(j in 1:length(otherc.drug)){
        #if(!is.na(otherc.drug[j])==TRUE){
        #   pmat3[pmat3[,tc] + pmat3[,as.character(paste("tox", otherc.drug[j]))] > E & pmat3[,dc] - 1>=1, dd] <- pmat3[pmat3[,tc] + pmat3[,as.character(paste("tox", otherc.drug[j]))] > E & pmat3[,dc] - 1>=1, dc] - 1}
        #}
        
        if(!is.na(otherc.drug[j])==TRUE){
          if(pmat3[j,tc]+pmat3[j, as.character(paste("tox", otherc.drug[j]))] > E & pmat3[j,dc] - 1 >= 1){pmat3[j,dd]<-pmat3[j, dc] - 1}
        }
      }
      #browser()
    }else{pmat3<-NULL}
    
    pmat<-rbind(pmat1,pmat2,pmat3)
    
    #if(i!=3){
    #da <- as.character(paste("d", i - 3))
    #ta <- as.character(paste("tox", i - 3))
    #sa <- as.character(paste("ssize", i - 3))
    #pmat3<-pmat[pmat[, "stop"] == 0 & pmat[, "desc"] == 1,]
    #if(dim(pmat3)[1]!=0){
    #otherc.drug<- sapply(1:nrow(pmat3), function(j) which(pmat3[,grep("d ",names(pmat3))][j,1:max(pmat3[,da])]==pmat3[j,db]))#column number
    #for(j in 1:length(otherc.drug)){if(length(otherc.drug[[j]])==0){otherc.drug[[j]]<-NA}}
    #otherc.drug<-unlist(otherc.drug)
    #pmat4<-pmat5<-NULL
    #for(j in 1:length(otherc.drug)){
    #           ifelse(!is.na(otherc.drug[j])==TRUE, ifelse(pmat3[j,tc] + pmat3[j,as.character(paste("tox", otherc.drug[j]))] > E, pmat4<-rbind(pmat4,pmat3[j,]), pmat5<-rbind(pmat5,pmat3[j,])), pmat5<-rbind(pmat5,pmat3[j,]))
    #}
    #            #browser()
    #pmat4null<-as.numeric(is.null(pmat4))+as.numeric(is.null(dim(pmat4)[1]))
    #dimpmat4null<-length(as.numeric(dim(pmat4)[1]))
    #pmat4 <- pmat4[rep(which(pmat4[, "stop"] == 0 & pmat4[, "desc"] == 1), each=B+1),]
    #if(pmat4null==0 & dimpmat4null>0){
    #pmat4[, tc] <- 0:B
    #pmat4[, sc] <- B
    #otherc.drug<- sapply(1:nrow(pmat3), function(j) which(pmat3[,grep("d ",names(pmat3))][j,1:max(pmat3[,db])]==pmat3[j,dc]))#column number
    #for(j in 1:length(otherc.drug)){if(length(otherc.drug[[j]])==0){otherc.drug[[j]]<-NA}}
    #otherc.drug<-unlist(otherc.drug)
    #                for(j in 1:length(otherc.drug)){
    #if(!is.na(otherc.drug[j])==TRUE){pmat4[pmat4[,tc] + pmat4[,as.character(paste("tox", otherc.drug[j]))] > E & pmat4[,dc] - 1>=1, dd] <- pmat4[pmat4[,tc] + pmat4[,as.character(paste("tox", otherc.drug[j]))] > E & pmat4[,dc] - 1>=1, dc] - 1}
    #}
    #pmat3[pmat3[,tc] + pmat[]]<-pmat3[]
    #}else{pmat4<-NULL}
    #pmat<-rbind(pmat1,pmat2, pmat4, pmat5)
    #}}else{pmat<-rbind(pmat1,pmat2)}
    
    
    #pmat[pmat[, tc] == 0 & pmat[, "desc"] == 0 & pmat[, dc] + 1 <= doses, dd] <- pmat[pmat[, tc] == 0 & pmat[,
    #"desc"] == 0 & pmat[, dc] + 1 <= doses, dc] + 1
    #pmat[pmat[, tc] == 1 & pmat[, "desc"] == 0 & pmat[, tb] == 0, dd] <- pmat[pmat[, tc] == 1 & pmat[, "desc"] ==
    #0 & pmat[, tb] == 0, dc]
    #pmat[pmat[, tc] == 1 & pmat[, "desc"] == 0 & pmat[, tb] == 1 & pmat[, dc] - 1 >= 1, dd] <- pmat[pmat[, tc] ==
    #1 & pmat[, "desc"] == 0 & pmat[, tb] == 1, dc] - 1
    #pmat[pmat[, tc] == 1 & pmat[, "desc"] == 0 & pmat[, tb] == 1 & pmat[, dc] - 1 >= 1, "desc"] <- 1
    #pmat[pmat[, tc] > 1 & pmat[, dc] - 1 >= 1, dd] <- pmat[pmat[, tc] > 1 & pmat[, dc] - 1 >= 1, dc] - 1
    #pmat[pmat[, tc] > 1 & pmat[, dc] - 1 >= 1, "desc"] <- 1
    
    excluding.dd <- names(pmat)[grepl("d ", names(pmat)) & names(pmat) != dd]
    cnt <- apply(pmat[!is.na(pmat[, dd]), dd] == pmat[!is.na(pmat[, dd]), excluding.dd], 1, sum, na.rm = T)
    pmat[!is.na(pmat[, dd]), dd][cnt > 1] <- NA
    pmat[is.na(pmat[, dd]), "stop"] <- 1
    stopped.pmat <- pmat[pmat$stop == 1, -2]
    if(dim(stopped.pmat)[1]>0){
      path.mat<-rbind(path.mat,stopped.pmat) # Edit from published threep3 - keeps record of all paths generated
      #browser()
      dose.mat <- stopped.pmat[, grep("d ", names(stopped.pmat))[1:(i - 1)]]
      tox.mat <- stopped.pmat[, grep("tox", names(stopped.pmat))[1:(i - 1)]]
      ssize.mat <- stopped.pmat[, grep("ssize", names(stopped.pmat))[1:(i - 1)]]
      prob.new <- apply(matrix(dbinom(as.matrix(tox.mat), as.matrix(ssize.mat),
                                      truep[as.matrix(dose.mat)]), nrow = nrow(dose.mat)),1, prod, na.rm = T)
      prob <- c(prob, prob.new)
      ssize.new <- apply(ssize.mat,1,sum)#rep(3 * (i - 1), nrow(stopped.pmat))
      ssize <- c(ssize, ssize.new)
      last.drug <- stopped.pmat[, dc]
      previous.drug <- stopped.pmat[, db]
      otherc.drug <- sapply(1:nrow(stopped.pmat), function(j) which(stopped.pmat[,grep("d ", names(stopped.pmat))][j,1:which(colnames(stopped.pmat[,grep("d ", names(stopped.pmat))])==db)]==stopped.pmat[j,dc]))
      #otherc.drug <- sapply(1:nrow(stopped.pmat), function(j) which(stopped.pmat[,grep("d",names(stopped.pmat))][j,1:max(stopped.pmat[,db])]==stopped.pmat[j,dc]))#column number
      for(j in 1:length(otherc.drug)){if(length(otherc.drug[[j]])==0){otherc.drug[[j]]<-NA}}
      otherc.drug<-unlist(otherc.drug)
      last.tox <- stopped.pmat[, tc]
      previous.tox <- stopped.pmat[, tb]
      otherc.tox<-rep(NA,length(otherc.drug))
      for(j in 1:length(otherc.tox)){ifelse(is.na(otherc.drug[j])==TRUE,otherc.tox[j]<-NA,otherc.tox[j]<-stopped.pmat[j,grep("tox", names(stopped.pmat))][,otherc.drug[j]])}
      #otherc.drug<- sapply(1:nrow(stopped.pmat), function(j) which(stopped.pmat[,grep("tox",names(stopped.pmat))][j,1:max(stopped.pmat[,dc])-1]==stopped.pmat[j,dc]))
      #for(i in 1:length(otherc.drug)){if(length(otherc.drug[[i]])==0){otherc.drug[[i]]<-NA}}
      #otherc.drug<-unlist(otherc.drug)
      last.ssize <- stopped.pmat[, sc]
      previous.ssize <- stopped.pmat[, sb]
      mtd.new <- rep(NA, nrow(stopped.pmat))
      #browser()
      mtd.new[last.tox < C & last.drug!=previous.drug] <- last.drug[last.tox < C & last.drug!=previous.drug]
      mtd.new[last.tox + previous.tox > E & previous.drug == last.drug] <- last.drug[last.tox + previous.tox > E & previous.drug == last.drug] - 1
      mtd.new[last.tox + previous.tox <= E & previous.drug == last.drug] <- last.drug[last.tox + previous.tox <= E & previous.drug == last.drug]
      mtd.new[last.tox > D & previous.drug != last.drug] <- last.drug[last.tox > D & previous.drug != last.drug] - 1
      for(j in 1:length(otherc.tox)){
        if(!is.na(otherc.tox[j])==TRUE){mtd.new[last.tox + otherc.tox <= E & last.drug == dose.mat[,otherc.drug[j]]] <- last.drug[last.tox + otherc.tox <= E & last.drug == dose.mat[,otherc.drug[j]]]}
      }
      #mtd.new[last.tox  1 & previous.drug != last.drug] <- last.drug[last.tox == 1 & previous.drug != last.drug]
      mtd <- c(mtd, mtd.new)
      exp <- exp + sapply(1:doses, function(j) {
        sum((stopped.pmat[,grep("ssize", names(stopped.pmat))]) * (stopped.pmat[, grep("d", names(stopped.pmat))] == j) * prob.new/ssize.new, na.rm = T) })
      dlt.no <- c(dlt.no, apply(stopped.pmat[, grep("tox",
                                                    names(stopped.pmat))], 1, sum, na.rm = T))
    }
  }
  path.mat<-cbind(path.mat,prob,ssize,mtd) # Edit from published threep3 - keeps record of all paths generated
  obj <- list(prob = prob, ssize = ssize, mtd = mtd, exp = exp,
              dlt.no = dlt.no, truep = truep, dose = dose, path.mat = path.mat)
  class(obj) <- "threep3"
  return(obj)
}


print.threep3<-function(x,tox.cutpoints=NULL,dose=NULL,...){
  if(is.null(dose)){
    dose<-1:length(x$truep)
  }
  # average sample size
  n.average<-weighted.mean(x$ssize,x$prob)
  n.min<-min(x$ssize[x$prob>0])
  n.max<-max(x$ssize[x$prob>0])
  tab0<-cbind(n.average,n.min,n.max)
  rownames(tab0)<-"Sample size"
  colnames(tab0)<-c("Mean","Minimum","Maximum")
  exp<-c(NA,x$exp)
  rec<-xtabs(x$prob~x$mtd)
  tab<-signif(rbind(exp,rec),3)
  rownames(tab)<-c("Experimentation proportion","Recommendation proportion")
  colnames(tab)<-c("No Dose",dose)
  names(dimnames(tab))<-c("","Doses")
  if(is.null(tox.cutpoints)){
    tox.cutpoints<-seq(0,1,by=0.2)
  } else {
    tox.cutpoints<-unique(c(0,tox.cutpoints,1))
  }
  exp.tox<-xtabs(x$exp~cut(x$truep,tox.cutpoints,include.lowest=T))
  rec.tox<-xtabs(rec[-1]~cut(x$truep,tox.cutpoints,include.lowest=T))
  tab2<-signif(rbind(exp.tox,rec.tox),3)
  rownames(tab2)<-c("Experimentation proportion","Recommendation proportion*")
  names(dimnames(tab2))<-c("","Probability of DLT")
  
  #print(tab0)
  #cat("\n")
  #print(tab)
  #cat("\n")
  #print(tab2)
  #cat("\n * Amongst those trials that recommend an MTD\n")
  return(list(tab0=tab0, tab=tab, tab2=tab2))
}


plot.threep3<-function(x,file=NULL,...){
  dose<-if(is.null(x$dose)) 1:length(x$truep) else x$dose
  dose.label<-if(is.null(x$dose)) "Dose levels" else "Dose"
  
  # sample size
  n.threep3<-x$ssize
  df.n.threep3<-data.frame(n=n.threep3,weight=x$prob)
  #a<-ggplot()+stat_bin(aes(x=n,y=100*..density..,weight=weight),data=df.n.threep3,binwidth=1)+xlab("Sample size")+ylab("Percent")+ggtitle("Sample size")
  a<-ggplot()+stat_count(aes(x=n,y=100*..count..,weight=weight),data=df.n.threep3, width=1)+xlab("Sample size")+ylab("Percent")+ggtitle("Sample size")
  
  # experimentation
  exp.threep3<-rep(dose,10000*x$exp)
  df.exp.threep3<-data.frame(exp=as.factor(exp.threep3))
#  b<-ggplot()+geom_histogram(aes(x=exp,y=..count../100),data=df.exp.threep3)+xlab(dose.label)+ylab("Percent")+ggtitle("Experimentation")
  b<-ggplot()+stat_count(aes(x=exp,y=..count../100),data=df.exp.threep3)+xlab(dose.label)+ylab("Percent")+ggtitle("Experimentation")
  
  # recommendation
  rec.threep3<-dose[x$mtd]
  df.rec.threep3<-data.frame(rec=factor(rec.threep3),weight=x$prob[x$mtd!=0])
  #c<-ggplot()+geom_histogram(aes(x=rec,y=100*..count..,weight=weight),data=df.rec.threep3)+xlab(dose.label)+ylab("Percent")+ggtitle("Recommendation")
  c<-ggplot()+stat_count(aes(x=rec,y=100*..count..,weight=weight),data=df.rec.threep3)+xlab(dose.label)+ylab("Percent")+ggtitle("Recommendation")
  
  # observed DLTs
  obs.threep3<-100*x$dlt.no/x$ssize
  #bw<-max(diff(range(obs.threep3[x$prob>0]))/30,1)
  #bw<-min(max(min(diff(sort(unique(obs.threep3[x$prob>0])))),1/max(x$ssize)),1)
  bw<-1
  df.obs.threep3<-data.frame(obs=obs.threep3,weight=x$prob,bw=bw)
  df.obs.threep3<-subset(df.obs.threep3,df.obs.threep3$weight>0)
  #d<-ggplot()+geom_histogram(aes(bw=bw,x=obs,y=100*..density..*bw,weight=weight),data=df.obs.threep3,binwidth=bw)+xlab("Percentage of subjects with DLTs")+ylab("Percent")+ggtitle("DLTs")
  suppressWarnings(d<-ggplot()+stat_count(aes(bw=bw,x=obs,y=100*..count..,weight=weight),data=df.obs.threep3,width=bw)+xlab("Percentage of subjects with DLTs")+ylab("Percent")+ggtitle("DLTs"))
  
  if(!is.null(file))
    pdf(paste(file,".pdf",sep=""),...)
  grid.newpage()
  pushViewport(viewport(layout=grid.layout(2,2)))
  vplayout<-function(x,y)	viewport(layout.pos.row=x,layout.pos.col=y)
  if(!is.null(a)) print(a,vp=vplayout(1,1))	
  print(b,vp=vplayout(1,2))
  print(c,vp=vplayout(2,1))	
  suppressWarnings(print(d,vp=vplayout(2,2))	)
  if(!is.null(file))
    dev.off()
}

stats.exp.fn<-function(obj){
  mat<-obj$path.mat
  prob<-mat[,"prob"]
  ssize<-mat[,"ssize"]
  mtd<-mat[,"mtd"]
  dose.mat<-mat[,grep("^d",colnames(mat))]
  tox.mat<-mat[,grep("^tox ",colnames(mat))]
  ssize.mat<-mat[,grep("^ssize ",colnames(mat))]
  dose.mat<-dose.mat[,-dim(dose.mat)[2]]
  N<-dim(dose.mat)[1]
  M<-dim(dose.mat)[2]/2
  fun1<-function(x,y){
    out1<-ifelse(mtd[x]==y,1,0)
    out1}
  mtd.mat<-outer(1:N, 1:M, fun1)
  mtd.prob<-sapply(1:M, function(x) sum(prob[which(mtd==x)]))
  out<-mclapply(1:N, function(x) lapply(1:M, function(y) {index<-which(dose.mat[x,]==y)
  numpat<-ifelse(length(index)==0,0,sum(ssize.mat[x,index]))
  dose.cond<-numpat/ssize[x]
  prob.cond<-ifelse(length(index)==0,0,sum(tox.mat[x,index])/numpat)
  c(numpat, dose.cond, prob.cond)
  }), mc.cores=if(.Platform$OS.type == "windows") {1} else {2})
  numpat.mat<-matrix(unlist(out)[seq(1,length(unlist(out)),by=3)],nrow=N,ncol=M,byrow=T)
  dose.cond.mat<-matrix(unlist(out)[seq(2,length(unlist(out)),by=3)],nrow=N,ncol=M,byrow=T)
  prob.cond.mat<-matrix(unlist(out)[seq(3,length(unlist(out)),by=3)],nrow=N,ncol=M,byrow=T)
  exp.xj<- round(t(numpat.mat)%*%prob,4)
  prob.xj<- round(t(dose.cond.mat)%*%prob,4)
  exp.tj<-round(diag(t(prob.cond.mat)%*%(numpat.mat*prob)),4)
  ETL<-obj$truep%*%(t(mtd.mat)%*%prob)/sum(t(mtd.mat)%*%prob)
  mean.dlt<-sum(exp.tj)
  EOTR<-mean.dlt/sum(exp.xj)
  tab2<-cbind(ETL,mean.dlt,EOTR)
  rownames(tab2)<-c("Results")
  colnames(tab2)<-c("ETL","Mean No. DLTs","EOTR")
  return(list(tab2=tab2))
}


intervalAplusB<-function(w=95, A=3, B=3, C=1, E=1, deesc=FALSE, method){
  v<-w/100
  if(deesc==FALSE){
    if(C>0){
      vals<-seq(0,C-1,by=1)
      set1<-matrix(NA, ncol=3, nrow=C)
      for(i in 1:length(vals)){
        set1[i,1]<-paste(vals[i],"/",A,sep="")
        set1[i,2:3]<-round(100*unlist(binom.confint(vals[i],A,v,method=method)[,5:6]),2)
        #set1[i,2:3]<-c(round(100*qbeta((1-v)/2,vals[i],A-vals[i]+1),2),round(100*qbeta((1+v)/2,vals[i]+1,A-vals[i]),2))
      }
    }else{set1<-NULL}
    set2<-matrix(NA, ncol=3, nrow=E-C+1)
    vals<-seq(C,E,by=1)
    for(i in 1:length(vals)){
      set2[i,1]<-paste(vals[i],"/",A+B,sep="")
      set2[i,2:3]<-round(100*unlist(binom.confint(vals[i],A+B,v,method=method)[,5:6]),2)
      #set2[i,2:3]<-c(round(100*qbeta((1-v)/2,vals[i],A+B-vals[i]+1),2),round(100*qbeta((1+v)/2,vals[i]+1,A+B-vals[i]),2))
    }
    set<-rbind(set1,set2)
  }else{
    set<-matrix(NA, ncol=3, nrow=E+1)
    for(i in 0:E){
      set[i+1,1]<-paste(i,"/",A+B,sep="")
      set[i+1,2:3]<-round(100*unlist(binom.confint(i,A+B,v,method=method)[,5:6]),2)
      #set[i+1,2:3]<-c(round(100*qbeta((1-v)/2,i,A+B-i+1),2),round(100*qbeta((1+v)/2,i+1,A+B-i),2))
    }
  }
  colnames(set)<-c("Data at MTD", "Lower Bound (%)", "Upper Bound (%)")
  set
}






foo<-function(p,A,B,C,D,E){
  pp<-inv.logit(p)
  p.esc<-pbinom(C-1,A,pp) + dbinom(seq(C,D,by=1),A,pp)%*%pbinom(E-seq(C,D,by=1),B,pp)
  (p.esc-0.5)^2
}

foofind<-function(A,B,C,D,E){
  loop.fin<-0
  answer<-99
  trials<-logit(seq(1e-10,1-1e-10,length=10))
  t1<-sapply(1:length(trials), function(z) foo(trials[z], A, B, C, D, E))
  x<-range(t1)
  y<-mean(t1)
  if(isTRUE(all.equal((x/y)[1],(x/y)[2],tolerance=.Machine$double.eps^0.5))){
    #if(length(unique(sapply(1:length(trials), function(z) foo(trials[z], A, B, C, D, E))))==1){
    answer<-NA
  }else{
    while(answer==99){
      loop.fin<-loop.fin+1
      inits<-logit(runif(1,1e-10,1-1e-10))
      fooscale<-foo(inits,A,B,C,D,E)
      fooopt<-try(optim(par=inits, fn=foo, method="BFGS", A=A, B=B, C=C, D=D, E=E, control=list(fnscale=fooscale, maxit=20000, abstol=.Machine$double.eps^6)), silent=TRUE)
      answer<-ifelse(class(fooopt)=="try-error", 99, ifelse(fooopt$conv!=0,99,ifelse(fooopt$value>0.01,99,inv.logit(fooopt$par))))
      if(loop.fin>50){answer<-NA}
    }
  }
  answer
}
