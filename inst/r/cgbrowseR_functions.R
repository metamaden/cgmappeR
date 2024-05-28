# Convenience functions for cgBrowseR


getAllDNseq <- function(windowrange.start,
                     windowrange.end,
                     windowrange.chr){
  
  seqall <- getSeq(BSgenome.Hsapiens.UCSC.hg19,
                   start=windowrange.start,
                   end=windowrange.end,
                   names=windowrange.chr)
  xseq.list <- unlist(strsplit(as.character(seqall),""))
  dndf <- data.frame(dnseq=paste0(xseq.list[1:length(xseq.list)-1],xseq.list[2:length(xseq.list)]),
                     seq.coor.start=seq(1,length(xseq.list)-1,1),
                     chr.coor.start=c(windowrange.start,windowrange.start+seq(1,length(xseq.list)-2,1)),
                     chr.coor=paste0(windowrange.chr,
                                     ":",
                                     c(windowrange.start,windowrange.start+seq(1,length(xseq.list)-2,1)),
                                     "-",
                                     c(windowrange.start,windowrange.start+seq(1,length(xseq.list)-2,1))+1),
                     stringsAsFactors = FALSE)
  dndf$val <- ifelse(dndf$dnseq %in% c("CG","GC"),1,0)
  return(dndf)
}

getCGtable <- function(windowrange.start,
                       windowrange.end,
                       windowrange.chr,
                       dndf,
                      epic450anno){

  # cg table to be returned in output
  ancpg <- epic450anno[epic450anno$chr==windowrange.chr &
                         (epic450anno$pos>=windowrange.start & 
                            epic450anno$pos<=windowrange.end),]
  
  dndf.cg <- dndf[dndf$val==1,c(1:4)]
  dndf.cg$cpg.id <- NA;dndf.cg$array <- NA
  
  x <- dndf.cg[dndf.cg$dnseq=="CG",]; 
  x <- x[x$chr.coor.start %in% ancpg$pos,]
  if(nrow(x)>0){
    for(i in 1:nrow(x)){
      posi <- x$chr.coor.start[i]
      cpgi <- ancpg[ancpg$pos==posi,]$Name
      arrayi <- ancpg[ancpg$pos==posi,]$array
      
      x[x$chr.coor.start==posi,]$cpg.id <- cpgi
      x[x$chr.coor.start==posi,]$array <- ancpg[cpgi,]$array 
    }
  }
  
  dndf.cpg <- x
  
  x <- dndf.cg[dndf.cg$dnseq=="GC",]; xpos1 <- x$chr.coor.start+1
  x <- x[which(xpos1 %in% ancpg$pos),]
  if(nrow(x)>0){
    for(i in 1:nrow(x)){
      posi <- x$chr.coor.start[i]+1
      cpgi <- ancpg[ancpg$pos==posi,]$Name
      arrayi <- ancpg[ancpg$pos==posi,]$array
      
      x[x$chr.coor.start==posi-1,]$cpg.id <- cpgi
      x[x$chr.coor.start==posi-1,]$array <- ancpg[cpgi,]$array 
    }
  }
  dndf.cpg <- rbind(dndf.cpg,x)
  dndf.cpg <- rbind(dndf.cpg,
                    dndf.cg[!dndf.cg$chr.coor.start %in% dndf.cpg$chr.coor.start,])
  dndf.cpg <- dndf.cpg[order(dndf.cpg$chr.coor.start),]
  return(dndf.cpg)
  
}

makeDtrackInfo <- function(dfi=read.csv("methyldf_test.csv"),
                           plottype=c("a","p"),
                           dtrackname="data track 1"){
  # NOTE: need 3 columns labeled 'cpg', 'sampGroups', and 'methyl'
  load("epic450anno.rda")
  
  if("sampID" %in% colnames(dfi)){
    dfi.sgrp <- unique(dfi$sampGroups)
    dfi.cpg <- unique(dfi$cpg)
    
    dfi.gr <- data.frame(matrix(nrow=length(dfi.cpg),ncol=1))
    colnames(dfi.gr) <- "cpg"; dfi.gr$cpg <- dfi.cpg
    
    annoix <- epic450anno[epic450anno$Name %in% dfi.cpg,]
    annoix <- annoix[order(match(annoix$Name,dfi.gr$cpg)),]
    # identical(annoix$Name,dfi.gr$cpg)
    
    for(i in 1:length(dfi.sgrp)){
      dfi.gi <- dfi[dfi$sampGroups==dfi.sgrp[i],]
      pati <- unique(dfi.gi$sampID)
      
      # for each patient, cbind the data
      for(j in 1:length(pati)){
        dfi.patj <- dfi.gi[dfi.gi$sampID==pati[j],]
        dfi.gr <- cbind(dfi.gr,dfi.patj$methyl)
        colnames(dfi.gr)[length(colnames(dfi.gr))] <- paste0(dfi.sgrp[i],";",pati[j])
        
      }
    }
    dfi.gr$start <- annoix$pos; dfi.gr$chr <- annoix$chr; dfi.gr$end <- dfi.gr$start+1
    gri <- makeGRangesFromDataFrame(dfi.gr,keep.extra.columns = T)
    
    #=================================================
    # for data track, inc. option for multiple groups
    if(length(dfi.sgrp)==1){
      plottracki <- DataTrack(gri, 
                              name = dtrackname,
                              type=plottype)
    } else{
      plottracki <- DataTrack(gri, 
                              name = dtrackname,
                              groups=gsub("\\;.*","",colnames(mcols(gri))[-1]),
                              legend=T,
                              type=plottype)
    }
  }
  else{
    dfi.cpg <- unique(dfi$cpg)
    anni <- epic450anno[epic450anno$Name %in% dfi.cpg,]
    anni <- anni[order(match(anni$Name,as.character(unique(dfi$cpg)))),]
    
    dfi.sgrp <- as.character(unique(dfi$sampGroups))
    
    # init df to be turned into GRanges object
    dfi.gr <- as.data.frame(matrix(nrow=length(dfi.cpg),ncol=1+length(dfi.sgrp))); 
    colnames(dfi.gr) <- c("cpg",dfi.sgrp)
    dfi.gr$cpg <- dfi.cpg;
    # transform data so every row is a unique cpg
    for(j in 2:(ncol(dfi.gr))){
      grpj <- colnames(dfi.gr)[j]
      for(x in 1:nrow(dfi.gr)){
        cpgx <- dfi.gr$cpg[x]
        dfi.gr[x,j] <- dfi[dfi$cpg==cpgx & dfi$sampGroups==grpj,]$methyl
      }
    }; 
    # start, end, and chr for granges object
    dfi.gr$start <- NA; dfi.gr$chr <- NA
    for(a in 1:nrow(dfi.gr)){
      aa <- anni[anni$Name==dfi.gr$cpg[a],]
      dfi.gr$start[a] <- aa$pos; dfi.gr$chr[a] <- aa$chr
    }; dfi.gr$end <- dfi.gr$start+1
    # convert to granges object
    gri <- makeGRangesFromDataFrame(dfi.gr,keep.extra.columns = T) # processed GRanges track to plot
    
    #=================================================
    # for data track, inc. option for multiple groups
    if(length(dfi.sgrp)==1){
      plottracki <- DataTrack(gri, 
                              name = dtrackname,
                              type=plottype)
    } else{
      plottracki <- DataTrack(gri, 
                              name = dtrackname,
                              groups=dfi.sgrp,
                              legend=T,
                              type=plottype)
    }
  }
  
  return(plottracki)
}
