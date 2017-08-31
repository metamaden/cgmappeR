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

