#!/usr/bin/env R

# group_plot_functions.R
#
#

libv <- c('GenomicRanges', 'Gviz')
sapply(libv, library, character.only = TRUE)

get_plot_data <- function(data){
  # get_plot_data
  # gets granges for plot, from data
  #
  # details:
  # colnames seq, start, stop reserved for genomic coordinates
  # remaining columns each correspond to 1 group of data to be plotted.
  # e.g.
  # chr    start     stop dataA dataB
  # 1   1  1700500  1700501  0.80  0.81
  # 2   1  1700555  1700556  0.31  0.29
  #
  #
  if(is(data, "character")){data <- read.csv(data, header=TRUE)}
  colnamesRequired <- c("chr", "start", "stop")
  colnamesFound <- intersect(colnames(data), colnamesRequired)
  if(length(colnamesFound)<length(colnamesRequired)){
    stop("Error, not all required colnames found in data.")
  }
  dataPre <- data
  plotData <- makeGRangesFromDataFrame(
    dataPre,
    keep.extra.columns=TRUE,
    ignore.strand=TRUE,
    seqnames.field="chr",
    start.field="start",
    end.field="stop",
    na.rm=TRUE
  )
  return(plotData)
}

get_group_plot <- function(data){
  # get_group_plot
  #
  # data: path to csv or data.frame object.
  #
  # details: 
  # colnames seq, start, stop reserved for genomic coordinates
  # remaining columns each correspond to 1 group of data to be plotted.
  # e.g.
  # chr    start     stop dataA dataB
  # 1   1  1700500  1700501  0.80  0.81
  # 2   1  1700555  1700556  0.31  0.29
  #
  # example
  # get_group_plot(data) |> plotTracks(type=c('p','l','b')) # makes group ideogram plot
  #
  plotData <- get_plot_data(data)
  groupVector <- mcols(plotData) |> colnames()
  dtTrack <- DataTrack(range = plotData, 
                       genome = "hg19", 
                       name = "random data",
                       groups = groupVector)
  return(dtTrack)
}