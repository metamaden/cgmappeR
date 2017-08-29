# gviz interactive shinyR browser

library(shiny)
library(shinythemes)
library(Gviz)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(org.Hs.eg.db)
library(BSgenome.Hsapiens.UCSC.hg19)
load("grgenes_symbols.rda");load("grgenes.rda")

ui <- fluidPage(theme=shinytheme("cerulean"),
                titlePanel(title=div(HTML(paste(h1("CpG MappeR"))))),
                sidebarPanel(
                  conditionalPanel(condition="input.conditionedPanels==2",
                                   helpText(h3("Genome Visualization")),
                                   helpText("This tab displays a genome visualization"),
                                   textInput("genesym2", label = h5("Gene Symbol"), value = ""),
                                   actionButton("loadgenecoordinates", "Load Gene Coordinates"),
                                   uiOutput("startcoor"),
                                   uiOutput("endcoor"),
                                   uiOutput("chrgene"),
                                   checkboxGroupInput("trackoptions",
                                               label=h5("Track options (default: Genome Coordinates; Ideogram; CG dinucleotides)"),
                                               choices=c("Ensembl Genes (BioMart/UCSC)" = 1,
                                                         "CpG Islands (UCSC)" = 2,
                                                         "GC Content (UCSC)" = 3)),
                                   actionButton("viewgenome", "View Genome"),
                                   downloadButton('dncgtable_download.csv', 'Download CG Coordinates')
                                   
                  )
                ),
                mainPanel(
                  tabsetPanel(
                    tabPanel("Genome Visualization",htmlOutput("txtgenomeviz"),plotOutput("genomeviz"),
                             htmlOutput("cgtable.title"),dataTableOutput("cgtable"), value=2)
                    , id = "conditionedPanels"
                  )
                )
)

server <- function(input, output) {
  
  # GET GENE COORDINATES/OTHER INFO
  geneinfoevent <- observeEvent(input$getgeneinfo,{
    
    # try to retrieve gene info
    if(input$genesym %in% xx){
      if(names(xx[xx==input$genesym]) %in% names(gx)){
        genesymchar <- as.character(input$genesym)
        grgene <<- gx[names(gx)==names(xx[xx==genesymchar])]
        output$geneinfo <- renderUI({
          str.con1 <- paste0("Gene Symbol: ",genesymchar)
          str.con2 <- paste(grgene,collapse="<br/>")
          HTML(paste(tags$span(style="color:green",str.con1)),"<br/>Gene Coordinates:<br/>",str.con2)
          
        })
      } else{
        output$geneinfo <- renderUI({
          str.error <- paste0("Error: Gene Info not found :(")
          HTML(paste(tags$span(style="color:red",str.error)))
        })  
      }
    } else{
      output$geneinfo <- renderUI({
        str.error <- paste0("Error: Gene Info not found :(")
        HTML(paste(tags$span(style="color:red",str.error)))
      })
    }
  })
  
  gettables <- observeEvent(input$loadgenecoordinates,{
    if(input$genesym2 %in% xx){
      if(names(xx[xx==input$genesym2]) %in% names(gx)){
        genesymchar2 <- as.character(input$genesym2)
        grgene2 <<- gx[names(gx)==names(xx[xx==genesymchar2])]
        
        startcoorobj <- start(grgene2)
        endcoorobj <- end(grgene2)
        chrgeneobj <- seqnames(grgene2)
        
        output$startcoor <- renderUI({
          numericInput("startcoorload", label = h5("Start Coordinate"), value = startcoorobj)
        })
        output$endcoor <- renderUI({
          numericInput("endcoorload", label = h5("End Coordinate"), value = endcoorobj)
        })
        output$chrgene <- renderUI({
          textInput("chrgeneload", label = h5("Gene Chromosome"), value = chrgeneobj)
        })
      } else{
        output$txtgenomeviz <- renderUI({
          str.error <- paste0("Error: Gene Info not found :(")
          HTML(paste(tags$span(style="color:red",str.error)))
        })
      }
    } else{
      output$txtgenomeviz <- renderUI({
        str.error <- paste0("Error: Gene Info not found :(")
        HTML(paste(tags$span(style="color:red",str.error)))
      })
    }
    
  })
  
  
  # GENERATE GENOME VISUALIZATION
  dataanalysis <- observeEvent(input$viewgenome,{
    
    if(length(input$startcoorload)>0 &
       length(input$endcoorload)>0 &
       length(input$chrgeneload)>0 &
       is.numeric(input$endcoorload) & 
       is.numeric(input$startcoorload) & 
       !input$endcoorload<input$startcoorload & 
       is.character(input$chrgeneload)){
      
      # assembling a dinucleotide df with all nucleotide pairs in selected sequence
      windowrange.start <- input$startcoorload
      windowrange.end <- input$endcoorload
      windowrange.chr <- input$chrgeneload
      
      gtrack <- GenomeAxisTrack()
      itrack <- IdeogramTrack(genome = "hg19", chromosome = windowrange.chr)
      
      tracklistusr <- list(gtrack,itrack); tracksizesusr <- c(1,1)
      #====================================
      # get sequence, make cg track/table
      #====================================
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
      #dndf[(dndf$start=="C" & dndf$end=="G")|
      #       (dndf$start=="G" & dndf$end=="C"),]$val <- 1
      #table(dndf$val)
      
      # make GRanges obj from CG dinucleotide df for gviz
      grcg <- GRanges(seqnames=windowrange.chr,
                      ranges=IRanges(start=dndf$chr.coor.start,
                                     end=dndf$chr.coor.start+1),mcols=dndf$val)
      dTrack.cg <- DataTrack(start=start(grcg),
                             end=end(grcg),
                             data=mcols(grcg)[,1],
                             chromosome=windowrange.chr,
                             strand=strand(grcg),
                             genome="hg19",
                             name="CGs",
                             type="gradient",
                             showColorBar=FALSE,
                             ncolor=2,
                             col.axis=NULL)
      tracklistusr <- c(tracklistusr,dTrack.cg); tracksizesusr <- c(tracksizesusr,1)
      
      # cg table to be returned in output
      dndf.cg <- dndf[dndf$val==1,c(1,2,4)]
      
      #==================================
      # Evaluate checkbox optional tracks
      #==================================
      # gene transcripts track
      if("1" %in% input$trackoptions){
        biomTrack <- BiomartGeneRegionTrack(genome = "hg19",
                                            chromosome = windowrange.chr, 
                                            start = windowrange.start-1000, 
                                            end = windowrange.end,
                                            name = "ENSEMBL",
                                            stacking="full")
        tracklistusr <- c(tracklistusr,biomTrack); tracksizesusr <- c(tracksizesusr,3)
      }
      
      # misc refgene tracks from ucsc
      if("2" %in% input$trackoptions){
        cpgIslands <- UcscTrack(genome = "hg19", chromosome = windowrange.chr,
                                track = "cpgIslandExt", from = windowrange.start, to = windowrange.end,
                                trackType = "AnnotationTrack", start = "chromStart",
                                end = "chromEnd", id = "name", shape = "box",
                                fill = "red", name = "CpG Islands")
        tracklistusr <- c(tracklistusr,cpgIslands); tracksizesusr <- c(tracksizesusr,1)
      }
      
      if("3" %in% input$trackoptions){
        gcContent <- UcscTrack(genome = "hg19", chromosome = windowrange.chr,
                               track = "GC Percent", table = "gc5Base", from = windowrange.start,
                               to = windowrange.end, trackType = "DataTrack", start = "start",
                               end = "end", data = "score", type = "hist",
                               windowSize = 1000, fill.histogram = "black",
                               col.histogram = "pink", ylim = c(30, 70), name = "GC Percent")
        tracklistusr <- c(tracklistusr,gcContent); tracksizesusr <- c(tracksizesusr,1)
      }
      
      #=============================
      # output to data analysis tab
      #=============================
      # print title for genome visualiZation
      output$txtgenomeviz <- renderUI({
        str.gviz <- paste0("Genome Visualization at: ",windowrange.chr,":",windowrange.start,"-",windowrange.end)
        HTML(paste(h4(str.gviz)))
      })
      
      # plot genome visualization
      output$genomeviz <- renderPlot(
        plotTracks(tracklistusr, 
                   from = windowrange.start, 
                   to = windowrange.end, 
                   cex = 0.8,
                   sizes=tracksizesusr)
      )
      
      # cg coordinates table title
      output$cgtable.title <- renderUI({
        str.gviz <- paste0("CG Coordinates Table at: ",windowrange.chr,":",windowrange.start,"-",windowrange.end)
        HTML(paste(h4(str.gviz)))
      })
      
      # cg coordinates data table
      output$cgtable <- renderDataTable(
        dndf.cg
      )
      
      
      # enable download for current selected table in csv format
      output$dncgtable_download.csv <- downloadHandler(
        filename = function() { paste0(input$cgtable,'.csv') },
        content = function(file) {
          write.csv(dndf.cg, file,row.names=FALSE)
        }
      )
      
    } else{
      output$txtgenomeviz <- renderUI({
        str.error <- paste0(h4("Error: Invalid chromosome and/or start, end coordinates :("))
        HTML(paste(str.error))
      })
      output$genomeviz <- renderPlot()
      output$cgtable.title <- renderUI()
      output$cgtable <- renderDataTable()
      
    }
  })
}



shinyApp(ui = ui, server = server)
