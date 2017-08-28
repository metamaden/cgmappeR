# gviz interactive shinyR browser

library(shiny)
library(shinythemes)
library(Gviz)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(org.Hs.eg.db)
library(BSgenome.Hsapiens.UCSC.hg19)
strack <- SequenceTrack(Hsapiens)

# get gene symbol mappings
x <- org.Hs.egSYMBOL
mapped_genes <- mappedkeys(x)
xx <- as.list(x[mapped_genes])
# get available gene coordinates
g = genes(TxDb.Hsapiens.UCSC.hg19.knownGene)
gx <- g[names(g) %in% names(xx)]

# GET GENE COORDINATES
grgene <- gx[names(gx)==names(xx[xx==genesym])]


ui <- fluidPage(theme=shinytheme("cerulean"),
                titlePanel(title=div(div(HTML(paste(tags$span(style="color:black",h1("CpG BrowseR"))))),
                                     header=h4("Enter Gene Info and Explore"))),
                sidebarPanel(
                  conditionalPanel(condition="input.conditionedPanels==1",
                                   helpText(h3("Gene Information")),
                                   helpText(h4("Enter your gene info:")),
                                   textInput("genesym", label = h5("Gene Symbol"), value = ""),
                                   actionButton("getgeneinfo", "Get Gene Information")
                  ),
                  conditionalPanel(condition="input.conditionedPanels==2",
                                   helpText(h3("Gene Visualization")),
                                   helpText("This tab displays a genome visualization"),
                                   textInput("genesym2", label = h5("Gene Symbol"), value = ""),
                                   actionButton("loadgenecoordinates", "Load Gene Coordinates"),
                                   uiOutput("startcoor"),
                                   uiOutput("endcoor"),
                                   uiOutput("chrgene"),
                                   actionButton("viewgenome", "View Genome")
                                   
                  )
                ),
                mainPanel(
                  tabsetPanel(
                    tabPanel("1. Gene Information",htmlOutput("geneinfo"), value=1), 
                    tabPanel("2. Genome Visualization",htmlOutput("txtgenomeviz"),
                             plotOutput("genomeviz"), value=2)
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
      
      gtrack <- GenomeAxisTrack()
      itrack <- IdeogramTrack(genome = unique(genome(grgene2)), chromosome = seqnames(grgene2))
      
      #gcContent <- UcscTrack(genome = unique(genome(grgene2)), chromosome = seqnames(grgene2),
      #                       track = "GC Percent", table = "gc5Base", from = input$startcoorload,
      #                       to = input$endcoorload, trackType = "DataTrack", start = "start",
      #                       end = "end", data = "score", type = "hist", window = -1,
      #                       windowSize = 1500, fill.histogram = "black",
      #                       col.histogram = "black", ylim = c(30, 70), name = "GC Percent")
      
      
      
      biomTrack <- BiomartGeneRegionTrack(genome = "hg19",
                                          chromosome = input$chrgeneload, 
                                          start = input$startcoorload, 
                                          end = input$endcoorload,
                                          name = "ENSEMBL")
      
      seqall <- getSeq(BSgenome.Hsapiens.UCSC.hg19,start=start(grgene2),end=end(grgene2),names=seqnames(grgene2))
      xseq.list <- unlist(strsplit(as.character(seqall),""))
      dndf <- data.frame(start=xseq.list[1:length(xseq.list)-1],
                         end=xseq.list[2:length(xseq.list)],
                         string.coor.start=seq(1,length(xseq.list)-1,1),
                         chr.coor.start=c(start(grgene2),start(grgene2)+seq(1,length(xseq.list)-2,1)))
      dndf.cg <- dndf[(dndf$start=="C" & dndf$end=="G")|
                        (dndf$start=="G" & dndf$end=="C"),]
      
      # output to data analysis tab
      output$txtgenomeviz <- renderUI({
        str.gviz <- paste0("Showing Genome Visualization for Gene: ",input$genesym2)
        HTML(paste(tags$span(style="color:blue",str.gviz)))
      })
      
      output$genomeviz <- renderPlot(
        plotTracks(list(gtrack,itrack,biomTrack,gcContent,strack), 
                   from = start(grgene), 
                   to = end(grgene), 
                   cex = 0.8)
      )
      
      
    } else{
      output$txtgenomeviz <- renderUI({
        str.error <- paste0("Error: Invalid chromosome and/or start, end coordinates :(")
        HTML(paste(tags$span(style="color:red",str.error)))
      })
    }
  }
  )
}

shinyApp(ui = ui, server = server)
