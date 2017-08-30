# gviz interactive shinyR browser

library(shiny)
library(shinythemes)
library(Gviz)
library(shinyWidgets)
#library(TxDb.Hsapiens.UCSC.hg19.knownGene)
#library(org.Hs.eg.db)
library(BSgenome.Hsapiens.UCSC.hg19)
load("grgenes_symbols.rda");load("grgenes.rda");load("epic450anno.rda")
source("cgbrowseR_functions.R")

ui <- fluidPage(theme=shinytheme("cerulean"),
                titlePanel(title=div(HTML(paste(h1("CGMappeR"))))),
                sidebarPanel(width=4,
                  conditionalPanel(condition="input.conditionedPanels==2",
                                   helpText("Visualize CG dinucleotides, with tracks indicating coverage by CpG array probes, gene transcripts, and genome features."),
                                   dropdownButton(helpText(h5("1. Enter a valid gene symbol, load its coordinates, and modify coordinates for the ideogram window.")),
                                                  helpText(h5("2. Using dropdown menus, select tracks to visualize and fine tune image dimensions.")),
                                                  helpText(h5("3. Click View Genome button to load the genome ideogram at the indicated coordinates. (Note: this may take awhile)")),
                                                  helpText(h5("4. View resultant ideogram (first tab) and CG dinucleotide table (second tab). Download the image by right-clicking, and download the table using the download button.")),
                                                  circle=FALSE,
                                                  label=h6("Instructions"),width="500px",size="xs"),
                                   HTML('<hr style="color: black;">'),
                                   textInput("genesym2", label = h5("Gene Symbol"), value = ""),
                                   actionButton("loadgenecoordinates", "Load Gene Coordinates"),
                                   HTML('<hr style="color: black;">'),
                                   fluidRow(column(5,dropdownButton(uiOutput("startcoor"),
                                                                    uiOutput("endcoor"),
                                                                    uiOutput("chrgene"),
                                                                    circle=FALSE,
                                                                    label=HTML(paste0(h6("Coordinates"))),
                                                                    width=20,
                                                                    size="sm")),
                                            column(7,uiOutput("windowsize"))),
                                   HTML('<hr style="color: black;">'),
                                   fluidRow(column(3,dropdownButton(checkboxGroupInput("arraytracks",
                                                                              label="",
                                                                              choices=c("Array CpGs (EPIC)" = 1,
                                                                                        "Array CpGs (HM450)" = 2,
                                                                                        "Promoter CpGs (TSS1500)" = 3,
                                                                                        "Promoter CpGs (TSS200)" = 4,
                                                                                        "Body CpGs (Body)" = 5,
                                                                                        "Body CpGs (1stExon)" = 6,
                                                                                        "Promoter CpGs (5'UTR)" = 7,
                                                                                        "NTR CpGs (3'UTR)" = 8,
                                                                                        "Body CpGs (ExonBnd)" = 9),
                                                                              selected=c("1","2","3","4","5","6","7","8","9")),
                                                           circle=FALSE,
                                                           label=HTML(paste0(h6("CpG"),"\n",h6("Tracks"))),
                                                           width=2,
                                                           size="sm")),
                                            column(4,dropdownButton(checkboxGroupInput("trackoptions",
                                                                              label=h5(""),
                                                                              choices=c(
                                                                                "Ensembl Genes (fast)" = 1,
                                                                                "UCSC RefGenes (slow)" = 2,
                                                                                "CpG Islands (slow)" = 3,
                                                                                "GC Content (slow)" = 4),
                                                                              selected=c("1")),
                                                           circle=FALSE,
                                                           label=HTML(paste0(h6("Genome"),"\n",h6("Tracks"))),
                                                           width=2,
                                                           size="sm")),
                                            column(4,dropdownButton(numericInput("gvizwidth", 
                                                                                 label = h5("Image Width (pixels)"), 
                                                                                 value = 1000),
                                                                    numericInput("gvizheight", 
                                                                                 label = h5("Image Height (pixels)"), 
                                                                                 value = 750),
                                                                    circle=FALSE,
                                                                    label=HTML(paste0(h6("Plot"),"\n",h6("Dimensions"))),
                                                                    width=2))),
                                   HTML('<hr style="color: purple;">'),
                                   actionButton("viewgenome", "View Genome"),
                                   downloadButton('dncgtable_download.csv', 'Download CG Table'),
                                   HTML('<hr style="color: purple;">'),
                                   dropdownButton(helpText(h4("Citations:")),
                                                  helpText(h6("This is a shiny app written in R. It relies heavily on several Bioconductor packages, including Gviz, BSgenome.Hsapiens.UCSC.hg19, org.Hs.eg.db, TxDb.Hsapiens.UCSC.hg19.knownGene, and manifests accessible in minfi. This app was designed using the shiny, shinythemes, and shinyWidgets packages.")),
                                                  helpText(h4("Disclaimer:")),
                                                  helpText(h6("CGMappeR, including its code and generated results, is free to use for research purposes. It is offered with absolutely no warranty or guarantee, and it is the responsibility of the user to verify and/or validate any findings from using CGMappeR.")),
                                                  helpText(h4("Thanks for your interest in this project, and happy mapping!")),
                                                  circle=FALSE,
                                                  label=h6("More Info"),
                                                  width="500px")
                  )
                ),
                mainPanel(
                  tabsetPanel(
                    tabPanel("Genome Visualization",htmlOutput("tracks.selected"),htmlOutput("txtgenomeviz"),plotOutput("genomeviz"),value=2),
                    tabPanel("Sequence CG Table",htmlOutput("cgtable.title"),dataTableOutput("cgtable"), value=2),
                    id = "conditionedPanels"
                  )
                )
)

server <- function(input, output) {
  output$windowsize <- renderUI({
    windowsize.current <- input$endcoorload-input$startcoorload
    paste0("\nWindow Size (bp):\n",windowsize.current)
  })
  
  output$tracks.selected <- renderUI({
    names.selectedtracks <- c("ideogram","genome_coordinates","CGs")
    arraytracknames <- c("Array CpGs (EPIC)",
                         "Array CpGs (HM450)",
                         "Promoter CpGs (TSS1500)",
                         "Promoter CpGs (TSS200)",
                         "Body CpGs (Body)",
                         "Body CpGs (1stExon)",
                         "Promoter CpGs (5'UTR)",
                         "NTR CpGs (3'UTR)",
                         "Body CpGs (ExonBnd)")
    genometracknames <- c("Ensembl Genes (fast)",
                          "UCSC RefGenes (slow)",
                          "CpG Islands (slow)",
                          "GC Content (slow)")
    names.selectedtracks <- paste(c(names.selectedtracks,
                                  arraytracknames[c(as.numeric(input$arraytracks))],
                                  genometracknames[c(as.numeric(input$trackoptions))]),
                                  collapse="; ")
    
    paste0("Tracks Selected:\n",names.selectedtracks)
  })
  
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
          numericInput("startcoorload", label = h5("Start"), value = startcoorobj)
        })
        output$endcoor <- renderUI({
          numericInput("endcoorload", label = h5("End"), value = endcoorobj)
        })
        output$chrgene <- renderUI({
          textInput("chrgeneload", label = h5("Chrom."), value = chrgeneobj)
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
    
    withProgress(message = 'Calculating CG tracks',
                 detail = 'This may take a while...', value = 0, {
                   for (i in 1:15) {
                     incProgress(1/15)
                     Sys.sleep(0.25)
                   }
                 })
    
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
      
      tracklistusr <- list(itrack,gtrack); tracksizesusr <- c(1,1)
      #====================================
      # get sequence, make cg track/table
      #====================================
      dndf.all <- getAllDNseq(windowrange.start = input$startcoorload,
                          windowrange.end = input$endcoorload,
                          windowrange.chr = input$chrgeneload)
      dndf.cg.return <- getCGtable(windowrange.start = input$startcoorload,
                                   windowrange.end = input$endcoorload,
                                   windowrange.chr = input$chrgeneload,
                                   dndf=dndf.all)
      
      # make GRanges obj from CG dinucleotide df for gviz
      grcg <- GRanges(seqnames=windowrange.chr,
                      ranges=IRanges(start=dndf.all$chr.coor.start,
                                     end=dndf.all$chr.coor.start+1),
                      mcols=dndf.all$val)
      
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
      
      #==================================
      # Evaluate checkbox optional tracks
      #==================================
      
      if(length(intersect(c("1","2","3","4","5"),input$arraytracks))>0){
        dndf.array <- dndf.cg.return[!is.na(dndf.cg.return$array),]
        epic450dndf <- epic450anno[dndf.array[!is.na(dndf.array$cpg.id),]$cpg.id,]
      }
      
      # Array tracks
      if("1" %in% input$arraytracks){
        ddf <- epic450dndf[grep("epic",epic450dndf$array),]
        if(nrow(ddf)>0){
        dTrack.cpg1 <- DataTrack(start=ddf$pos,
                                 end=ddf$pos+1,
                               data=rep(1,nrow(ddf)),
                               chromosome=windowrange.chr,
                               genome="hg19",
                               name="EPIC CpGs",
                               type="gradient",
                               showColorBar=FALSE,
                               ncolor=2,
                               col.axis=NULL,
                               gradient=rep("red",nrow(ddf)))
        tracklistusr <- c(tracklistusr,dTrack.cpg1); tracksizesusr <- c(tracksizesusr,1)
        }
      }
      if("2" %in% input$arraytracks){
        ddf <- epic450dndf[grep("hm450",epic450dndf$array),]
        if(nrow(ddf)>0){
        dTrack.cpg2 <- DataTrack(start=ddf$pos,
                                 end=ddf$pos+1,
                                 data=rep(1,nrow(ddf)),
                                 chromosome=windowrange.chr,
                                 genome="hg19",
                                 name="HM450 CpGs",
                                 type="gradient",
                                 showColorBar=FALSE,
                                 ncolor=2,
                                 col.axis=NULL,
                                 gradient=rep("purple",nrow(ddf)))
        tracklistusr <- c(tracklistusr,dTrack.cpg2); tracksizesusr <- c(tracksizesusr,1)
        }
      }
      if("3" %in% input$arraytracks){
        ddf <- epic450dndf[grep("TSS1500",epic450dndf$UCSC_RefGene_Group),]
        if(nrow(ddf)>0){
          dTrack.cpg3 <- DataTrack(start=ddf$pos,
                                   end=ddf$pos+1,
                                   data=rep(1,nrow(ddf)),
                                   chromosome=windowrange.chr,
                                   genome="hg19",
                                   name="TSS1500 CpGs",
                                   type="gradient",
                                   showColorBar=FALSE,
                                   ncolor=2,
                                   col.axis=NULL,
                                   gradient=rep("darkolivegreen",nrow(ddf)))
          tracklistusr <- c(tracklistusr,dTrack.cpg3); tracksizesusr <- c(tracksizesusr,1)
        }
        
      }
      if("4" %in% input$arraytracks){
        ddf <- epic450dndf[grep("TSS200",epic450dndf$UCSC_RefGene_Group),]
        if(nrow(ddf)>0){
          dTrack.cpg4 <- DataTrack(start=ddf$pos,
                                   end=ddf$pos+1,
                                   data=rep(1,nrow(ddf)),
                                   chromosome=windowrange.chr,
                                   genome="hg19",
                                   name="TSS200 CpGs",
                                   type="gradient",
                                   showColorBar=FALSE,
                                   ncolor=2,
                                   col.axis=NULL,
                                   gradient=rep("darkolivegreen2",nrow(ddf)))
        tracklistusr <- c(tracklistusr,dTrack.cpg4); tracksizesusr <- c(tracksizesusr,1)
      }
      }
      if("5" %in% input$arraytracks){
        ddf <- epic450dndf[grep("Body",epic450dndf$UCSC_RefGene_Group),]
        if(nrow(ddf)>0){
        dTrack.cpg5 <- DataTrack(start=ddf$pos,
                                 end=ddf$pos+1,
                                 data=rep(1,nrow(ddf)),
                                 chromosome=windowrange.chr,
                                 genome="hg19",
                                 name="Body CpGs",
                                 type="gradient",
                                 showColorBar=FALSE,
                                 ncolor=2,
                                 col.axis=NULL,
                                 gradient=rep("gold",nrow(ddf)))
        tracklistusr <- c(tracklistusr,dTrack.cpg5); tracksizesusr <- c(tracksizesusr,1)
        }
      }
      if("6" %in% input$arraytracks){
        ddf <- epic450dndf[grep("1stExon",epic450dndf$UCSC_RefGene_Group),]
        if(nrow(ddf)>0){
        dTrack.cpg6 <- DataTrack(start=ddf$pos,
                                 end=ddf$pos+1,
                                 data=rep(1,nrow(ddf)),
                                 chromosome=windowrange.chr,
                                 genome="hg19",
                                 name="1stExon CpGs",
                                 type="gradient",
                                 showColorBar=FALSE,
                                 ncolor=2,
                                 col.axis=NULL,
                                 gradient=rep("orange",nrow(ddf)))
        tracklistusr <- c(tracklistusr,dTrack.cpg6); tracksizesusr <- c(tracksizesusr,1)
        }
      }
      if("7" %in% input$arraytracks){
        ddf <- epic450dndf[grep("5'UTR",epic450dndf$UCSC_RefGene_Group),]
        if(nrow(ddf)>0){
        dTrack.cpg7 <- DataTrack(start=ddf$pos,
                                 end=ddf$pos+1,
                                 data=rep(1,nrow(ddf)),
                                 chromosome=windowrange.chr,
                                 genome="hg19",
                                 name="5'UTR CpGs",
                                 type="gradient",
                                 showColorBar=FALSE,
                                 ncolor=2,
                                 col.axis=NULL,
                                 gradient=rep("darkslategray1",nrow(ddf)))
        tracklistusr <- c(tracklistusr,dTrack.cpg7); tracksizesusr <- c(tracksizesusr,1)
        }
      }
      if("8" %in% input$arraytracks){
        ddf <- epic450dndf[grep("3'UTR",epic450dndf$UCSC_RefGene_Group),]
        if(nrow(ddf)>0){
        dTrack.cpg8 <- DataTrack(start=ddf$pos,
                                 end=ddf$pos+1,
                                 data=rep(1,nrow(ddf)),
                                 chromosome=windowrange.chr,
                                 genome="hg19",
                                 name="3'UTR CpGs",
                                 type="gradient",
                                 showColorBar=FALSE,
                                 ncolor=2,
                                 col.axis=NULL,
                                 gradient=rep("dodgerblue2",nrow(ddf)))
        tracklistusr <- c(tracklistusr,dTrack.cpg8); tracksizesusr <- c(tracksizesusr,1)
        }
      }
      if("9" %in% input$arraytracks){
        ddf <- epic450dndf[grep("ExonBnd",epic450dndf$UCSC_RefGene_Group),]
        if(nrow(ddf)>0){
        dTrack.cpg9 <- DataTrack(start=ddf$pos,
                                 end=ddf$pos+1,
                                 data=rep(1,nrow(ddf)),
                                 chromosome=windowrange.chr,
                                 genome="hg19",
                                 name="ExonBnd CpGs",
                                 type="gradient",
                                 showColorBar=FALSE,
                                 ncolor=2,
                                 col.axis=NULL,
                                 gradient=rep("lightcyan3",nrow(ddf)))
        tracklistusr <- c(tracklistusr,dTrack.cpg9); tracksizesusr <- c(tracksizesusr,1)
        }
      }
      
      
      # gene transcripts track
      withProgress(message = 'Retrieving genome tracks',
                   detail = 'This may take a while...', value = 0, {
                     for (i in 1:15) {
                       incProgress(1/15)
                       Sys.sleep(0.25)
                     }
                   })
      
      if("1" %in% input$trackoptions){
        biomTrack <- BiomartGeneRegionTrack(genome = "hg19",
                                            chromosome = windowrange.chr, 
                                            start = windowrange.start-1000, 
                                            end = windowrange.end,
                                            name = "ENSEMBL",
                                            stacking="full")
        tracklistusr <- c(tracklistusr,biomTrack); tracksizesusr <- c(tracksizesusr,3)
      }
      if("2" %in% input$trackoptions){
        refGenes <- UcscTrack(genome = "hg19", chromosome = windowrange.chr,
                              track = "xenoRefGene", from = windowrange.start, to = windowrange.end,
                              trackType = "GeneRegionTrack", rstarts = "exonStarts",
                              rends = "exonEnds", gene = "name", symbol = "name2",
                              transcript = "name", strand = "strand", fill = "#8282d2",
                              stacking = "full", name = "RefSeq")
        tracklistusr <- c(tracklistusr,refGenes); tracksizesusr <- c(tracksizesusr,3)
      }
      
      # misc refgene tracks from ucsc
      if("3" %in% input$trackoptions){
        cpgIslands <- UcscTrack(genome = "hg19", chromosome = windowrange.chr,
                                track = "cpgIslandExt", from = windowrange.start, to = windowrange.end,
                                trackType = "AnnotationTrack", start = "chromStart",
                                end = "chromEnd", id = "name", shape = "box",
                                fill = "red", name = "CpG Islands")
        tracklistusr <- c(tracklistusr,cpgIslands); tracksizesusr <- c(tracksizesusr,1)
      }
      
      if("4" %in% input$trackoptions){
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
        str.gviz <- paste0("Coordinates: ",windowrange.chr,":",windowrange.start,"-",windowrange.end)
        HTML(paste(h4(str.gviz)))
      })
      
      # plot genome visualization
      output$genomeviz <- renderPlot({
        withProgress(message = 'Plotting genome tracks',
                     detail = 'This may take a while...', value = 0, {
                       for (i in 1:15) {
                         incProgress(1/15)
                         Sys.sleep(0.25)
                       }
                     })
        plotTracks(tracklistusr, 
                   from = windowrange.start, 
                   to = windowrange.end, 
                   cex = 0.8,
                   sizes=tracksizesusr)
      },
      height=input$gvizheight,width=input$gvizwidth
      )
      
      # cg coordinates table title
      output$cgtable.title <- renderUI({
        str.gviz <- paste0("CG Coordinates Table at: ",windowrange.chr,":",windowrange.start,"-",windowrange.end)
        HTML(paste(h4(str.gviz)))
      })
      
      # cg coordinates data table
      output$cgtable <- renderDataTable({
        dndf.cg.return
      }
      )
      
      
      # enable download for current selected table in csv format
      output$dncgtable_download.csv <- downloadHandler(
        filename = function() { paste0(input$cgtable,'.csv') },
        content = function(file) {
          write.csv(dndf.cg.return, file,row.names=FALSE)
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
