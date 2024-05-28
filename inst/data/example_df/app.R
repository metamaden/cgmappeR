#====================================
# CgmappeR:
# gviz interactive shinyR browser
# author/maintainer: Sean Maden
#====================================

# shiny and cgmappeR dependencies:
library(shiny)
library(shinythemes)
library(Gviz)
library(shinyWidgets)

# need these 2 modules for TxDb track:
library(GenomicFeatures) 
library(TxDb.Hsapiens.UCSC.hg19.knownGene)

# need this for genome sequence info:
library(BSgenome.Hsapiens.UCSC.hg19)

# load objects and source scripts
load("grgenes_symbols.rda");load("grgenes.rda");load("epic450anno.rda")
source("cgbrowseR_functions.R")

ui <- fluidPage(theme=shinytheme("cerulean"),
                titlePanel(title=div(HTML(paste(h1("CGMappeR (hg19) v.1.2.0"))))),
                sidebarPanel(width=4,
                             helpText("Visualize CG dinucleotides, with tracks indicating coverage by CpG array probes, gene transcripts, and genome features."),
                             fluidRow(column(4,conditionalPanel(condition="input.conditionedPanels==2",
                                   dropdownButton(helpText(h5("1. Enter a gene to load coordinates into the UI, and optionally modify coordinates as desired.")),
                                                  helpText(h5("2. Upload user methylation data if desired. Note instructions for formatting, and parameter options in the UI menus.")),
                                                  helpText(h5("3. Select any additional annotation and genome feature tracks to display. You can also add a custom cursor track to highlight a specific range within the coordinate window being viewed.")),
                                                  helpText(h5("4. Click 'View Genome' button to: View resultant ideogram (first tab), get CG dinucleotide table (second tab), get CpG probe annotations table (third tab), and/or get sequence in the selected window (fourth tab). Download the image by right-clicking, and download the tables using the download buttons.")),
                                                  circle=FALSE,
                                                  label=h6("Instructions"),width="500px",size="xs"))),
                           column(4,dropdownButton(helpText(h4("Citations:")),
                                                 helpText(h6("This is a shiny app written in R. It relies heavily on several Bioconductor packages, including Gviz, BSgenome.Hsapiens.UCSC.hg19, org.Hs.eg.db, TxDb.Hsapiens.UCSC.hg19.knownGene, and manifests accessible in minfi. This app was designed using the shiny, shinythemes, and shinyWidgets packages.")),
                                                 helpText(h4("Disclaimer:")),
                                                 helpText(h6("CGMappeR, including its code and generated results, is free to use for research purposes. It is offered with absolutely no warranty or guarantee, and it is the responsibility of the user to verify and/or validate any findings from using CGMappeR.")),
                                                 helpText(h4("Thanks for your interest in this project, and happy mapping!")),
                                                 a(img(src="github-octocat.png",height=100,width=175,units="px"),alt="View CgmappeR on GitHub!",href="https://github.com/metamaden/cgmappeR"),
                                                 #uiOutput("example"), # octo-cat Image, link to GitHub page...
                                                 circle=FALSE,
                                                 label=h6("More Info"),
                                                 width="500px"))),
                                   HTML('<hr style="color: black;">'),
                                   helpText(h2("1. Load Gene Rgn./Enter Window Coord.")),
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
                                   helpText(h2("2. Methylation Data")),
                                   fluidRow(column(5,dropdownButton(helpText(h4("Methylation Data Display")),
                                                  helpText(h6("You can now upload your own methylation data, to display in cgmappeR with optional added tracks.")),
                                                  helpText(h4("Directions to Upload")),
                                                  helpText(h6("Select a CSV (comma-separated file) to upload. The file needs 3 columns with the following titles: (1) 'cpg'; (2) 'sampGroups'; (3) 'methyl'. See the following example:")),
                                                  img(src = "mdatExample.PNG"), # image needs to be in 'www' subdir
                                                  helpText(h6("Note that you can upload a single CSV (plot a single track) or multiple CSVs (plot two tracks). If one CSV is uploaded, multiple sample group levels will automatically be detected, if present, and indicated in different colors.")),
                                                  helpText(h6("Finally, you can optionally include a column called 'sampID' to delineate multiple samples/patients in the same treatment group. These will be processed automatically.")),
                                                  fileInput(inputId="file1",label="Choose Methyl. Data CSV 1", # fileInput is the methyl data file (#1)
                                                           accept = c("text/csv","text/comma-separated-values,text/plain",".csv")),
                                                  fileInput(inputId="file2",label="Choose Methyl. Data CSV 2", # fileInput is the methyl data file (#2)
                                                            accept = c("text/csv","text/comma-separated-values,text/plain",".csv")),
                                                  fluidRow(column(4,checkboxGroupInput("graph1para",
                                                                                     label="CSV1 parameters",
                                                                                     choices=c("Size Up-Scale" = 1,
                                                                                               "Auto-detect Sample Groups" = 2,
                                                                                               "Dotplot+Conf.Int." = 3),
                                                                                     selected=c("1","2","3"))), # checkbox to specify type(s) of graph parameters (CSV1)
                                                           column(4,checkboxGroupInput("graph2para",
                                                                                label="CSV2 parameters",
                                                                                choices=c("Size Up-Scale" = 1,
                                                                                          "Auto-detect Sample Groups" = 2,
                                                                                          "Dotplot+Conf.Int." = 3), # checkbox to specify type(s) of graph parameters (CSV2)
                                                                                selected=c("1","2","3")))),
                                                  circle=FALSE,
                                                  label=h6("Upload Data"),
                                                  width="500px")),
                                            column(7,uiOutput("datsummary"))),
                                   HTML('<hr style="color: black;">'),
                                   helpText(h2("3. Genome and Methyl. Array Tracks")),
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
                                                                              selected=c("1","2")),
                                                           circle=FALSE,
                                                           label=HTML(paste0(h6("CpG"),"\n",h6("Tracks"))),
                                                           width=2,
                                                           size="sm")),
                                            column(4,dropdownButton(checkboxGroupInput("trackoptions",
                                                                              label=h5(""),
                                                                              choices=c(
                                                                                "TxDb Annotation (v.fast)" = 5,
                                                                                "Ensembl Genes (fast)" = 1,
                                                                                "UCSC RefGenes (slow)" = 2,
                                                                                "CpG Islands (slow)" = 3,
                                                                                "GC Content (slow)" = 4),
                                                                              selected=c("5")),
                                                           circle=FALSE,
                                                           label=HTML(paste0(h6("Genome"),"\n",h6("Tracks"))),
                                                           width=2,
                                                           size="sm")),
                                            column(4,dropdownButton(checkboxInput("cursortrack", 
                                                                                  "Cursor Track", 
                                                                                  value = FALSE, 
                                                                                  width = NULL),
                                                                    numericInput("cursorstart", 
                                                                                 label = h5("Start Coor."), 
                                                                                 value = 0),
                                                                    numericInput("cursorend", 
                                                                                 label = h5("End Coor."), 
                                                                                 value = 0),
                                                                    circle=FALSE,
                                                                    label=HTML(paste0(h6("Cursor"))),
                                                                    width=2))),
                                   dropdownButton(numericInput("gvizwidth", 
                                                               label = h5("Image Width (pixels)"), 
                                                               value = 1000),
                                                  numericInput("gvizheight", 
                                                               label = h5("Image Height (pixels)"), 
                                                               value = 750),
                                                  circle=FALSE,
                                                  label=HTML(paste0(h6("Plot Dimensions"))),
                                                  width=2),
                                   HTML('<hr style="color: purple;">'),
                                   helpText(h2("4. Get Data/Results")),
                                   actionButton("viewgenome", "View Genome",icon=icon("angle-double-right")),
                                   downloadButton('dncgtable_download.csv', 'Download CG Table'),
                                   downloadButton('cpgprobetable_download.csv', 'Download CpG Probe Table'),
                                   HTML('<hr style="color: purple;">')
                  ),
                mainPanel(
                  tabsetPanel(
                    tabPanel("Ideogram",htmlOutput("tracks.selected"),htmlOutput("txtgenomeviz"),plotOutput("genomeviz"),value=2),
                    tabPanel("CG Table",htmlOutput("cgtable.title"),dataTableOutput("cgtable"), value=2),
                    tabPanel("CpG Probes Table",htmlOutput("cpgtable.title"),dataTableOutput("cpgtable"), value=2),
                    tabPanel("Sequence",
                             htmlOutput("txtgenomeviz.seq"),
                             wellPanel(id = "tPanel",style = "overflow-y:scroll; max-height: 600px",
                                       htmlOutput("txtseq")), 
                             value=2),
                    id = "conditionedPanels"
                  )
                )
)

server <- function(input, output) {
  
  #========================
  # WINDOW SIZE CALCULATOR
  #========================
  output$windowsize <- renderUI({
    windowsize.current <- input$endcoorload-input$startcoorload
    paste0("\nWindow Size (bp):\n",windowsize.current)
  })
  
  #===================
  # DATASET SUMMARY 1
  #===================
  output$datsummary <- renderUI({
    status.dialog <- list()
    if(!is.null(input$file1)|!is.null(input$file2)){
       if(!is.null(input$file1)){
         mi <- read.csv(input$file1$datapath)
         dimi <- paste0("CSV1 info:\nncol=",ncol(mi),",\nnrow=",nrow(mi),"\ngroups=",paste0(unique(mi$sampGroups),collapse=";"))
         status.dialog <- append(status.dialog,dimi); names(status.dialog)[length(status.dialog)] <- "file1_info" 
       }
      if(!is.null(input$file2)){
        mi <- read.csv(input$file2$datapath)
        dimi <- paste0("CSV2 info:\nncol=",ncol(mi),",\nnrow=",nrow(mi),"\ngroups=",paste0(unique(mi$sampGroups),collapse=";"))
        status.dialog <- append(status.dialog,dimi); names(status.dialog)[length(status.dialog)] <- "file2_info" 
      }
      
      paste0(status.dialog,collapse="\n")
    } else{
      return(paste0("No data selected."))
    }
    
  })
  
  # selected tracks object for plotting...excluding methyl. data table
  output$tracks.selected <- renderUI({
    names.selectedtracks <- c("ideogram","genome_coordinates","CGs")
    cursortrack <- c("Cursor")
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
                          "GC Content (slow)",
                          "TxDb Track (v. fast)")
    names.selectedtracks <- paste(c(names.selectedtracks,
                                    cursortrack[c(as.numeric(input$cursortrack))],
                                  arraytracknames[c(as.numeric(input$arraytracks))],
                                  genometracknames[c(as.numeric(input$trackoptions))]),
                                  collapse="; ")
    
    paste0("Tracks Selected:",names.selectedtracks)
  })
  
  # tab show genome coordinates
  output$txtgenomeviz.seq <- renderUI({
    paste0("Coordinates: ",input$chrgeneload,":",input$startcoorload,"-",input$endcoorload)
  })
  
  #=============================
  # TAB: GENOMIC SEQUENCE DISPLAY
  #=============================
  output$txtseq <- renderUI({
    if(input$startcoorload<=input$endcoorload){
      txtseqtxt <- as.character(getSeq(BSgenome.Hsapiens.UCSC.hg19,
                                       start=input$startcoorload,
                                       end=input$endcoorload,
                                       names=input$chrgeneload))
      paste(unlist(strsplit(txtseqtxt,"")),collapse='')
    }
  },outputArgs=list(inline=TRUE))
  
  #==================================
  # GET GENE COORDINATES/OTHER INFO
  #==================================
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
  
  #==================================
  # POPULATE GENE COORDINATE INTERFACE
  #==================================
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
  
  #================================
  # GENERATE GENOME VISUALIZATION
  #================================
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
      
      #===================================
      # ASSEMBLE THE DINUCLEOTIDE (CG) DF
      #===================================
      windowrange.start <- input$startcoorload
      windowrange.end <- input$endcoorload
      windowrange.chr <- input$chrgeneload
      
      gtrack <- GenomeAxisTrack()
      itrack <- IdeogramTrack(genome = "hg19", chromosome = windowrange.chr)
      
      tracklistusr <- list(itrack,gtrack); tracksizesusr <- c(1,1)
      #====================================
      # GET SEQUENCE AND MAKE CG TRACK TABLE
      #====================================
      
      dndf.all <- getAllDNseq(windowrange.start = input$startcoorload,
                          windowrange.end = input$endcoorload,
                          windowrange.chr = input$chrgeneload)
      dndf.cg.return <- getCGtable(windowrange.start = input$startcoorload,
                                   windowrange.end = input$endcoorload,
                                   windowrange.chr = input$chrgeneload,
                                   dndf=dndf.all,
                                   epic450anno = epic450anno)
      
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
                             showAxis=F)
      
      tracklistusr <- c(tracklistusr,dTrack.cg); tracksizesusr <- c(tracksizesusr,1)
      
      #==================================
      # Evaluate checkbox optional tracks
      #==================================
      
      # cursor track
      if(input$cursortrack){
        cursortrack <- DataTrack(start=input$cursorstart,
                                 end=input$cursorend,
                                 data=1,
                                 chromosome=windowrange.chr,
                                 genome="hg19",
                                 name="Cursor",
                                 type="gradient",
                                 showColorBar=FALSE,
                                 ncolor=2,
                                 gradient="black",
                                 showAxis=F)
        tracklistusr <- c(tracklistusr,cursortrack); tracksizesusr <- c(tracksizesusr,1)
      }
      
      # cpg array tracks
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
                               gradient=rep("red",nrow(ddf)),
                               showAxis=F)
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
                                 gradient=rep("purple",nrow(ddf)),
                                 showAxis=F)
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
                                   gradient=rep("darkolivegreen",nrow(ddf)),
                                   showAxis=F)
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
                                   gradient=rep("darkolivegreen2",nrow(ddf)),
                                   showAxis=F)
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
                                 gradient=rep("gold",nrow(ddf)),
                                 showAxis=F)
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
                                 gradient=rep("orange",nrow(ddf)),
                                 showAxis=F)
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
                                 gradient=rep("darkslategray1",nrow(ddf)),
                                 showAxis=F)
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
                                 gradient=rep("dodgerblue2",nrow(ddf)),
                                 showAxis=F)
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
                                 gradient=rep("lightcyan3",nrow(ddf)),
                                 showAxis=F)
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
      
      if("5" %in% input$trackoptions){
      
        txTr <- GeneRegionTrack(TxDb.Hsapiens.UCSC.hg19.knownGene, 
                                genome = "hg19",
                                chromosome = windowrange.chr, 
                                start = windowrange.start-1000, 
                                end = windowrange.end,
                                name="TxDb Gene")
        # note: use slightly up-scaled ratio for feature-rich regions..
        tracklistusr <- c(tracklistusr,txTr); tracksizesusr <- c(tracksizesusr,5)
      }
      
      
      if("1" %in% input$trackoptions){
        biomTrack <- BiomartGeneRegionTrack(genome = "hg19",
                                            chromosome = windowrange.chr, 
                                            start = windowrange.start-1000, 
                                            end = windowrange.end,
                                            name = "ENSEMBL",
                                            stacking="full",showAxis=F)
        tracklistusr <- c(tracklistusr,biomTrack); tracksizesusr <- c(tracksizesusr,3)
      }
      if("2" %in% input$trackoptions){
        refGenes <- UcscTrack(genome = "hg19", chromosome = windowrange.chr,
                              track = "xenoRefGene", from = windowrange.start, to = windowrange.end,
                              trackType = "GeneRegionTrack", rstarts = "exonStarts",
                              rends = "exonEnds", gene = "name", symbol = "name2",
                              transcript = "name", strand = "strand", fill = "#8282d2",
                              stacking = "full", name = "RefSeq",showAxis=F)
        tracklistusr <- c(tracklistusr,refGenes); tracksizesusr <- c(tracksizesusr,3)
      }
      
      # misc refgene tracks from ucsc
      if("3" %in% input$trackoptions){
        cpgIslands <- UcscTrack(genome = "hg19", chromosome = windowrange.chr,
                                track = "cpgIslandExt", from = windowrange.start, to = windowrange.end,
                                trackType = "AnnotationTrack", start = "chromStart",
                                end = "chromEnd", id = "name", shape = "box",
                                fill = "red", name = "CpG Islands",showAxis=F)
        tracklistusr <- c(tracklistusr,cpgIslands); tracksizesusr <- c(tracksizesusr,1)
      }
      
      if("4" %in% input$trackoptions){
        gcContent <- UcscTrack(genome = "hg19", chromosome = windowrange.chr,
                               track = "GC Percent", table = "gc5Base", from = windowrange.start,
                               to = windowrange.end, trackType = "DataTrack", start = "start",
                               end = "end", data = "score", type = "hist",
                               windowSize = 1000, fill.histogram = "black",
                               col.histogram = "pink", ylim = c(30, 70), name = "GC Percent",
                               showAxis=F)
        tracklistusr <- c(tracklistusr,gcContent); tracksizesusr <- c(tracksizesusr,1)
      }
      
      #=============================
      # OUTPUT TO DATA ANALYSIS TAB
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
        
        # DATA SET 1: for methyl. data csv file, store as data frame
        methyldatainput1 <-reactive({ 
          
          if (is.null(input$file1))
            return(NULL)                
          
          data1<-read.csv(input$file1$datapath)
          data1
        })
        
        # DATA SET 2: for methyl. data csv file, store as data frame
        methyldatainput2 <-reactive({ 
          
          if (is.null(input$file2))
            return(NULL)                
          
          data2<-read.csv(input$file2$datapath)
          data2
        })
        
        
        #======================
        # UPLOADED DATA TRACKS
        #======================
        if (!is.null(methyldatainput1())|!is.null(methyldatainput2())){
          
          if(!is.null(methyldatainput1())){
            methyldata1.track <- makeDtrackInfo(dfi=methyldatainput1(),dtrackname = "usr input 1",plottype=c("a","p"))
            tracklistusr <- c(tracklistusr,methyldata1.track); tracksizesusr <- c(tracksizesusr,3)
          }
          if(!is.null(methyldatainput2())){
            methyldata2.track <- makeDtrackInfo(dfi=methyldatainput2(),dtrackname = "usr input 1",plottype=c("a","p"))
            tracklistusr <- c(tracklistusr,methyldata2.track); tracksizesusr <- c(tracksizesusr,3)
          }
          
          # plot with data track(s)
          #plotTracks(tracklistusr, 
          #           from = windowrange.start, 
          #           to = windowrange.end, 
          #           cex = 0.8,
          #           sizes=tracksizesusr,
          #           groups=unique(methyldatainput1()$sampGroups),
          #           legend=T)
          plotTracks(tracklistusr, 
                     from = windowrange.start, 
                     to = windowrange.end, 
                     cex = 0.8,
                     sizes=tracksizesusr)
        
          } else{
            # plot without data track(s)
            plotTracks(tracklistusr, 
                     from = windowrange.start, 
                     to = windowrange.end, 
                     cex = 0.8,
                     sizes=tracksizesusr)
        }
      },
      height=input$gvizheight,width=input$gvizwidth
      )
      
      #==================================
      # CG COORDINATES TABLE AND TITLE
      #==================================
      output$cgtable.title <- renderUI({
        str.gviz <- paste0("CG Coordinates Table at: ",windowrange.chr,":",windowrange.start,"-",windowrange.end)
        HTML(paste(h4(str.gviz)))
      })
      output$cgtable <- renderDataTable({
        dndf.cg.return
      })
      
      # cpg probe table and title
      output$cpgtable.title <- renderUI({
        str.gviz <- paste0("CpG Probe Annotations at: ",windowrange.chr,":",windowrange.start,"-",windowrange.end)
        HTML(paste(h4(str.gviz)))
      })
      output$cpgtable <- renderDataTable({
        epic450anno[dndf.cg.return[!is.na(dndf.cg.return$cpg.id),]$cpg.id,]
      })
      
      
      # enable download for current selected CG dn. table in csv format
      output$dncgtable_download.csv <- downloadHandler(
        filename = function() { paste0(input$cgtable,'.csv') },
        content = function(file) {
          write.csv(dndf.cg.return, file,row.names=FALSE)
        }
      )
      # enable download for current selected CpG probe table in csv format
      output$cpgprobetable_download.csv <- downloadHandler(
        filename = function() { paste0(input$cpgtable,'.csv') },
        content = function(file) {
          write.csv(epic450anno[dndf.cg.return[!is.na(dndf.cg.return$cpg.id),]$cpg.id,], file,row.names=FALSE)
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
