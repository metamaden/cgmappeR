# CGMappeR v.1.2.0
R shiny app to map and visualize genomic CG dinucleotides and CpG probes.

## Download/Install instructions
1. Install the latest versions of [R](https://cran.r-project.org/) and [R Studio](https://www.rstudio.com/products/rstudio/download/).

2. Install essential packages off CRAN:
`install.packages(c("shiny","shinythemes","shinyWidgets"))`

3. Install essential packages off Bioconductor:
```
source("https://bioconductor.org/biocLite.R");   
biocLite(c("Gviz", "BSgenome.Hsapiens.UCSC.hg19", "org.Hs.eg.db", "TxDb.Hsapiens.UCSC.hg19.knownGene","devtools"))
## At the interactive prompt, enter 'n'
```

4. Install and expand the latest version of the cgmappeR app locally, from GitHub.

5. Open R Studio and change current session to the main cgmappeR directory. Open the "app.R" script, and click "Run App" button.

## GUI 

Iteratively map CG dinucleotides and Illumina CpG probe locations in genome ideograms. 

<img src="https://github.com/metamaden/cgmappeR/blob/master/readme_content/cgbrowseR_tp53exe.JPG" width="500">

## Instructions
1. Enter a valid gene symbol, load its coordinates, and modify coordinates for the ideogram window.

<img src="https://github.com/metamaden/cgmappeR/blob/master/readme_content/readme_instructions1.JPG" height="350">

2. Using dropdown menus, select tracks to visualize and fine tune image dimensions. Optionally add a custom cursor using custom coordinates. (Note: CpG tracks are only shown if CpGs overlap with view window)

<img src="https://github.com/metamaden/cgmappeR/blob/master/readme_content/readme_instructions2.JPG" height="350">

3. Click View Genome button to load the genome ideogram at the indicated coordinates. (Note: this may take awhile)

Loading Gviz:

<img src="https://github.com/metamaden/cgmappeR/blob/master/readme_content/readme_instructions3.JPG" width="300">

4. View resultant ideogram (first tab), CG dinucleotide table (second tab), CpG probe annotations table (third tab), and/or sequence in the selected window (fourth tab). 

View Ideogram Tab:

<img src="https://github.com/metamaden/cgmappeR/blob/master/readme_content/readme_instructions4.JPG" width="450">

View CG Dinucleotides Table Tab:

<img src="https://github.com/metamaden/cgmappeR/blob/master/readme_content/readme_instructions5.JPG" width="450">

View CpG Probes Table Tab:

<img src="https://github.com/metamaden/cgmappeR/blob/master/readme_content/readme_instructions6.JPG" width="450">

View Sequence Tab:

<img src="https://github.com/metamaden/cgmappeR/blob/master/readme_content/readme_instructions7.JPG" width="450">

Download the image by right-clicking, and download the tables using the download buttons.

<img src="https://github.com/metamaden/cgmappeR/blob/master/readme_content/readme_instructions8.JPG" width="450">


## Citations

This is a shiny app written in R. It relies heavily on several Bioconductor packages, including Gviz, BSgenome.Hsapiens.UCSC.hg19, org.Hs.eg.db, TxDb.Hsapiens.UCSC.hg19.knownGene, and manifests accessible in minfi. This app was designed using the shiny, shinythemes, and shinyWidgets packages.

## Disclaimer

CGMappeR, including its code and generated results, is free to use for research purposes. It is offered with absolutely no warranty or guarantee, and it is the responsibility of the user to verify and/or validate any findings from using CGMappeR.

### Thanks for your interest in this project, and happy mapping!

#
## NEWS/TODO
### 1/20/18
Ideas for additional features: CRXCG probes (cross-reactive probes) or optionally enable subtraction of CRXCG CpGs from viewed CpGs? Enable visualization of Methylation (eg. Beta) values from custom specified matrix. Clean up interface: Remove meaningless axis labels, touch-up proportions and consider options for automating proportion adjustment. Enable viewing of ENCODE enhancer data from UCSC or consider options for downloading objects locally.

### 3/27/18
Massive v.1.2.0 update is live! New menu options, cleaner display (no meaningless y-axes), options to upload and visualize methylation data, faster TxDb gene mapping option, instructions for data upload, new default selections for tracks. 




