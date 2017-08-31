# CGMappeR
R shiny app to map and visualize genomic CG dinucleotides and CpG probes.

## GUI 

Iteratively map CG dinucleotides and Illumina CpG probe locations in genome ideograms. 

![Example: TP53 gene region CGs/CpGs](https://github.com/metamaden/cgmappeR/blob/master/readme_content/cgbrowseR_tp53exe.JPG)
<img src="https://github.com/metamaden/cgmappeR/blob/master/readme_content/cgbrowseR_tp53exe.JPG" width="2">
## Instructions
1. Enter a valid gene symbol, load its coordinates, and modify coordinates for the ideogram window.
![Example: Load Gene Coordinates](https://github.com/metamaden/cgmappeR/blob/master/readme_content/readme_instructions1.JPG)

2. Using dropdown menus, select tracks to visualize and fine tune image dimensions. Optionally add a custom cursor using custom coordinates. (Note: CpG tracks are only shown if CpGs overlap with view window)
![Example: Select Tracks](https://github.com/metamaden/cgmappeR/blob/master/readme_content/readme_instructions2.JPG)

3. Click View Genome button to load the genome ideogram at the indicated coordinates. (Note: this may take awhile)

Loading Gviz
![Example: Loading Gviz](https://github.com/metamaden/cgmappeR/blob/master/readme_content/readme_instructions3.JPG)

4. View resultant ideogram (first tab), CG dinucleotide table (second tab), CpG probe annotations table (third tab), and/or sequence in the selected window (fourth tab). 

View Ideogram Tab
![Example: View Ideogram Tab](https://github.com/metamaden/cgmappeR/blob/master/readme_content/readme_instructions4.JPG)

View CG Dinucleotides Table Tab
![Example: View CG Dinucleotides Table Tab](https://github.com/metamaden/cgmappeR/blob/master/readme_content/readme_instructions5.JPG)

View CpG Probes Table Tab
![Example: View CpG Probes Table Tab](https://github.com/metamaden/cgmappeR/blob/master/readme_content/readme_instructions6.JPG)

View Sequence Tab
![Example: View Sequence Tab](https://github.com/metamaden/cgmappeR/blob/master/readme_content/readme_instructions7.JPG)


Download CpG Probes Table
Download the image by right-clicking, and download the tables using the download buttons.
![Example: Download CpG Probes Table](https://github.com/metamaden/cgmappeR/blob/master/readme_content/readme_instructions8.JPG)


## Citations

This is a shiny app written in R. It relies heavily on several Bioconductor packages, including Gviz, BSgenome.Hsapiens.UCSC.hg19, org.Hs.eg.db, TxDb.Hsapiens.UCSC.hg19.knownGene, and manifests accessible in minfi. This app was designed using the shiny, shinythemes, and shinyWidgets packages.

## Disclaimer

CGMappeR, including its code and generated results, is free to use for research purposes. It is offered with absolutely no warranty or guarantee, and it is the responsibility of the user to verify and/or validate any findings from using CGMappeR.

### Thanks for your interest in this project, and happy mapping!

#
