# CGMappeR
R shiny app to map and visualize genomic CG dinucleotides and CpG probes.

## GUI 

Iteratively map CG dinucleotides and Illumina CpG probe locations in genome ideograms. 

![Example: TP53 gene region CGs/CpGs](cgbrowseR_tp53eg.JPG)

## Instructions
1. Enter a valid gene symbol, load its coordinates, and modify coordinates for the ideogram window.

2. Using dropdown menus, select tracks to visualize and fine tune image dimensions. Optionally add a custom cursor using custom coordinates. (Note: CpG tracks are only shown if CpGs overlap with view window)

3. Click View Genome button to load the genome ideogram at the indicated coordinates. (Note: this may take awhile)

4. View resultant ideogram (first tab), CG dinucleotide table (second tab), CpG probe annotations table (third tab), and/or sequence in the selected window (fourth tab). Download the image by right-clicking, and download the tables using the download buttons.

## Citations

This is a shiny app written in R. It relies heavily on several Bioconductor packages, including Gviz, BSgenome.Hsapiens.UCSC.hg19, org.Hs.eg.db, TxDb.Hsapiens.UCSC.hg19.knownGene, and manifests accessible in minfi. This app was designed using the shiny, shinythemes, and shinyWidgets packages.

## Disclaimer

CGMappeR, including its code and generated results, is free to use for research purposes. It is offered with absolutely no warranty or guarantee, and it is the responsibility of the user to verify and/or validate any findings from using CGMappeR.

### Thanks for your interest in this project, and happy mapping!

#
