<title>"cgmappeR ReadMe" </title>

<!-- badges: start -->

[![Zotero]()
[![R build status]()

<!-- badges: end -->

# cgmappeR

An R/Shiny dashboard to map and visualize cytosine-guanine (CG/CpG) dinucleotides and microarray-targeted probes.

# Cite

To cite this tool, you may use the following: 

```
```

# Tutorial

## Setup

Clone the latest version of `cgmappeR` from your terminal with:

```
git clone https://github.com/metamaden/cgmappeR"
```

Now install all dependencies by calling the script `cgmapper/inst/r/install.R`:

```
R ./inst/install.R
```

This installs dependencies from CRAN and Bioconductor. These include `shiny`, `shinythemes`, `shinyWidgets`, `devtools`, `Gviz`, `BSgenome.Hsapiens.UCSC.hg19`, `org.Hs.eg.db`, and `TxDb.Hsapiens.UCSC.hg19.knownGene`.

## Run

<video width="550" controls>
    <source src="./inst/mp4/cgmappeR_run1_clip2.mp4" type="video/mp4">
    Your browser does not support the video tag.
</video>

Run the dashboard from shell with:

```
Rscript ./cgmappeR/inst/run.R
```

<video width="550" controls>
    <source src="./inst/mp4/cgmappeR_run1_clip2.mp4" type="video/mp4">
    Your browser does not support the video tag.
</video>

The `cgmappeR` UI is reactive to user input, and will start with no input. Specify inputs such as the coordinate ranges and your input data using the menu options to the left.

Iteratively map CG dinucleotides and Illumina CpG probe locations in genome ideograms. 

<img src="https://github.com/metamaden/cgmappeR/blob/master/readme_content/cgbrowseR_tp53exe.JPG" width="500">

## Query

<video width="550" controls>
    <source src="./inst/mp4/cgmappeR_query1_clip1.mp4" type="video/mp4">
    Your browser does not support the video tag.
</video>

Enter a valid gene symbol, load its coordinates, and modify coordinates for the ideogram window.

<img src="https://github.com/metamaden/cgmappeR/blob/master/readme_content/readme_instructions1.JPG" height="350">

Using dropdown menus, select tracks to visualize and fine tune image dimensions. Optionally add a custom cursor using custom coordinates. (Note: CpG tracks are only shown if CpGs overlap with view window)

<img src="https://github.com/metamaden/cgmappeR/blob/master/readme_content/readme_instructions2.JPG" height="350">

Click View Genome button to load the genome ideogram at the indicated coordinates.

<img src="https://github.com/metamaden/cgmappeR/blob/master/readme_content/readme_instructions3.JPG" width="300">

View resultant ideogram (first tab), CG dinucleotide table (second tab), CpG probe annotations table (third tab), and/or sequence in the selected window (fourth tab). 

## Dashboard navigation

The lefthand menu contains filter, query, import, and export options. 

The different tabs will switch the display between complementary information for a specified view window, including ideogram plot, dinucleotide coordinates table, probe annotations, and genome sequence.

View Ideogram Tab:

<img src="https://github.com/metamaden/cgmappeR/blob/master/readme_content/readme_instructions4.JPG" width="450">

View CG Dinucleotides Table Tab:

<img src="https://github.com/metamaden/cgmappeR/blob/master/readme_content/readme_instructions5.JPG" width="450">

View CpG Probes Table Tab:

<img src="https://github.com/metamaden/cgmappeR/blob/master/readme_content/readme_instructions6.JPG" width="450">

View Sequence Tab:

<img src="https://github.com/metamaden/cgmappeR/blob/master/readme_content/readme_instructions7.JPG" width="450">

## View management

Use the coordinates and other filters at the lefthand menu in order to modify the view window, such as shift upstream, downstream, expand the view window, or narrow the view window.

## Import

You may view your own probe data in the `cgmappeR` ideogram by selecting from the left menu.

## Export

You may save the probe annotations, dinucleotide coordinates, and sequence using the bottom-left menu options.

Download the image by right-clicking, and download the tables using the download buttons.

<img src="https://github.com/metamaden/cgmappeR/blob/master/readme_content/readme_instructions8.JPG" width="450">
