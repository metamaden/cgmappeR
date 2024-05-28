#!/usr/bin/env R

# Author: Sean Maden
# 
# Install dependencies for cgmappeR
#

# dependencies
install.packages(c("shiny","shinythemes","shinyWidgets","devtools"))
options(install.packages.compile.from.source = "always")
BiocManager::install(
  c("Gviz", "org.Hs.eg.db",
    "BSgenome.Hsapiens.UCSC.hg19", 
    "TxDb.Hsapiens.UCSC.hg19.knownGene"),
  force = TRUE)

# repo

