#!/usr/bin/env R

# Author: Sean Maden
#
# Run the cgmappeR dashboard.
#
setwd("inst")
shiny::runApp('r/app.R')