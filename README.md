# Shiny app for exploration of single Cell data

This is a shiny / RStudio app for the exploration and analysis of single cell data. It is currently being developed based on 10X datat (one experiment).
The original version of this app is CellView (https://github.com/mohanbolisetty/CellView), but it was substantially modified. It just helped me get started.

## Trello board:

 https://trello.com/b/aZ4iLvAA
 
 Here, I orgnaize the development tasks of the application. Could probably also be done within Gitlab, but I am more used to working with Trello...
 
## Installation:

```R
requiredPackages= c("shiny", "shinyTree", "plotly", "shinythemes", "ggplot2", "DT", "pheatmap", "threejs", "sm", "RColorBrewer", "mclust", "reshape", "cellrangerRkit", "SCORPIUS", "knitr", "kableExtra", "shinyWidgets", "scater")
install.packages("rafalib")
rafalib::install_bioc(requiredPackages)
```

## Running the app

Simply open server.R and clikck the Run button in RStudio.

## Description of repository

* contributions
  User contributions can be implemented here. Files like ui.R, server.R are being read by the ui.R and server.R files in the root directory. See "Develoment guide" for further information.
  
* aad0501_Table_S3.xlsx
  Excel file with cell type specific genes (from: Dissecting the multicellular ecosystem of metastatic melanoma by single-cell RNA-seq)
  http://science.sciencemag.org/content/suppl/2016/04/07/352.6282.189.DC1
  
* geneLists.R
  script to create geneLists.RData from proteinatlas.tsv(.zip), used for gene-classes in gene selection tab.
  
* privatePlotFunctions.R
  externalized functions for plotting.

* reactives.R
  Global reactive variables used in server.R
  
* report.Rmd(.html)
  global report template (needs some heavy recreation and plugin(ization))
  
* serverFunctions.R
  other externalized functions used in server.R
  
* shinyTreeExample.R 
  test script to test shinyTree
  
* tabs.R
  ui elements used in ui.R that were not externalized.


## Development guide

being developed

