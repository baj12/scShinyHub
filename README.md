# Shiny app for exploration of single Cell data

This is a shiny / RStudio app for the exploration and analysis of single cell data. It is currently being developed based on 10X datat (one experiment).
The original version of this app is CellView (https://github.com/mohanbolisetty/CellView), but it was substantially modified. It just helped me get started.

## Trello board:

 https://trello.com/b/aZ4iLvAA
 
 Here, I orgnaize the development tasks of the application. Could probably also be done within Gitlab, but I am more used to working with Trello...
 
## Installation:

Download git repository from gitlab.

```R
requiredPackages= c("shiny", "shinyTree", "plotly", "shinythemes", "ggplot2", "DT", "pheatmap", "threejs", "sm", "RColorBrewer", "mclust", "reshape", "cellrangerRkit", "SCORPIUS", "knitr", "kableExtra", "shinyWidgets", "scater")
install.packages("rafalib")
rafalib::install_bioc(requiredPackages)
```

## Running the app

Simply open the file server.R in RStudio and clikck the Run button in RStudio.

## Description of repository

* contributions
  User contributions can be implemented here. Files like ui.R, server.R are being read by the ui.R and server.R files in the root directory. See "Develoment guide" for further information.
  
* aad0501_Table_S3.xlsx
  Excel file with cell type specific genes (from: Dissecting the multicellular ecosystem of metastatic melanoma by single-cell RNA-seq)
  http://science.sciencemag.org/content/suppl/2016/04/07/352.6282.189.DC1

* Examples
  directory containing example data coming from the CellView package.
  
* geneLists.R
  script to create geneLists.RData from proteinatlas.tsv(.zip), used for gene-classes in gene selection tab.
  
* privatePlotFunctions.R
  externalized functions for plotting.

* reactives.R
  Global reactive variables used in server.R
  
* report.Rmd(.html)
  global report template (needs some heavy recreation and plugin(ization))
  
* serverFunctions.R
  other externalized functions used in server.R, coming from the original CellView package.
  
* shinyTreeExample.R 
  test script to test shinyTree
  
* tabs.R
  ui elements used in ui.R that were not externalized (namely for the cell and gene selection, and the basic stats on the side panel).

* ui.R
  main ui creation for the shiny app. Handles import of ui.R from contributions
  
* reformtExamples.R
  script to reformat the examples from CellView to a structure that can be used by this app.

## Development guide

One of the challenges when dealing with scRNA seq experients and R is the vast number of data structures that is being used by the different packages. One task that the base version of this app is doing is to provide reactives that convert from the original input data to the required data structure. This has been done so far for the cellrange package (https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/rkit) and the scater (https://bioconductor.org/packages/release/bioc/html/scater.html) packages.


### components/plug-ins

Plugins are components under the contributions directory that follow some basic conventions. They will be loaded on run-time to the application by parsing the child direcotries. They abide to the following basic structure. There is also a Dummy component that holds some commented out examples.

Directory names:

`contributions/NAME/component/`

`NAME` should be replaced by the contributor's name and `component` by a somewhat meaningful name of the component.

Files that are loaded from these directories:

* outputs.R
store the output.R server side code here.

* reactives.R
define any new reactives here.
In order to be able to use the same functions/values in the report and the reactive we have to separate the reactive part from the calculation part. We there for create a function directly before the reactive that does the calculations, given the reactive variable that it uses as arguments. Here is a short example:

medianUMIfunc <- function(log2cpm){
    umiC = colSums(log2cpm, na.rm = TRUE)
    if(DEBUG)cat(file=stderr(), "medianUMI:done\n")
    return(median(t(umiC)))
  }
  
medianUMI <- reactive({
  if(DEBUG)cat(file=stderr(), "medianUMI\n")
  log2cpm = log2cpm()
  if(is.null(log2cpm) ){
    if(DEBUG)cat(file=stderr(), "medianUMI:NULL\n")
    return(0)
  }
  return(medianUMIfunc(log2cpm))
})

In the above example, medianUMI only uses reactives that are not available in the markdown and assign them to "real" variables and checks their conditions. We then return the actual value that is needed based only on non-reactives. This function (medainUMIfunc) can be used in the markdown to produce the desired result based solely on non-reactive values that are made available using params.

The function has to be in reactive.R such that we can import it in the R markdown.

* report.Rmd 
  see below

* ui.R
User interface code goes here.

Each of outputs.R or reactives.R can have a variable called:

`myHeavyCalculations`

that holds a list of names and functions(reactives) that can be executed by pressing the "Force Calculations" button.

E.g.:
```R
myHeavyCalculations = list(c("pcaName", "pca"),
                           c("kmClusteringName", "kmClustering"))
```
Here, pcaName and kmClusteringName are the names displayed on the console/gui element during execution and pca and kmClustering are the names of the reactives that are called.

See also below Modules for global variables / reactives that can be used in any component.

### Modules

#### conventions used here.

* Shiny modules can be considered as functions or self-contained, integratable shiny applications. They contain of a UI element and the server side logic. There are recommendations for the naming of shiny modules, which state that the UI function should end in UI, Input, or Output, and that the server side side function should be the base name. I find it more convenient to have the server side function end in "Server".

* access to global data:

The following variables are defined as reactives on the base level and can be access from anywhere without the need for installing a plugin.
Be aware that the standard parameters are set for project that I am currently working on and might not reflect your needs.

* `useGenes`
  genes selected using the user interface. 

* `useCells`
  cells selected using the user interface. 
  
* `gbm`
  data structure from cellranger that holds the original unmodified count data. Contains only the genes/cells that have been selected.
  
* `featureDataReact`
  feature data with gene name mappings. Contains only the genes/cells that have been selected.

* `log2cpm`
  log2 transformed data coming from the origianl CellView application. Contains only the genes/cells that have been selected.

* `tnse.data`
  tnse coordinates from run_tsne from the cellranger package.  Contains only the genes/cells that have been selected.

* `inputData`
  base data coming from Rds file that is loaded initially, no genes/cells removed.

* `medianENSG`
  median of the number of genes with at least one read of all samples
  
* `medianUMI`
  median expression of all samples
  
* `gbm_log`
  log 10 transformed data from cellranger package. (It uses normalize_barcode_sums_to_median.) Contains only the genes/cells that have been selected.

* `pca`
  results from run_pca from the cellranger package. 

* `kmClustering`
  results from run_kmeans_clustering from the cellranger package with 10 clusters

### tool tips

use tipify to surround the element you want to provide further information for.

### reports from contributions / plugins

* report.Rmd

The file report.Rmd is searched for under the contributions directory and all files are included in the html report via the child attribute

Please make sure that the code chunks are named and unique. Otherwise we might run into problems with duplicate chunk names. Unique names should begin with a unique ID, e.g. the naming that follows the directory structure used for the contributions.

all inputs are available for the report. When running with DEBUG=TRUE a report.RData file is generated on the desktop that can be used to develop the report as it contains all available variables. Inputs from plugins are also available.

