# jscode <- '
# $(function() {
# var $els = $("[data-proxy-click]");
# $.each(
# $els,
# function(idx, el) {
# var $el = $(el);
# var $proxy = $("#" + $el.data("proxyClick"));
# $el.keydown(function (e) {
# if (e.keyCode == 13) {
# $proxy.click();
# }
# });
# }
# );
# });
# '
# 
source('modulesUI.R')
# this is where the general tabs are defined:
if(file.exists('defaultValues.R')){
  source('defaultValues.R')
}else{
  defaultValueSingleGene = ""
  defaultValueMultiGenes = ""
  defaultValueRegExGene = ""
}
# defaultValueSingleGene='malat1'
# defaultValueMultiGenes='MALAT1, B2M, TMSB4X, EEF1A1, FTH1'
# defaultValueRegExGene='' # tip: '^CD7$|^KIT$; genes with min expression

inputTab = tabItem(tabName = "input",
                   fluidRow(div(h3('CellView2'), align = 'center')),
                   br(),
                   fluidRow(div(
                     h5(
                       'This app is designed for exploratory data analysis of
                       processed RNA-Seq data of single cell experiments.'
                     ),
                     align = 'center'
                   )),
                   br(),
                   br(),
                   fluidRow(column(
                     5,
                     offset = 4,
                     fileInput(
                       'file1',
                       'Choose .Rds file to upload',
                       accept = c(
                         '.Rds',
                         'text/comma-separated-values',
                         'text/tab-separated-values',
                         'text/plain',
                         '.csv',
                         '.tsv'
                       )
                     )))
                   # br(),
                   # br(),
                   # fluidRow(
                   #   column(3, offset = 1,
                   #          textInput('defaultValueSingleGene', 'default value for single genes input', value = defaultValueSingleGene))
                   #   
                   # ),
                   # br(),
                   # fluidRow(
                   #   column(6, offset = 1,
                   #          textInput('defaultValueMultiGenes', 'default value for multi genes input', value = defaultValueMultiGenes))
                   #   
                   # )
                   
)


geneSelectionTab  = tabItem(tabName = "geneSelection",
                            fluidRow(div(h3('Gene selection'), align = 'center')),
                            br(),
                            fluidRow(div(
                              h4(
                                'Here we filter out genes'
                              ),
                              align = 'center'
                            )),
                            fluidRow(
                              column(3, offset = 1,
                                     textInput('selectIds', 'regular expression for selection of genes to be removed', value = '^MT-|^RP')),
                              column(5,
                                     h4('GeneList Selection'),
                                     shinyTree("geneListSelection", checkbox = TRUE)
                              ),
                              column(2,
                                     h4('Minoverall expression'),
                                     numericInput('minGenesGS', 'Min # of UMIs over all cells', 2, min=2, max = 1000000)
                              )
                            ),
                            fluidRow(
                              column(6, offset=1,
                                     textInput("genesKeep", "genes to keep"))
                            ),
                            br(),
                            fluidRow(
                              h3('Genes kept, with mean Expression, and number of cells expressing min 1', align = "center"),
                              br(),
                              h4('Selected genes'),
                              column(12, offset=0,
                                     textOutput("gsSelectedGenes", inline = FALSE)
                              ),
                              br(),
                              column(10, offset = 1,
                                     DT::dataTableOutput('selectedGenesTable')
                              )
                            ),
                            br(),
                            fluidRow(
                              h3('Genes removed, with mean Expression, and number of cells expressing min 1', align = "center"),
                              h4('Selected genes'),
                              br(),
                              textOutput("gsrmGenes", inline = FALSE)
                            ),br(),
                            fluidRow(
                              column(10, offset = 1,
                                     DT::dataTableOutput('removedGenesTable'))
                            )
)


cellSelectionTab  = tabItem(tabName = "cellSelection",
                            fluidRow(div(h3('Cell selection'), align = 'center')),
                            br(),
                            fluidRow(div(
                              h4(
                                'Here we filter out cells'
                              ),
                              align = 'center'
                            )),fluidRow(
                              column(6, offset=1,
                                     tipify(textInput('minExpGenes', 'List of genes with minimal expression', value = defaultValueRegExGene),
                                            title = "<h3>Cells must have one or more</h3> <ul><li>These cells must have at least one of those genes expressed</li> </ul> ", 
                                            options = list("width"='300px', 'placement'='right', 'max-width'='350px',
                                                           'data-container'="body", container="body")) # tool tip: '^CD7$|^KIT$
                              )
                            ),
                            fluidRow(
                              column(5, offset=1,
                                     numericInput('minGenes', 'Min # of UMIs', 2, min=2, max = 1000000)
                              ),
                              column(5,
                                     numericInput('maxGenes', 'Max # of UMIs', 1000000, min=10, max = 1000000)
                              )),br(),
                            fluidRow(
                              column(6, offset=1,
                                     textInput("cellSelectionComment", "Comment for selection of cells"))
                            ),
                            fluidRow(
                              column(6, offset=1,
                                     tipify(textInput("cellPatternRM", "cells to be filtered out by pattern"),
                                            title = "regular expression for cells to be removed (e.g. '-1' will remove all cells from sample 1")
                              )
                            ),
                            fluidRow(
                              column(6, offset=1,
                                     tipify(textInput("cellKeep", "cells to keep"),
                                            title="comma separated list of cells (with min expression) that should be kept")
                              )
                            ),fluidRow(
                              column(6, offset=1,
                                     tipify(textInput("cellKeepOnly", "cells to keep; remove others"),
                                            title="comma separated list of cells (with min expression) that should be kept and anything else removed")
                              )
                            ),
                            fluidRow(
                              column(10, offset = 1,
                                     tipify(textInput("cellsFiltersOut","Cells to be removed", width = '100%'),
                                            title = "comma separted list of cell names to be explicitly removed"))
                            ),br(),
                            fluidRow(column(
                              10,offset = 1,
                              tableSelectionUi("cellSelectionMod")
                            )
                            ),br()
                            
                            
)





