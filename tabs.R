jscode <- '
$(function() {
var $els = $("[data-proxy-click]");
$.each(
$els,
function(idx, el) {
var $el = $(el);
var $proxy = $("#" + $el.data("proxyClick"));
$el.keydown(function (e) {
if (e.keyCode == 13) {
$proxy.click();
}
});
}
);
});
'



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
                     )
                   ))
                   
)


geneSelectionTab  = tabItem(tabName = "geneSelection",
                            fluidRow(div(h3('Gene selection'), align = 'center')),
                            br(),
                            fluidRow(div(
                              h5(
                                'Here we filter out genes'
                              ),
                              align = 'center'
                            )),
                            fluidRow(
                              textInput('selectIds', 'regular expression for selection of genes to be removed', value = '^MT-|^RP'))
                            ,fluidRow(
                              h4('GeneList Selection')
                              ,shinyTree("geneListSelection", checkbox = TRUE)
                            ),fluidRow(
                              h4('Minoverall expression'),
                              numericInput('minGenesGS', 'Min # of UMIs over all cells', 1, min=0, max = 1000000)
                            ),
                            br(),
                            fluidRow(
                              h4('Selected genes'),br(),
                              textOutput("gsSelectedGenes", inline = FALSE)
                            ),br(),
                            fluidRow(
                              h4('Genes kept, with mean Expression, and number of cells expressing min 1', align = "center"),
                              br(),
                              DT::dataTableOutput('selectedGenesTable')
                            ),
                            br(),
                            fluidRow(
                              h4('Selected genes'),br(),
                              textOutput("gsrmGenes", inline = FALSE)
                            ),br(),
                            fluidRow(
                              h4('Genes removed, with mean Expression, and number of cells expressing min 1', align = "center"),
                              column(4, offset = 3,
                                     DT::dataTableOutput('removedGenesTable'))
                            )
)

cellSelectionTab  = tabItem(tabName = "cellSelection",
                            fluidRow(div(h3('10X selector'), align = 'center')),
                            br(),
                            fluidRow(div(
                              h5(
                                'Here we filter out cells'
                              ),
                              align = 'center'
                            )),fluidRow(
                              column(
                                2,
                                textInput('minExpGenes', 'List of genes with minimal expression', value = '^CD7$|^KIT$')
                              )
                            ),
                            br(),
                            fluidRow(
                              column(
                                2,
                                numericInput('minGenes', 'Min # of UMIs', 1, min=0, max = 1000000)
                              )),
                            br(),
                            fluidRow(
                              column(
                                2,
                                numericInput('maxGenes', 'Max # of UMIs', 1000000, min=10, max = 1000000)
                              ))
                            
                            
)





