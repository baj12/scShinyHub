# to select clusters from the list of available knn clusters

clusterUI <- function(id){
  if(DEBUG)cat(file=stderr(), paste("clusterUI: ", NS(id)("clusters"), "\n"))
  ns <- NS(id)
  tagList(fluidRow(
    column(4,
           uiOutput(ns("clusters"))),
    column(4,
           selectInput(
             ns('dimension_x'),
             label = 'X',
             choice = c('tsne1', 'tsne2', 'tsne3'),
             selected = 'tsne1'
           )),
    column(4,
           selectInput(
             ns('dimension_y'),
             label = 'Y',
             choice = c('tsne1', 'tsne2', 'tsne3'),
             selected = 'tsne2'
           ))),
    plotOutput(ns('clusterPlot'), brush = brushOpts(id = ns('b1'))) %>% withSpinner(),
    fluidRow(
      checkboxInput(ns("showCells"), "show cell names", FALSE),
      verbatimTextOutput(ns('cellSelection'))
    )
  )
}

tableSelectionUi <- function(id){
  ns <- NS(id)
  tagList(
    fluidRow(div(
      h5("Selected cell names to be copied"), align='left')
    ),
    fluidRow(
      verbatimTextOutput(ns('cellSelection'))
    ),
    fluidRow(
      h4('Cells', offset = 1),
      checkboxInput(ns("selectAll"), "Select all rows", FALSE),br(),
      DTOutput(ns('cellNameTable')) %>% withSpinner()
    )
  )
}