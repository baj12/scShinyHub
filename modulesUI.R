# to select clusters from the list of available knn clusters
DEBUG=TRUE

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
             choice = c('V1', 'V2', 'V3'),
             selected = 'V1'
           )),
    column(4,
           selectInput(
             ns('dimension_y'),
             label = 'Y',
             choice = c('V1', 'V2', 'V3'),
             selected = 'V2'
           ))),
    plotOutput(ns('clusterPlot'), brush = brushOpts(id = ns('b1')))
  )
}

