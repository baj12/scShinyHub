
# list of menu Items
menuList =  list(
  menuItem("Data Exploration", tabName = "expore", startExpanded = FALSE,
           menuSubItem("Expression", tabName = "expression"),
           menuSubItem("Panel plot", tabName = "panelPlot"),
           menuSubItem("Scater QC", tabName = "scaterQC")
  )
  
)

source("modulesUI.R")
# list of tab Items
tabList = list(
  expressionTab = tabItem("expression",
                          fluidRow(div(
                            p(strong('\tInformation:')),
                            tags$ul(
                              tags$li(
                                strong('Clustering'),
                                ':Clustering was performed with t-SNE followed by identification using DBSCAN'
                              ),
                              tags$li(
                                strong('Cluster 0'),
                                ':Cells that cannot be assigned to any cluster'
                              ),
                              tags$li(
                                strong('3D Plot'),
                                ':Enter gene name to visualize expression in a single cell'
                              ),
                              tags$li(
                                strong('2D Plot'),
                                ':Pick a cluster, highlight cells of interest to download gene expression matrix'
                              )
                            )
                          )),
                          br(),
                          br(),
                          fluidRow(
                            column(6,
                                   fluidRow(
                                     column(
                                       2,
                                       textInput('gene_id', 'Enter gene', value = 'CD7')
                                     )
                                     ,
                                     
                                     column(
                                       2,
                                       div(
                                         align = "center",
                                         style = "margin-center:50px; margin-top:25px",
                                         downloadButton("downloadExpression", "Download Expression")
                                       )
                                     )
                                   ),
                                   plotlyOutput('tsne_plt'))
                            ,column(6,
                                    clusterUI("expclusters")
                            )
                          ),
                          br(),
                          fluidRow(column(12, 
                                          plotOutput('gene_vio_plot')))
                          
  ),
  
  panelPlotTab = tabItem("panelPlot",
                         tags$ul(
                           tags$li(
                             strong('Panel plot'),
                             ':Select a cluster. Enter',
                             strong('ONE'),
                             'or',
                             strong('MULTIPLE'),
                             'gene ids to visualize expression in all clusters'
                           )
                           
                         ),
                         fluidRow(
                           column(2,
                                  uiOutput("clusters4")),
                           column(
                             2,
                             selectInput(
                               'dimension_x4',
                               label = 'X',
                               choice = c('V1', 'V2', 'V3'),
                               selected = 'V1'
                             )
                           ),
                           column(
                             2,
                             selectInput(
                               'dimension_y4',
                               label = 'Y',
                               choice = c('V1', 'V2', 'V3'),
                               selected = 'V2'
                             )
                           ),
                           column(
                             2,
                             
                             textInput('panelplotids', 'Comma seperated gene names', value = 'CD7')
                           )
                         ),
                         fluidRow(column(
                           12, 
                           plotOutput('panelPlot')
                         ))
                         
  ),
  
  scaterQCTab = tabItem("scaterQC",
                        tags$ul(
                          tags$li(
                            strong('Scater QC plots')
                          ),
                          fluidRow(
                            column(
                              10, offset = 1,
                              plotOutput('scaterQC')
                            )
                          )
                          
                        )
  )
  
)


