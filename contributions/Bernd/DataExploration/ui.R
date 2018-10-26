
# list of menu Items
menuList =  list(
  menuItem("Data Exploration", tabName = "expore", startExpanded = FALSE,
           menuSubItem("Expression", tabName = "expression"),
           menuSubItem("Panel plot", tabName = "panelPlot"),
           menuSubItem("Coeff. Variance", tabName = "coefVar"),
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
                                       4,
                                       textInput('gene_id', 'Enter gene', value = defaultValueSingleGene)
                                     )
                                     ,
                                     
                                     column(
                                       2,
                                       div(
                                         align = "right",
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
                                          plotOutput('gene_vio_plot') %>% withSpinner()
                          ))
                          
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
                               choice = c('tsne1', 'tsne2', 'tsne3'),
                               selected = 'tsne1'
                             )
                           ),
                           column(
                             2,
                             selectInput(
                               'dimension_y4',
                               label = 'Y',
                               choice = c('tsne1', 'tsne2', 'tsne3'),
                               selected = 'tsne2'
                             )
                           ),
                           column(
                             2,
                             
                             textInput('panelplotids', 'Comma seperated gene names', value = defaultValueMultiGenes)
                           )
                         ),
                         fluidRow(column(
                           12, 
                           plotOutput('panelPlot') %>% withSpinner()
                         ))
                         
  ),
  coefVarTab = tabItem("coefVar",
                       tags$h3("Coefficient of Variance"),
                       tags$p("plot genes based coefficient of variance. If no gene is given, plot the the first 50 genes ordered by increasing CV."),
                       tags$p("One gene given, plot heatmap with most similar genes."),
                       tags$p("More than one gene given, plot heatmap of these genes."),
                       tags$p("All of the abobe are only clustered by cells, based on the here visible genes."),
                       tags$h4("Histogram of coefficients"),
                       fluidRow(column(12,
                                       plotOutput('cvHist') %>% withSpinner()
                       )),
                       tags$h4("Heatmap of sorted coeficients"),
                       fluidRow(column(
                         6, offset = 1,
                         textInput('cvHeatmap_geneids', 'Comma seperated gene names(none = first 100)', 
                                   value = '')
                       )),
                       fluidRow(column(12,
                                       plotOutput('cvHeatMap') %>% withSpinner()
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
                              imageOutput('scaterQC') %>% withSpinner() # PNG output with temp file
                            )
                          )
                          
                        )
  )
  
)


