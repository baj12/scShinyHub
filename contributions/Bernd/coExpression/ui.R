
source("modulesUI.R")

# list of menu Items
menuList =  list(
  menuItem("Co-expression", tabName = "coexpression", startExpanded = FALSE,
           menuSubItem("All clusters", tabName = "coexpressionAll"),
           menuSubItem("Selected", tabName = "coexpressionSelected"),
           menuSubItem("binarized", tabName = "coexpressionBinarized")
  )
)



# list of tab Items
tabList = list(
  coexpressionAllTab = tabItem("coexpressionAll",
                               fluidRow(
                                 column(
                                   2,
                                   
                                   textInput('heatmap_geneids', 'Comma seperated gene names', value = 'CD7,Cd3d,Genes,Umi')
                                 )
                               ),
                               
                               fluidRow(column(
                                 10, offset = 1,
                                 plotOutput('heatmap')
                               ))
  ),
  
  coexpressionSelectedTab = tabItem("coexpressionSelected",
                                    tags$ul(
                                      tags$li(
                                        strong('Subclustering'),
                                        ':Select a group of cells in plot1 based based on a single gene expression. Enter multiple gene ids to assess the co-expression of genes in these cells'
                                      )
                                      
                                    ),
                                    fluidRow(
                                      column(
                                        2,
                                        textInput('gene_id_sch', 'Enter gene', value = 'CD7')
                                      ),
                                      column(10,
                                             clusterUI("selected"))
                                    ),
                                    #   column(2,
                                    #          uiOutput("clusters2")),
                                    #   column(
                                    #     2,
                                    #     selectInput(
                                    #       'dimension_x2',
                                    #       label = 'X',
                                    #       choice = c('V1', 'V2', 'V3'),
                                    #       selected = 'V1'
                                    #     )
                                    #   ),
                                    #   column(
                                    #     2,
                                    #     selectInput(
                                    #       'dimension_y2',
                                    #       label = 'Y',
                                    #       choice = c('V1', 'V2', 'V3'),
                                    #       selected = 'V2'
                                    #     )
                                    #   )
                                    # ),
                                    # fluidRow(column(
                                    #   5, offset = 1,
                                    #   plotOutput('clusterPlot2', brush = brushOpts(id =
                                    #                                                  'scb1'))
                                    # )),
                                    fluidRow(column(
                                      2,
                                      
                                      textInput('heatmap_geneids2', 'Comma seperated gene names', value = 'CD7,Cd3d')
                                    )),
                                    fluidRow(column(
                                      10, offset = 1,
                                      plotOutput('selectedHeatmap')
                                    ))
  ),
  binarizeTab = tabItem("coexpressionBinarized",
                        tags$ul(
                          tags$li(
                            strong('Binary Expression'),
                            ':Select a cluster. Enter',
                            strong('ONE'),
                            'or',
                            strong('MULTIPLE'),
                            'gene ids to assess the co-expression of genes in these cells. Highlighted cells have all genes expressed as determined by a GMM'
                          )
                          
                        ),
                        fluidRow(
                          column(2,
                                 uiOutput("clusters3")),
                          column(
                            2,
                            selectInput(
                              'dimension_x3',
                              label = 'X',
                              choice = c('V1', 'V2', 'V3'),
                              selected = 'V1'
                            )
                          ),
                          column(
                            2,
                            selectInput(
                              'dimension_y3',
                              label = 'Y',
                              choice = c('V1', 'V2', 'V3'),
                              selected = 'V2'
                            )
                          ),
                          column(
                            2,
                            
                            textInput('mclustids', 'Comma seperated gene names', value = 'CD7')
                          )
                        ),
                        fluidRow(column(
                          10, offset = 1,
                          plotOutput('plotCoExpression')
                        )),
                        fluidRow(
                          div(
                            align = "center",
                            style = "margin-center:50px; margin-top:25px",
                            downloadButton(
                              "downloadExpressionOnOff",
                              "Download Expression +ve Cells in cluster"
                            )
                            
                          )
                        ),
                        br(),
                        br(),
                        br(),
                        fluidRow(
                          h4('Positive Cells in all clusters', align = "center"),
                          column(6, offset = 3,
                                 DT::dataTableOutput('onOffTable'))
                        )
  )
  
  
)


