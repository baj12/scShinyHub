
source("modulesUI.R")

# list of menu Items
menuList =  list(
  menuItem("Co-expression", tabName = "coexpression", startExpanded = FALSE,
           menuSubItem("All clusters", tabName = "coexpressionAll"),
           menuSubItem("Selected", tabName = "coexpressionSelected"),
           menuSubItem("binarized", tabName = "coexpressionBinarized"),
           menuSubItem("Co-expression Violin plot", tabName = "CoExpressionViolin")
  )
)



# list of tab Items
tabList = list(
  coexpressionAllTab = tabItem("coexpressionAll",
                               fluidRow(
                                 column(
                                   7, offset = 1,
                                   
                                   textInput('heatmap_geneids', 'Comma seperated gene names', value = defaultValueMultiGenes)
                                 )
                               ),
                               
                               fluidRow(column(
                                 10, offset = 1,
                                 plotOutput('heatmap') %>% withSpinner()
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
                                        textInput('gene_id_sch', 'Enter gene', value = defaultValueSingleGene)
                                      ),
                                      column(10,
                                             clusterUI("selected"))
                                    ),
                                    fluidRow(column(
                                      6, offset = 1,
                                      
                                      textInput('heatmap_geneids2', 'Comma seperated gene names', value = defaultValueMultiGenes)
                                    )),
                                    fluidRow(column(
                                      10, offset = 1,
                                      plotOutput('selectedHeatmap') %>% withSpinner()
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
                              choice = c('tsne1', 'tsne2', 'tsne3'),
                              selected = 'tsne1'
                            )
                          ),
                          column(
                            2,
                            selectInput(
                              'dimension_y3',
                              label = 'Y',
                              # TODO  get from tsne columns
                              choice = c('tsne1', 'tsne2', 'tsne3'),
                              selected = 'tsne2'
                            )
                          ),
                          column(
                            6,
                            
                            textAreaInput('mclustids', 'Comma seperated gene names', value = defaultValueMultiGenes)
                          )
                        ),
                        fluidRow(column(
                          10, offset = 1,
                          plotOutput('plotCoExpression') %>% withSpinner()
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
                                 DT::dataTableOutput('onOffTable') %>% withSpinner()
                          )
                        )
  ),
  tabList = list(
    expressionTab = tabItem("CoExpressionViolin",
                            fluidRow(div(
                              p(strong('\tInformation:')),
                              tags$ul(
                                tags$li(
                                  strong('Violin plot'),
                                  'for each cell we count how many of the genes specified have an expression larger or equal than the minimum exprssion.\nThese counts are then divided up for any variable that can be used as a factor (has less than 20 levels).'
                                )
                              )
                            )),
                            br(),
                            br(),
                            
                            fluidRow(
                              column(
                                3,
                                
                                textInput('geneGrpVioIds', 'Comma seperated gene names', value = defaultValueMultiGenes)
                              ),    
                              column(3,
                                     selectInput(
                                       'dimension_xVioiGrp',
                                       label = 'X',
                                       choice = c('tsne1', 'tsne2', 'tsne3'),
                                       selected = 'tsne1'
                                     )
                              ),
                              column(3,
                                     numericInput("coEminExpr", "min expression of genes:",
                                                  10, min = 1, max = 100000)
                              )
                            ),
                            br(),
                            fluidRow(column(12, 
                                            plotOutput('geneGrp_vio_plot') %>% withSpinner()
                            ))
                            
    )
  )
  
  
  
)


