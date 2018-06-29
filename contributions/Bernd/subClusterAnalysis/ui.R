menuList =   list(    
  menuItem("Subcluster analysis", tabName = "subcluster", startExpanded = FALSE,
           menuSubItem("DGE analysis", tabName = "dge")
  )
)

tabList = list(
  dgeTab = tabItem("dge",
                   tags$ul(
                     tags$li(
                       strong('Subclustering'),
                       ':Select a group of cells in plot1 and a different group of cells in plot2 for identifying differential features between these subclusters'
                     )
                     
                   ),
                   
                   fluidRow(
                     column(2,
                            uiOutput("clusters1") %>% withSpinner()
                            ),
                     column(
                       2,
                       selectInput(
                         'dimension_x1',
                         label = 'X',
                         choice = c('tsne1', 'tsne2', 'tsne3'),
                         selected = 'tsne1'
                       )
                     ),
                     column(
                       2,
                       selectInput(
                         'dimension_y1',
                         label = 'Y',
                         choice = c('tsne1', 'tsne2', 'tsne3'),
                         selected = 'tsne2'
                       )
                     )
                   ),
                   
                   fluidRow(
                     column(6,
                            plotOutput('dge_plot1', brush = brushOpts(id =
                                                                        "db1")) %>% withSpinner()
                            ),
                     column(6,
                            plotOutput('dge_plot2', brush = brushOpts(id =
                                                                        'db2')) %>% withSpinner()
                            )
                   ),
                   fluidRow(column(11, offset = 1,
                     h4('Selected genes'),br(),
                     textOutput("crSelectedGenes", inline = FALSE)
                     
                   )),br(),
                   fluidRow(column(11, offset = 1,
                     h4('Top Differentially Expressed Genes', offset = 1),
                     DT::dataTableOutput('dge') %>% withSpinner()
                   )),br(),
                   fluidRow(
                     div(
                       align = "right",
                       style = "margin-right:15px; margin-bottom:10px",
                       downloadButton("download_dge_table", "Download DGE Table")
                     )
                     
                   )
  )
)