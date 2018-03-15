menuList =  list(
  menuItem("General QC", tabName = "generalQC", startExpanded = FALSE,
           menuSubItem("UMI histogram", tabName = "umiHist"),
           menuSubItem("PC variance", tabName = "variancePC"),
           menuSubItem("TSNE plot", tabName = "tsnePlot")
  )
)

tabList = list(
  tabItem("umiHist",
                        tags$h3("Histogram of UMI counts"),
                        fluidRow(column(
                          10, offset = 1,
                          plotOutput('plotUmiHist')
                        ))
  ),
  
  tabItem("variancePC",
                          tags$h3("Variance of PCs"),
                          fluidRow(column(
                            10, offset = 1,
                            plotOutput('variancePCA')
                          ))
  ),
  tsnePlotTab = tabItem(tabName = "tsnePlot",
                        fluidRow(div(h3('TSNE Plot'), align = 'center')),
                        br(),
                        
                        fluidRow(plotlyOutput('tsne_main'))
  )
  
)

