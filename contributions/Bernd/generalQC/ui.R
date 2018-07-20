menuList =  list(
  menuItem("General QC", tabName = "generalQC", startExpanded = FALSE,
           menuSubItem("UMI histogram", tabName = "umiHist"),
           menuSubItem("Sample histogram", tabName = "sampleHist"),
           menuSubItem("PC variance", tabName = "variancePC"),
           menuSubItem("TSNE plot", tabName = "tsnePlot")
  )
)

tabList = list(
  tabItem("umiHist",
                        tags$h3("Histogram of UMI counts"),
                        fluidRow(column(
                          10, offset = 1,
                          plotOutput('plotUmiHist') %>% withSpinner()
                        ))
  ),
  
  tabItem("sampleHist",
          tags$h3("Histogram of cells per sample"),
          fluidRow(column(
            10, offset = 1,
            plotOutput('plotSampleHist') %>% withSpinner()
          ))
  ),
  
  tabItem("variancePC",
                          tags$h3("Variance of PCs"),
                          fluidRow(column(
                            10, offset = 1,
                            plotOutput('variancePCA') %>% withSpinner()
                          ))
  ),
  tsnePlotTab = tabItem(tabName = "tsnePlot",
                        fluidRow(div(h3('TSNE Plot'), align = 'center')),
                        br(),
                        fluidRow(column(12,
                                        numericInput("tsneDim","Tsne dimensions", 3, min=3,max=5))),
                        fluidRow(column(12,
                                        numericInput("tsnePerplexity","Tsne tsnePerplexity", 30, min=1,max=100))),
                        fluidRow(column(12,
                                        numericInput("tsneTheta","Tsne tsneTheta", 0.5, min=0.1,max=10))),
                        fluidRow(column(12,
                                        numericInput("tsneSeed","Tsne tsneSeed", 1, min=1,max=10000))),
                        fluidRow(column(12,
                                        plotlyOutput('tsne_main') %>% withSpinner()
                                        )
                        ),
                        fluidRow(column(
                          10,offset = 1,
                          tableSelectionUi("cellSelectionTSNEMod")
                        )
  )
  
)
)

