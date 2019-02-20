# menu list
# defines the main entry 
menuList =  list(
  menuItem("Trajectories", tabName = "TrajectoryList", startExpanded = FALSE,
           menuSubItem("Scorpius", tabName = "scorpiusTab")
           , menuSubItem("ELPIGraph", tabName = "elpiGraphTab")
  )
)


# tabs with the actual content
tabList = list(
  crHeatMapTab = tabItem("scorpiusTab",
                         tags$h3("trajectory in 2D space"),
                         fluidRow(column(12,offset = 1,
                                         checkboxInput("scorpiusCalc", "calculate", FALSE)))
                         ,
                         fluidRow(
                           column(4,
                                  selectInput(
                                    'dimScorpiusX',
                                    label = 'Component 1',
                                    choices = c('tsne1', 'tsne2', 'tsne3'),
                                    selected = 'tsne1'
                                  )),
                           column(4,
                                  selectInput(
                                    'dimScorpiusY',
                                    label = 'Component 2',
                                    choices = c('tsne1', 'tsne2', 'tsne3'),
                                    selected = 'tsne2'
                                  )),
                           column(4,
                                  selectInput(
                                    'dimScorpiusCol',
                                    label = 'Color by',
                                    choices = c('sample', 'tsne1', 'tsne2', 'tsne3'),
                                    selected = 'sample'
                                  ))),
                         fluidRow(column(12,
                                         tipify(
                                           downloadButton("downLoadTraj", "Download trajectory"),
                                           "<h3>download trajectory to csv file</h3>"
                                         ),fileInput(
                                           'trajInputFile',
                                           'Choose .csv file with trajectory informaiton',
                                           accept = c(
                                             '.csv',
                                             'text/comma-separated-values',
                                             'text/tab-separated-values',
                                             'text/plain',
                                             '.csv',
                                             '.tsv'
                                           )
                                         )
                         )  ),
                         fluidRow(column(12,
                                         plotOutput('scropius_trajectory_plot', height = '672px') #%>% withSpinner()
                         )  ),
                         # tags$h3("Heatmap "),
                         fluidRow(column(12,
                                         pHeatMapUI("scorpiusHeatmapPlotModule") %>% withSpinner()
                                         # imageOutput('scorpiusHeatmapPlotModule', height = '672px') #%>% withSpinner() 
                         )),
                         tags$h3("table"),
                         fluidRow(column(
                           10,offset = 1,
                           tableSelectionUi("scorpiusTableMod")
                         ))
  )
  ,
  elpiTab = tabItem("elpiGraphTab",
                    tags$h3("trajectory by ElpiGraph"),
                    fluidRow(column(12,offset = 1,
                                    checkboxInput("elpiCalc", "calculate", FALSE)))
                    ,
                    fluidRow(
                      column(4,
                             selectInput(
                               'dimElpi',
                               label = 'Dimensions to use',
                               choices = c('elpiPCA'),
                               selected = 'elpiPCA'
                             ))
                      ,
                      column(4,
                             selectInput(
                               'ElpiMethod',
                               label = 'Method to use',
                               choices = c('computeElasticPrincipalCurve',
                                           'computeElasticPrincipalTree',
                                           'computeElasticPrincipalCircle'),
                               selected = 'computeElasticPrincipalTree'
                             )),
                      column(4,
                             numericInput(
                               inputId = 'elpiNumNodes',
                               label = 'Number of nodes',
                               value = 60,
                               min = 10, max = 100, step = 1
                             )),
                      column(4,
                             numericInput(
                               inputId = 'elpinReps',
                               label = 'Number of repeats',
                               value = 1,
                               min = 1, max = 50, step = 1
                             )),
                      column(4,
                             numericInput(
                               inputId = 'elpiProbPoint',
                               label = 'probability of inclusing of a single point for each computation',
                               value = 0.6,
                               min = 0.1, max = 1, step = 0.1
                             ))),
                    fluidRow(column(12,
                                    plotOutput('elpi_plot', height = '672px') #%>% withSpinner()
                    ))  ,
                    # tags$h3("Heatmap "),
                    fluidRow(column(12,
                                    plotOutput('elpi_histo', height = '672px') #%>% withSpinner()
                    ))
                    
  )
  
)

# declare heavy calculations
myHeavyCalculations = NULL


