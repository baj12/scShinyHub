# menu list
# defines the main entry 
menuList =  list(
  menuItem("Trajectories", tabName = "TrajectoryList", startExpanded = FALSE,
           menuSubItem("Scorpius", tabName = "scorpiusTab")
  )
)


# tabs with the actual content
tabList = list(
  crHeatMapTab = tabItem("scorpiusTab",
                         tags$h3("trajectory in 2D space"),
                         fluidRow(12,checkboxInput("scorpiusCalc", "calculate", FALSE))
                         ,
                         fluidRow(12,
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
  
                         fluidRow(
                           plotOutput('scropius_trajectory_plot') #%>% withSpinner()
                         )  ,
                         tags$h3("Heatmap "),
                         fluidRow(
                           imageOutput('scorpiusHeatmapPlot') #%>% withSpinner() 
                         )
  )
)

# declare heavy calculations
myHeavyCalculations = NULL


