
 # list of menu Items
menuList =  list()

# E.g.
# menuList =  list(
#   menuItem("CellRangerTools", tabName = "cellRanger", startExpanded = FALSE,
#            menuSubItem("pheatmap", tabName = "crHeatMap")
#   )
# )


# list of tab Items
tabList = list()


# E.g.
# tabList = list(
#   crHeatMapTab = tabItem("crHeatMap",
#                          tags$h3("Heatmap plot"),
#                          fluidRow(
#                            plotOutput('crHeat_plot1', brush = brushOpts(id =
#                                                                           "crh1"))
#                          ),
#                          column(2,
#                                 uiOutput("clusters5")),
#                          DT::dataTableOutput('crPrioGenes')
#   ),
#   pheatmapTab = tabItem("cellRanger",
#                         tags$h3("CellRanger tools"),
#                         fluidRow(column(
#                           10, offset = 1
#                           # ,
#                           # plotOutput('variancePCA')
#                         ))
#   )
# )
