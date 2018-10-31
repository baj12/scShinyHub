
menuList <- list(
  menuItem("CellRangerTools",
    tabName = "cellRanger", startExpanded = FALSE,
    menuSubItem("pheatmap", tabName = "crHeatMap")
  )
)


tabList <- list(
  crHeatMapTab = tabItem(
    "crHeatMap",
    tags$h3("Heatmap plot"),
    fluidRow(column(
      8,
      pHeatMapUI("crHeatmapPlotModule") %>% withSpinner()
      # plotOutput('crHeat_plot1', height='auto', brush = brushOpts(id =
      #                                                "crh1")) %>% withSpinner()
    )),
    tableSelectionUi("crPrioGenesTableModule") %>% withSpinner()
    # column(
    #   2,
    #   uiOutput("clusters5")
    # ),
    # DT::dataTableOutput("crPrioGenes") %>% withSpinner()
  )
)