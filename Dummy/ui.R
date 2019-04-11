# menu list
# defines the main entry
menuList <- list(
  menuItem("DummyTools",
    tabName = "dummyMenuList", startExpanded = FALSE,
    # id = "dummyID",
    menuSubItem("dummyTab", tabName = "dummyTab")
  )
)


# tabs with the actual content
tabList <- list(
  crHeatMapTab = tabItem(
    "dummyTab",
    tags$h3("Dummy plot"),
    fluidRow(
      plotOutput("Dummy_plot") # %>% withSpinner()
    ),
    fluidRow(
      imageOutput("DummySavedPlot") # %>% withSpinner()
    )
  )
)

