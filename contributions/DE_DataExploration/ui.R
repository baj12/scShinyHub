
# list of menu Items
menuList <- list(
  menuItem("Data Exploration",
    tabName = "expore", startExpanded = FALSE,
<<<<<<< HEAD:contributions/DE_DataExploration/ui.R
    menuSubItem("Expression", tabName = "DE_expression"),
=======
    menuSubItem("Expression", tabName = "expression"),
>>>>>>> 5086e5a710a9fa9022afdac260eee376be79abfb:contributions/DE_DataExploration/ui.R
    menuSubItem("Panel plot", tabName = "DE_panelPlot"),
    menuSubItem("Scater QC", tabName = "DE_scaterQC")
  )
)

source("modulesUI.R")
# list of tab Items
tabList <- list(
  expressionTab = shinydashboard::tabItem(
    "DE_expression",
    shiny::fluidRow(div(
      htmltools::p(strong("\tInformation:")),
      htmltools::tags$ul(
        tags$li(
          strong("Clustering"),
          ":Clustering was performed with t-SNE followed by identification using DBSCAN"
        ),
        tags$li(
          strong("Cluster 0"),
          ":Cells that cannot be assigned to any cluster"
        ),
        tags$li(
          strong("3D Plot"),
          ":Enter gene name to visualize expression in a single cell"
        ),
        tags$li(
          strong("2D Plot"),
          ":Pick a cluster, highlight cells of interest to download gene expression matrix"
        )
      )
    )),
    br(),
    br(),
    fluidRow(
      column(
        6,
        fluidRow(
          column(
            4,
            textInput("DE_gene_id", "Enter gene", value = defaultValueSingleGene)
          )
        ),
        jqui_resizable(plotlyOutput("DE_tsne_plt"))
      ), column(
        6,
        clusterUI("DE_expclusters")
      )
    ),
    br(),
    fluidRow(column(
      12,
      jqui_resizable( plotOutput("DE_gene_vio_plot") %>% withSpinner())
    ))
  ),

  DE_panelPlotTab = tabItem(
    "DE_panelPlot",
    tags$ul(
      tags$li(
        strong("Panel plot"),
        ":Select a cluster. Enter",
        strong("ONE"),
        "or",
        strong("MULTIPLE"),
        "gene ids to visualize expression in all clusters"
      ),
      tags$li("If the x-axis is a categorical value and the y-axis is UMI.counts the y-axis related to the count for that gene. Otherwise, all genes are used")
    ),
    fluidRow(
      column(
        2,
        uiOutput("DE_clusterSelectionPanelPlot")
      ),
      column(
        2,
        selectInput(
          "DE_dim_x",
          label = "X",
          choice = c("tsne1", "tsne2", "tsne3"),
          selected = "tsne1"
        )
      ),
      column(
        2,
        selectInput(
          "DE_dim_y",
          label = "Y",
          choice = c("tsne1", "tsne2", "tsne3"),
          selected = "tsne2"
        )
      ),
      column(
        2,

        textInput("DE_panelplotids", "Comma seperated gene names", value = defaultValueMultiGenes)
      )
    ),
    fluidRow(column(
      12,
      jqui_resizable(plotOutput("DE_panelPlot") )
    ))
  ),

  DE_scaterQCTab = tabItem(
    "DE_scaterQC",
    tags$ul(
      tags$li(
        strong("Scater QC plots")
      ),
      fluidRow(
        column(
          10,
          offset = 1,
          imageOutput("DE_scaterQC") %>% withSpinner() # PNG output with temp file
        )
      )
    )
  )
)
