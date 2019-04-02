
# list of menu Items
menuList <- list(
  menuItem("Data Exploration",
    tabName = "expore", startExpanded = FALSE,
    menuSubItem("Expression", tabName = "expression"),
    menuSubItem("Panel plot", tabName = "panelPlot"),
    menuSubItem("Scater QC", tabName = "scaterQC")
  )
)

source("modulesUI.R")
# list of tab Items
tabList <- list(
  expressionTab = shinydashboard::tabItem(
    "expression",
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
            textInput("gene_id", "Enter gene", value = defaultValueSingleGene)
          ),

          column(
            2,
            div(
              align = "right",
              style = "margin-center:50px; margin-top:25px",
              downloadButton("downloadExpression", "Download Expression")
            )
          )
        ),
        jqui_resizable(plotlyOutput("tsne_plt"))
      ), column(
        6,
        clusterUI("expclusters")
      )
    ),
    br(),
    fluidRow(column(
      12,
      jqui_resizable( plotOutput("gene_vio_plot") %>% withSpinner())
    ))
  ),

  panelPlotTab = tabItem(
    "panelPlot",
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
        uiOutput("clusterSelectionPanelPlot")
      ),
      column(
        2,
        selectInput(
          "dimension_x4",
          label = "X",
          choice = c("tsne1", "tsne2", "tsne3"),
          selected = "tsne1"
        )
      ),
      column(
        2,
        selectInput(
          "dimension_y4",
          label = "Y",
          choice = c("tsne1", "tsne2", "tsne3"),
          selected = "tsne2"
        )
      ),
      column(
        2,

        textInput("panelplotids", "Comma seperated gene names", value = defaultValueMultiGenes)
      )
    ),
    fluidRow(column(
      12,
      jqui_resizable(plotOutput("panelPlot") )
    ))
  ),

  scaterQCTab = tabItem(
    "scaterQC",
    tags$ul(
      tags$li(
        strong("Scater QC plots")
      ),
      fluidRow(
        column(
          10,
          offset = 1,
          imageOutput("scaterQC") %>% withSpinner() # PNG output with temp file
        )
      )
    )
  )
)
