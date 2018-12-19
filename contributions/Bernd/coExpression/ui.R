
source("modulesUI.R")

# list of menu Items
menuList <- list(
  menuItem("Co-expression",
    tabName = "coexpression", startExpanded = FALSE,
    menuSubItem("All clusters", tabName = "coexpressionAll"),
    menuSubItem("Selected", tabName = "coexpressionSelected"),
    # menuSubItem("binarized", tabName = "coexpressionBinarized"),
    menuSubItem("Co-expression Violin plot", tabName = "CoExpressionViolin"),
    menuSubItem("SOM cluster", tabName = "SOMcluster")
  )
)



# list of tab Items
tabList <- list(
  coexpressionAllTab = tabItem(
    "coexpressionAll",
    fluidRow(
      column(
        7,
        offset = 1,

        textInput("heatmap_geneids", "Comma seperated gene names", value = defaultValueMultiGenes)
      )
    ),
    fluidRow(column(
      10,
      offset = 1,
      pHeatMapUI("coExpHeatmapModule") %>% withSpinner()
    ))
  ),

  coexpressionSelectedTab = tabItem(
    "coexpressionSelected",
    tags$ul(
      tags$li(
        strong("Subclustering"),
        ":Select a group of cells in plot1 based based on a single gene expression. Enter multiple gene ids to assess the co-expression of genes in these cells"
      )
    ),
    fluidRow(
      column(
        2,
        textInput("gene_id_sch", "Enter gene", value = defaultValueSingleGene)
      ),
      column(
        10,
        clusterUI("selected")
      )
    ),
    fluidRow(column(
      6,
      offset = 1,

      textInput("heatmap_geneids2", "Comma seperated gene names", value = defaultValueMultiGenes)
    )),
    fluidRow(column(
      10,
      offset = 1,
      pHeatMapUI("heatmapSelectedModule") %>% withSpinner()
    )), br(),
    fluidRow(
      column(
        3,
        numericInput("coEtgMinExpr", "min UMI count per gene:",
          1,
          min = 0, max = 100000
        )
      ), column(
        3,
        numericInput("coEtgPerc", "min percentage of cells expressing a genes:",
          60,
          min = 1, max = 100
        )
      )
    ),
    tableSelectionUi("topExpGenes") %>% withSpinner()


    # ,
    # fluidRow(column(
    #   10, offset = 1,
    #   plotOutput('selectedHeatmap') %>% withSpinner()
    # ))
  ),
  # binarizeTab = tabItem(
  #   "coexpressionBinarized",
  #   tags$ul(
  #     tags$li(
  #       strong("Binary Expression"),
  #       ":Select a cluster. Enter",
  #       strong("ONE"),
  #       "or",
  #       strong("MULTIPLE"),
  #       "gene ids to assess the co-expression of genes in these cells. Highlighted cells have all genes expressed as determined by a GMM"
  #     )
  #   ),
  # fluidRow(
  #   column(
  #     2,
  #     uiOutput("clusters3")
  #   ),
  #   column(
  #     2,
  #     selectInput(
  #       "dimension_x3",
  #       label = "X",
  #       choices = c("tsne1", "tsne2", "tsne3"),
  #       selected = "tsne1"
  #     )
  #   ),
  #   column(
  #     2,
  #     selectInput(
  #       "dimension_y3",
  #       label = "Y",
  #       # TODO  get from tsne columns
  #       choices = c("tsne1", "tsne2", "tsne3"),
  #       selected = "tsne2"
  #     )
  #   ),
  #   column(
  #     6,
  #
  #     textAreaInput("mclustids", "Comma seperated gene names", value = defaultValueMultiGenes)
  #   )
  # ),
  # fluidRow(column(
  #   10,
  #   offset = 1,
  #   plotOutput("plotCoExpression") %>% withSpinner()
  # )),
  # fluidRow(
  #   div(
  #     align = "center",
  #     style = "margin-center:50px; margin-top:25px",
  #     downloadButton(
  #       "downloadExpressionOnOff",
  #       "Download Expression +ve Cells in cluster"
  #     )
  #   )
  # ),
  #   br(),
  #   br(),
  #   br(),
  #   fluidRow(
  #     h4("Positive Cells in all clusters", align = "center"),
  #     column(6,
  #       offset = 3,
  #       DT::dataTableOutput("onOffTable") %>% withSpinner()
  #     )
  #   )
  # ),

  expressionTab = tabItem(
    "CoExpressionViolin",
    fluidRow(div(
      p(strong("\tInformation:")),
      tags$ul(
        tags$li(
          strong("Violin plot"),
          "for each cell we count how many of the genes specified have an expression larger or equal than the minimum exprssion.\nThese counts are then divided up for any variable that can be used as a factor (has less than 20 levels)."
        )
      )
    )),
    br(),
    br(),
    tipify(
      checkboxInput("showPermutations", "show Permutations", FALSE),
      "check this if you are working on the cell/gene selection to avoid certain calculations"
    ),
    fluidRow(
      column(
        3,

        textInput("geneGrpVioIds", "Comma seperated gene names", value = defaultValueMultiGenes)
      ),
      column(
        3,
        selectInput(
          "dimension_xVioiGrp",
          label = "X",
          choices = c("tsne1", "tsne2", "tsne3"),
          selected = "tsne1"
        )
      ),
      column(
        3,
        numericInput("coEminExpr", "min expression of genes:",
          1,
          min = 1, max = 100000
        )
      )
    ),
    br(),
    fluidRow(column(
      12,
      plotOutput("geneGrp_vio_plot") %>% withSpinner()
    ))
  ),
  tabList = tabItem(
    "SOMcluster",
    fluidRow(div(p(strong("\t Self organizing map and heatmap of cluster with gene of interest")))),
    br(),
    fluidRow(column(
      3,
      numericInput("dimSOM", "number of nodes per dimension",
        20,
        min = 1, max = 100
      ),
      textInput("geneSOM", "Gene of interest", value = defaultValueSingleGene)
    )),
    br(),
    fluidRow(column(
      10,
      offset = 1,
      pHeatMapUI("heatmapSOM") %>% withSpinner()
    )),
    br()
  )
)
